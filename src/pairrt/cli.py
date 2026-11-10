"""pairrt コマンドのエントリポイント。"""

import argparse
import os
import sys
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path

import pandas as pd

from pairrt.io import load_spots
from pairrt.output import summary, write_excluded, write_main
from pairrt.pair_eval import run
from pairrt.preprocess import FilterOptions, preprocess
from pairrt.table import (
    GenericColumns, RawTable, TableError, apply_use_column, check_added_columns, generic_frame, msdial_spots,
    read_table, use_mask, write_table,
)
from pairrt.tags import DEFAULT_RT_TAG, TagFile, merge_tags, read_tags, render_tags, save_with_backup

# 汎用 TSV の列名オプション(argparse の属性名 → オプション名)
_COLUMN_OPTIONS = {"rt_column": "--rt-column", "c_column": "--c-column", "db_column": "--db-column",
                   "ontology_column": "--ontology-column"}
# 名前を使う処理のオプション。汎用 TSV では効果がない(--keep-ontology-null は汎用 TSV でも効く)
_NAME_OPTIONS = {"keep_low_score": "--keep-low-score", "keep_no_ms2": "--keep-no-ms2",
                 "keep_wo_ms2": "--keep-wo-ms2", "keep_riken": "--keep-riken",
                 "no_oxygen_tag": "--no-oxygen-tag", "no_ether_split": "--no-ether-split"}
_TAG_OPTIONS = {"rt_tag": "--rt-tag", "tag_confirmed": "--tag-confirmed"}


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="pairrt",
        description="同じ脂質クラスの中で RT 順序が矛盾するアノテーションを見つける。"
                    ".arf2 では MS-DIAL のタグを付け、MS-DIAL の txt と汎用 TSV では判定結果の列を足した表を出す。",
    )
    p.add_argument("input", type=Path,
                   help="MS-DIAL5 の AlignmentResult (.arf2)、MS-DIAL のアライメント結果の txt、"
                        "または RT・C・DB・Ontology の列を持つ表 (.txt / .tsv はタブ区切り、.csv はカンマ区切り)")
    p.add_argument("-o", "--out-dir", type=Path,
                   help="出力先フォルダ(既定: 入力と同じフォルダ。MS-DIAL はここの <入力名>_tags.xml を読む)")
    p.add_argument("--tsv", action="store_true",
                   help="判定結果の表 <入力名>_rtcurated.tsv も出力する(txt / TSV の入力では常に出す)")
    p.add_argument("--write-excluded", action="store_true",
                   help="RT 順序で自動除外したものの一覧 <入力名>_rt_excluded.tsv も出力する")
    p.add_argument("--overwrite", action="store_true",
                   help="出力の表 (_rtcurated.tsv / _rt_excluded.tsv) がすでにあれば上書きする"
                        "(既定ではエラーで止める。_tags.xml は常に古いものを退避してから書く)")
    p.add_argument("--rt-tag", type=int, choices=[2, 3, 4, 5], default=None,
                   help=f"RT 順序で外れたスポットに付けるタグ ID(既定 {DEFAULT_RT_TAG} = Overannotation。.arf2 のみ)")
    p.add_argument("--tag-confirmed", action="store_true",
                   help="ペア評価で残ったスポットに Confirmed (1) を付け足す(既定では 1 に触れない。.arf2 のみ)")
    p.add_argument("--n-miss-threshold", type=int, default=4, help="除外候補になる n_miss の下限(既定 4)")
    p.add_argument("--miss-percent-threshold", type=float, default=10.0,
                   help="除外候補になる miss%% の下限(既定 10)")
    p.add_argument("--keep-low-score", action="store_true", help="low score を除外しない")
    p.add_argument("--keep-no-ms2", action="store_true", help="no MS2 を除外しない")
    p.add_argument("--keep-wo-ms2", action="store_true", help="w/o MS2 を除外しない")
    p.add_argument("--keep-riken", action="store_true", help="RIKEN を除外しない")
    p.add_argument("--keep-ontology-null", action="store_true", help="Ontology なしを除外しない")
    p.add_argument("--no-oxygen-tag", action="store_true", help="Ontology に酸素数タグ (;O2 など) を付けない")
    p.add_argument("--no-ether-split", action="store_true",
                   help="Ether 系の O- と P- を分けずに同じ Ontology で判定する(既定では EtherPE_O / EtherPE_P に分ける)")
    p.add_argument("--exclude-ontology", action="append", default=[], metavar="NAME",
                   help="ペア評価から外す Ontology(複数回指定可)")
    defaults = GenericColumns()
    for option, default in (("--rt-column", defaults.rt), ("--c-column", defaults.c),
                            ("--db-column", defaults.db), ("--ontology-column", defaults.ontology)):
        p.add_argument(option, metavar="NAME", help=f"汎用 TSV の列名(既定 {default}。大文字小文字は区別しない)")
    p.add_argument("--use-column", metavar="NAME",
                   help="この列が FALSE / 0 / no の行を判定から完全に外す(txt / 汎用 TSV のみ)")
    p.add_argument("--header-row", type=int, metavar="N",
                   help="列名の行をファイルの N 行目(1 始まり)に指定する(txt / TSV / CSV のみ。既定: 自動で探す)")
    return p


def _utf8_for_git_bash() -> None:
    """Git Bash は出力を UTF-8 として表示するので、Windows 既定の cp932 のままだと文字化けする。"""
    if os.environ.get("MSYSTEM"):
        for stream in (sys.stdout, sys.stderr):
            if hasattr(stream, "reconfigure"):
                stream.reconfigure(encoding="utf-8")


def _given(args, options: dict[str, str]) -> list[str]:
    """指定された(既定値のままでない)オプションの名前を返す。"""
    return [name for attr, name in options.items() if getattr(args, attr) not in (None, False)]


def _error(message: str) -> int:
    print(f"エラー: {message}", file=sys.stderr)
    return 2


def _refuse_overwrite(args, outputs: list[Path]) -> bool:
    """出力の表がすでにあれば、--overwrite が無い限り何も書かずに止める(_tags.xml は退避するので対象外)。"""
    existing = [p for p in outputs if p.exists()]
    if existing and not args.overwrite:
        names = ", ".join(str(p) for p in existing)
        _error(f"出力先にすでにファイルがあります: {names}(上書きするなら --overwrite を付けてください。"
               "何も書き込んでいません)")
        return True
    return False


def main(argv: list[str] | None = None) -> int:
    _utf8_for_git_bash()
    args = _parser().parse_args(argv)
    src: Path = args.input

    suffix = src.suffix.lower()
    if suffix not in (".arf2", ".txt", ".tsv", ".csv"):
        return _error(f".arf2 / .txt / .tsv / .csv のファイルを指定してください: {src}")
    if not src.is_file():
        return _error(f"ファイルが見つかりません: {src}")

    opts = FilterOptions(
        low_score=not args.keep_low_score,
        no_ms2=not args.keep_no_ms2,
        wo_ms2=not args.keep_wo_ms2,
        riken=not args.keep_riken,
        ontology_null=not args.keep_ontology_null,
        oxygen_tag=not args.no_oxygen_tag,
        ether_split=not args.no_ether_split,
    )
    out_dir: Path = args.out_dir or src.parent
    if suffix == ".arf2":
        return _main_arf2(args, src, out_dir, opts)
    return _main_table(args, src, out_dir, opts)


def _main_arf2(args, src: Path, out_dir: Path, opts: FilterOptions) -> int:
    wrong = _given(args, {**_COLUMN_OPTIONS, "use_column": "--use-column", "header_row": "--header-row"})
    if wrong:
        return _error(f"{', '.join(wrong)} は txt / TSV / CSV の入力でだけ使えます")
    rt_tag = DEFAULT_RT_TAG if args.rt_tag is None else args.rt_tag
    tsv = out_dir / f"{src.stem}_rtcurated.tsv"
    excluded = out_dir / f"{src.stem}_rt_excluded.tsv"
    if _refuse_overwrite(args, [tsv] * args.tsv + [excluded] * args.write_excluded):
        return 2

    tags_path = out_dir / f"{src.stem}_tags.xml"
    try:
        existing = read_tags(tags_path)
    except (ET.ParseError, ValueError, TypeError) as e:
        return _error(f"既存のタグファイル {tags_path} を読めませんでした(何も書き込んでいません): {e}")

    try:
        spots = load_spots(src)
    except Exception as e:  # 壊れたファイルなど、読み込みの失敗はすべてここで報告する
        return _error(f"{src} を読み込めませんでした: {e}")

    result = run(preprocess(spots, opts), args.n_miss_threshold, args.miss_percent_threshold,
                 args.exclude_ontology)
    print(summary(result))

    if args.tsv:
        if not _try_write(write_main, result, tsv):
            return 2
        print(f"表: {tsv}")
    if args.write_excluded:
        if not _try_write(write_excluded, result, excluded):
            return 2
        print(f"除外一覧: {excluded}")

    peaks, both = merge_tags(existing.peaks, result, rt_tag, confirm=args.tag_confirmed)
    content = render_tags(TagFile(existing.definitions, peaks))
    try:
        backup = save_with_backup(tags_path, content, datetime.now())
    except OSError as e:
        return _error(f"{tags_path} に書き込めませんでした: {e}")
    if backup is not None:
        print(f"既存のタグファイルを退避しました: {backup}")
    print(f"タグ: {tags_path}")
    if both:
        print(f"注意: Confirmed (1) とタグ {rt_tag} の両方が付いたスポットが {len(both)} 件あります"
              f"(例: ID {', '.join(map(str, both[:10]))})")
    return 0


def _main_table(args, src: Path, out_dir: Path, opts: FilterOptions) -> int:
    wrong = _given(args, _TAG_OPTIONS)
    if wrong:
        return _error(f"{', '.join(wrong)} は .arf2 の入力でだけ使えます(txt / TSV からはタグファイルを作りません)")
    tsv = out_dir / f"{src.stem}_rtcurated.tsv"
    excluded = out_dir / f"{src.stem}_rt_excluded.tsv"
    if _refuse_overwrite(args, [tsv] + [excluded] * args.write_excluded):
        return 2
    try:
        table = read_table(src, header_row=args.header_row)
        wrong = _given(args, _COLUMN_OPTIONS) if table.kind == "msdial_txt" else []
        if wrong:
            return _error(f"{', '.join(wrong)} は汎用 TSV でだけ使えます({src} は MS-DIAL の txt です)")
        check_added_columns(table)
        df = _table_frame(args, table, opts)
    except TableError as e:
        return _error(f"{src}: {e}")

    result = run(df, args.n_miss_threshold, args.miss_percent_threshold, args.exclude_ontology)
    print(summary(result, total_label="行数" if table.kind == "generic" else "検出スポット数"))

    if not _try_write(lambda r, p: write_table(table, r, p), result, tsv):
        return 2
    print(f"表: {tsv}")
    if args.write_excluded:
        if not _try_write(lambda r, p: write_table(table, r, p, rt_order_only=True), result, excluded):
            return 2
        print(f"除外一覧: {excluded}")
    return 0


def _table_frame(args, table: RawTable, opts: FilterOptions) -> pd.DataFrame:
    """txt / 汎用 TSV から、pair_eval.run に渡す表を作る。"""
    if table.kind == "msdial_txt":
        spots = msdial_spots(table)
        if spots.attrs["ids_from_rows"]:
            print("注意: Alignment ID が重複しているので、同点のときに外す順番には行番号を使います", file=sys.stderr)
        df = preprocess(spots, opts)
        return df if args.use_column is None else apply_use_column(df, use_mask(table.frame, args.use_column))

    ignored = _given(args, _NAME_OPTIONS)
    if ignored:
        print(f"注意: 汎用 TSV では名前を使わないので、{', '.join(ignored)} は効果がありません", file=sys.stderr)
    defaults = GenericColumns()
    cols = GenericColumns(args.rt_column or defaults.rt, args.c_column or defaults.c,
                          args.db_column or defaults.db, args.ontology_column or defaults.ontology)
    use = None if args.use_column is None else use_mask(table.frame, args.use_column)
    return generic_frame(table, cols, use, keep_ontology_null=args.keep_ontology_null)


def _try_write(writer, df, path: Path) -> bool:
    try:
        writer(df, path)
    except OSError as e:
        print(f"エラー: {path} に書き込めませんでした(Excel などで開いていませんか): {e}",
              file=sys.stderr)
        return False
    return True
