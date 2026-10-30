"""rtcurate コマンドのエントリポイント。"""

import argparse
import os
import sys
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path

from rtcurate.io import load_spots
from rtcurate.output import summary, write_excluded, write_main
from rtcurate.pair_eval import run
from rtcurate.preprocess import FilterOptions, preprocess
from rtcurate.tags import DEFAULT_RT_TAG, TagFile, merge_tags, read_tags, render_tags, save_with_backup


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="rtcurate",
        description="MS-DIAL5 の AlignmentResult (.arf2) で、RT 順序が矛盾する脂質アノテーションに"
                    "MS-DIAL のタグを付ける。",
    )
    p.add_argument("input", type=Path, help="MS-DIAL5 の AlignmentResult (.arf2)")
    p.add_argument("-o", "--out-dir", type=Path,
                   help="出力先フォルダ(既定: 入力と同じフォルダ。MS-DIAL はここの <入力名>_tags.xml を読む)")
    p.add_argument("--tsv", action="store_true", help="判定結果の表 <入力名>_rtcurated.tsv も出力する")
    p.add_argument("--write-excluded", action="store_true",
                   help="RT 順序で自動除外したものの一覧 <入力名>_rt_excluded.tsv も出力する")
    p.add_argument("--rt-tag", type=int, choices=[2, 3, 4, 5], default=DEFAULT_RT_TAG,
                   help="RT 順序で外れたスポットに付けるタグ ID(既定 5 = Overannotation)")
    p.add_argument("--tag-confirmed", action="store_true",
                   help="ペア評価で残ったスポットに Confirmed (1) を付け足す(既定では 1 に触れない)")
    p.add_argument("--n-miss-threshold", type=int, default=4, help="除外候補になる n_miss の下限(既定 4)")
    p.add_argument("--miss-percent-threshold", type=float, default=10.0,
                   help="除外候補になる miss%% の下限(既定 10)")
    p.add_argument("--keep-low-score", action="store_true", help="low score を除外しない")
    p.add_argument("--keep-no-ms2", action="store_true", help="no MS2 を除外しない")
    p.add_argument("--keep-wo-ms2", action="store_true", help="w/o MS2 を除外しない")
    p.add_argument("--keep-riken", action="store_true", help="RIKEN を除外しない")
    p.add_argument("--keep-ontology-null", action="store_true", help="Ontology なしを除外しない")
    p.add_argument("--no-oxygen-tag", action="store_true", help="Ontology に酸素数タグ (;O2 など) を付けない")
    p.add_argument("--p-tag", action="store_true", help='名前に "p-" があれば Ontology の先頭に付ける')
    p.add_argument("--o-tag", action="store_true", help='名前に "O-" があれば Ontology の末尾に付ける')
    p.add_argument("--exclude-ontology", action="append", default=[], metavar="NAME",
                   help="ペア評価から外す Ontology(複数回指定可)")
    return p


def _utf8_for_git_bash() -> None:
    """Git Bash は出力を UTF-8 として表示するので、Windows 既定の cp932 のままだと文字化けする。"""
    if os.environ.get("MSYSTEM"):
        for stream in (sys.stdout, sys.stderr):
            if hasattr(stream, "reconfigure"):
                stream.reconfigure(encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    _utf8_for_git_bash()
    args = _parser().parse_args(argv)
    src: Path = args.input

    if src.suffix.lower() != ".arf2":
        print(f"エラー: .arf2 ファイルを指定してください: {src}", file=sys.stderr)
        return 2
    if not src.is_file():
        print(f"エラー: ファイルが見つかりません: {src}", file=sys.stderr)
        return 2

    out_dir: Path = args.out_dir or src.parent
    tags_path = out_dir / f"{src.stem}_tags.xml"
    try:
        existing = read_tags(tags_path)
    except (ET.ParseError, ValueError, TypeError) as e:
        print(f"エラー: 既存のタグファイル {tags_path} を読めませんでした(何も書き込んでいません): {e}",
              file=sys.stderr)
        return 2

    try:
        spots = load_spots(src)
    except Exception as e:  # 壊れたファイルなど、読み込みの失敗はすべてここで報告する
        print(f"エラー: {src} を読み込めませんでした: {e}", file=sys.stderr)
        return 2

    opts = FilterOptions(
        low_score=not args.keep_low_score,
        no_ms2=not args.keep_no_ms2,
        wo_ms2=not args.keep_wo_ms2,
        riken=not args.keep_riken,
        ontology_null=not args.keep_ontology_null,
        oxygen_tag=not args.no_oxygen_tag,
        p_tag=args.p_tag,
        o_tag=args.o_tag,
    )
    result = run(preprocess(spots, opts), args.n_miss_threshold, args.miss_percent_threshold,
                 args.exclude_ontology)
    print(summary(result))

    if args.tsv:
        tsv = out_dir / f"{src.stem}_rtcurated.tsv"
        if not _try_write(write_main, result, tsv):
            return 2
        print(f"表: {tsv}")
    if args.write_excluded:
        excluded = out_dir / f"{src.stem}_rt_excluded.tsv"
        if not _try_write(write_excluded, result, excluded):
            return 2
        print(f"除外一覧: {excluded}")

    peaks, both = merge_tags(existing.peaks, result, args.rt_tag, confirm=args.tag_confirmed)
    content = render_tags(TagFile(existing.definitions, peaks))
    try:
        backup = save_with_backup(tags_path, content, datetime.now())
    except OSError as e:
        print(f"エラー: {tags_path} に書き込めませんでした: {e}", file=sys.stderr)
        return 2
    if backup is not None:
        print(f"既存のタグファイルを退避しました: {backup}")
    print(f"タグ: {tags_path}")
    if both:
        print(f"注意: Confirmed (1) とタグ {args.rt_tag} の両方が付いたスポットが {len(both)} 件あります"
              f"(例: ID {', '.join(map(str, both[:10]))})")
    return 0


def _try_write(writer, df, path: Path) -> bool:
    try:
        writer(df, path)
    except OSError as e:
        print(f"エラー: {path} に書き込めませんでした(Excel などで開いていませんか): {e}",
              file=sys.stderr)
        return False
    return True
