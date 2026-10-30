"""rtcurate コマンドのエントリポイント。"""

import argparse
import sys
from pathlib import Path

from rtcurate.io import load_spots
from rtcurate.output import summary, write_excluded, write_main
from rtcurate.pair_eval import run
from rtcurate.preprocess import FilterOptions, preprocess


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="rtcurate",
        description="MS-DIAL5 の AlignmentResult (.arf2) で、RT 順序が矛盾する脂質アノテーションに印を付ける。",
    )
    p.add_argument("input", type=Path, help="MS-DIAL5 の AlignmentResult (.arf2)")
    p.add_argument("-o", "--output", type=Path,
                   help="出力 TSV(既定: 入力と同じフォルダの <入力名>_rtcurated.tsv)")
    p.add_argument("--write-excluded", action="store_true",
                   help="RT 順序で自動除外したものの一覧 <入力名>_rt_excluded.tsv も出力する")
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


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    src: Path = args.input

    if src.suffix.lower() != ".arf2":
        print(f"エラー: .arf2 ファイルを指定してください: {src}", file=sys.stderr)
        return 2
    if not src.is_file():
        print(f"エラー: ファイルが見つかりません: {src}", file=sys.stderr)
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

    output = args.output or src.with_name(f"{src.stem}_rtcurated.tsv")
    if not _try_write(write_main, result, output):
        return 2
    print(summary(result))
    print(f"出力: {output}")

    if args.write_excluded:
        excluded = output.with_name(f"{src.stem}_rt_excluded.tsv")
        if not _try_write(write_excluded, result, excluded):
            return 2
        print(f"除外一覧: {excluded}")
    return 0


def _try_write(writer, df, path: Path) -> bool:
    try:
        writer(df, path)
    except OSError as e:
        print(f"エラー: {path} に書き込めませんでした(Excel などで開いていませんか): {e}",
              file=sys.stderr)
        return False
    return True


if __name__ == "__main__":
    sys.exit(main())
