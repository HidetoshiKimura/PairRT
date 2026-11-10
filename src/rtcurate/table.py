"""MS-DIAL の txt と汎用 TSV の読み込み・変換・書き出し。

設計: docs/2026-11-09-table-input-design.md。入力の値は出力でそのまま書き戻すので、
すべての列を文字列のまま読み、判定に使う列だけを別に数値へ変換する。
"""

import csv
import os
from dataclasses import dataclass

import pandas as pd

MSDIAL_HEADER_SCAN = 20  # MS-DIAL の txt の列名の行を探す範囲(先頭からの行数)
_MSDIAL_REQUIRED = {"Average Rt(min)", "Metabolite name", "Ontology"}
_ENCODINGS = ("utf-8-sig", "cp932")  # Excel で日本語を含む表を保存すると cp932 になる


class TableError(ValueError):
    """入力の表の誤り。CLI はこれを受けて終了コード 2 にする。"""


@dataclass
class RawTable:
    kind: str  # "msdial_txt" | "generic"
    frame: pd.DataFrame  # 列名は元のまま(重複あり)、値はすべて str


def _read_rows(path: str | os.PathLike) -> list[list[str]]:
    for encoding in _ENCODINGS:
        try:
            with open(path, encoding=encoding, newline="") as f:
                return list(csv.reader(f, delimiter="\t"))
        except UnicodeDecodeError:
            continue
    raise TableError(f"文字コードを読み取れません({' / '.join(_ENCODINGS)} のどちらでもありません)")


def _is_msdial_header(row: list[str]) -> bool:
    return bool(row) and row[0].strip() == "Alignment ID" and _MSDIAL_REQUIRED <= set(row)


def read_table(path: str | os.PathLike) -> RawTable:
    """タブ区切りの表を読み、MS-DIAL の txt か汎用 TSV かを見分けて返す。"""
    rows = _read_rows(path)
    while rows and not any(cell.strip() for cell in rows[-1]):  # ファイル末尾の空行
        rows.pop()
    if not rows:
        raise TableError("表が空です")

    header_at = next((i for i, row in enumerate(rows[:MSDIAL_HEADER_SCAN]) if _is_msdial_header(row)),
                     None)
    kind = "generic" if header_at is None else "msdial_txt"
    header_at = header_at or 0
    header = rows[header_at]

    data = rows[header_at + 1:]
    for offset, row in enumerate(data):
        if len(row) != len(header):
            line = header_at + offset + 2  # ファイル上の行番号(1 始まり)
            raise TableError(f"{line} 行目の列数({len(row)})が列名の行の列数({len(header)})と違います")
    return RawTable(kind, pd.DataFrame(data, columns=header, dtype=object))


def find_column(frame: pd.DataFrame, name: str) -> int:
    """大文字小文字を区別せずに列を探し、列の位置(0 始まり)を返す。"""
    hits = [i for i, col in enumerate(frame.columns) if str(col).lower() == name.lower()]
    if not hits:
        raise TableError(f"列 '{name}' がありません")
    if len(hits) > 1:
        found = ", ".join(f"'{frame.columns[i]}'" for i in hits)
        raise TableError(f"列 '{name}' に当たる列が複数あります: {found}")
    return hits[0]
