"""MS-DIAL の txt と汎用 TSV の読み込み・変換・書き出し。

設計: docs/2026-11-09-table-input-design.md。入力の値は出力でそのまま書き戻すので、
すべての列を文字列のまま読み、判定に使う列だけを別に数値へ変換する。
"""

import csv
import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from rtcurate.preprocess import is_null_ontology

MSDIAL_HEADER_SCAN = 20  # MS-DIAL の txt の列名の行を探す範囲(先頭からの行数)
_MSDIAL_REQUIRED = {"Average Rt(min)", "Metabolite name", "Ontology"}
_ENCODINGS = ("utf-8-sig", "cp932")  # Excel で日本語を含む表を保存すると cp932 になる
_CSV_FIELD_LIMIT = 2**31 - 1


class TableError(ValueError):
    """入力の表の誤り。CLI はこれを受けて終了コード 2 にする。"""


@dataclass
class RawTable:
    kind: str  # "msdial_txt" | "generic"
    frame: pd.DataFrame  # 列名は元のまま(重複あり)、値はすべて str


def _split_rows(f, csv_format: bool) -> list[list[str]]:
    if not csv_format:
        # タブと改行だけで区切る。MS-DIAL は値を引用符で囲まないので、" も値の一部としてそのまま読む
        # (csv モジュールの引用符の解釈を使うと、" で始まる値から次の " までの行がまとめられて消える)
        return [line.split("\t") for line in f.read().split("\n")]
    # CSV は値を " で囲む書き方が一般的なので、囲まれた中のカンマは値の一部として読む
    old_limit = csv.field_size_limit(_CSV_FIELD_LIMIT)  # MS/MS スペクトルの列は長いことがある
    try:
        return list(csv.reader(f))
    except csv.Error as e:
        raise TableError(f"CSV として読めません: {e}") from e
    finally:
        csv.field_size_limit(old_limit)


def _read_rows(path: str | os.PathLike, csv_format: bool) -> list[list[str]]:
    with open(path, "rb") as f:
        if f.read(2) in (b"\xff\xfe", b"\xfe\xff"):
            raise TableError("UTF-16 のファイルには対応していません(Excel なら「テキスト(タブ区切り)」"
                             "か UTF-8 で保存し直してください)")
    for encoding in _ENCODINGS:
        try:
            with open(path, encoding=encoding, newline="" if csv_format else None) as f:
                return _split_rows(f, csv_format)
        except UnicodeDecodeError:
            continue
    raise TableError(f"文字コードを読み取れません({' / '.join(_ENCODINGS)} のどちらでもありません)")


def _is_msdial_header(row: list[str]) -> bool:
    return bool(row) and row[0].strip() == "Alignment ID" and _MSDIAL_REQUIRED <= set(row)


def read_table(path: str | os.PathLike, header_row: int | None = None) -> RawTable:
    """表を読み、MS-DIAL の txt か汎用 TSV かを見分けて返す。`.csv` はカンマ区切り、それ以外はタブ区切り。

    header_row(1 始まり)を渡すと、その行を列名の行にする。渡さなければ自動で探す。
    """
    rows = _read_rows(path, csv_format=Path(path).suffix.lower() == ".csv")
    while rows and not any(cell.strip() for cell in rows[-1]):  # ファイル末尾の空行
        rows.pop()
    if not rows:
        raise TableError("表が空です")

    if header_row is not None:
        if not 1 <= header_row <= len(rows):
            raise TableError(f"--header-row {header_row} は範囲外です(1〜{len(rows)} 行目で指定してください)")
        header_at = header_row - 1
        kind = "msdial_txt" if _is_msdial_header(rows[header_at]) else "generic"
    else:
        found = next((i for i, row in enumerate(rows[:MSDIAL_HEADER_SCAN]) if _is_msdial_header(row)), None)
        kind = "generic" if found is None else "msdial_txt"
        header_at = found or 0
    header = rows[header_at]

    data = rows[header_at + 1:]
    for offset, row in enumerate(data):
        for j, cell in enumerate(row):  # 出力はタブ区切りなので、タブや改行を含む値は書けない(CSV でだけ起こる)
            if "\t" in cell or "\n" in cell or "\r" in cell:
                column = header[j] if j < len(header) else f"{j + 1} 列目"
                raise TableError(f"列 '{column}' の {offset + 1} 件目の値にタブか改行が含まれています")
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


# MS-DIAL の txt の列 → load_spots と同じ列名
_MSDIAL_COLUMNS = {
    "Alignment ID": "MasterAlignmentID",
    "Metabolite name": "name",
    "Average Mz": "mz",
    "Average Rt(min)": "rt",
    "Ontology": "Ontology",
    "Adduct type": "Adduct_type",
}
_TRUE = {"true", "1", "yes"}
_FALSE = {"false", "0", "no"}
_SHOWN = 5  # エラーに出す値の件数


def _exact_column(frame: pd.DataFrame, name: str) -> pd.Series:
    hits = [i for i, col in enumerate(frame.columns) if col == name]
    if not hits:
        raise TableError(f"列 '{name}' がありません")
    if len(hits) > 1:
        raise TableError(f"列 '{name}' が複数あります")
    return frame.iloc[:, hits[0]].astype(str).str.strip()


def _bad_values(column: str, values: pd.Series, bad: pd.Series, hint: str) -> TableError:
    """bad が True の値を、データの何件目か(1 始まり)と一緒に最初の数件だけ並べたエラーを作る。"""
    shown = [f"{i + 1} 件目 '{v}'" for i, v in list(values[bad].items())[:_SHOWN]]
    return TableError(f"列 '{column}' に{hint}: {', '.join(shown)}(全 {int(bad.sum())} 件)")


def _to_float(values: pd.Series, column: str) -> pd.Series:
    """null と空は NaN、それ以外で数値にできない値はエラー。"""
    missing = values.isin(["", "null"])
    numbers = pd.to_numeric(values.where(~missing), errors="coerce")
    bad = numbers.isna() & ~missing
    if bad.any():
        raise _bad_values(column, values, bad, "数値として読めない値があります")
    return numbers.astype(float)


def msdial_spots(table: RawTable) -> pd.DataFrame:
    """MS-DIAL の txt から、load_spots と同じ形の表(1 行 = 1 スポット)を作る。"""
    cols = {new: _exact_column(table.frame, old) for old, new in _MSDIAL_COLUMNS.items()}

    ids = pd.to_numeric(cols["MasterAlignmentID"], errors="coerce")
    bad = ids.isna() | (ids != ids.round())
    if bad.any():
        raise _bad_values("Alignment ID", cols["MasterAlignmentID"], bad, "整数でない値があります")
    ids = ids.astype(int)
    # Pos と Neg をまとめた表などで重複していたら、行番号(0 始まり)を ID の代わりに使う。
    # txt からはタグを作らないので、ID は同点のときに外す順番にしか使わない
    ids_from_rows = bool(ids.duplicated().any())
    if ids_from_rows:
        ids = pd.Series(np.arange(len(ids)), index=ids.index)

    df = pd.DataFrame({
        "MasterAlignmentID": ids,
        "name": cols["name"],
        "mz": _to_float(cols["mz"], "Average Mz"),
        "rt": _to_float(cols["rt"], "Average Rt(min)"),
        "Ontology": cols["Ontology"],
        "Adduct_type": cols["Adduct_type"],
    }, index=table.frame.index)
    for col in ("name", "Ontology", "Adduct_type"):
        df[col] = df[col].where(df[col] != "null", "")
    df.attrs["ids_from_rows"] = ids_from_rows
    return df


def use_mask(frame: pd.DataFrame, column: str) -> pd.Series:
    """--use-column の列を真偽として読む。TRUE/1/yes が真、FALSE/0/no が偽。それ以外はエラー。"""
    pos = find_column(frame, column)
    name = frame.columns[pos]
    values = frame.iloc[:, pos].astype(str).str.strip().str.lower()
    bad = ~values.isin(_TRUE | _FALSE)
    if bad.any():
        raise _bad_values(name, frame.iloc[:, pos], bad, " TRUE/FALSE、1/0、yes/no 以外の値があります")
    return values.isin(_TRUE)


def apply_use_column(df: pd.DataFrame, use: pd.Series) -> pd.DataFrame:
    """use が偽の行を判定から外す(理由の最後に use_column を足す)。"""
    out = df.copy()
    off = ~use.to_numpy(dtype=bool)
    reasons = out["exclusion_reason"].astype(str)
    out.loc[off, "exclusion_reason"] = [f"{r};use_column" if r else "use_column" for r in reasons[off]]
    out.loc[off, "inclusion"] = False
    out.loc[off, "curation"] = False
    return out


@dataclass
class GenericColumns:
    rt: str = "RT"
    c: str = "C"
    db: str = "DB"
    ontology: str = "Ontology"


def _checked_numbers(values: pd.Series, checked: pd.Series, column: str, integer: bool) -> pd.Series:
    """checked の行だけ数値にする。有限の数値でない値(C・DB は整数値でない値も)はエラー。"""
    numbers = pd.to_numeric(values.where(checked), errors="coerce").astype(float)
    bad = checked & ~np.isfinite(numbers)
    if integer:
        bad |= checked & np.isfinite(numbers) & (numbers != numbers.round())
    if bad.any():
        kind = "整数" if integer else "数値"
        raise _bad_values(column, values, bad,
                          f"{kind}として読めない値があります(値が無い場合は空欄にしてください)")
    return numbers.where(checked)


def generic_frame(table: RawTable, cols: GenericColumns, use: pd.Series | None,
                  keep_ontology_null: bool) -> pd.DataFrame:
    """汎用 TSV から、pair_eval.run にそのまま渡せる表を作る。名前を使った処理は何もしない。"""
    frame = table.frame
    names = {}
    values = {}
    for key, name in (("rt", cols.rt), ("c", cols.c), ("db", cols.db), ("ontology", cols.ontology)):
        pos = find_column(frame, name)
        names[key] = frame.columns[pos]
        values[key] = frame.iloc[:, pos].astype(str).str.strip()

    ontology = values["ontology"]
    null_ontology = ontology.map(is_null_ontology).astype(bool)
    use_ok = pd.Series(True, index=frame.index) if use is None else pd.Series(
        use.to_numpy(dtype=bool), index=frame.index)
    missing = (values["rt"] == "") | (values["c"] == "") | (values["db"] == "")

    # 判定から外れる行(use が偽、Ontology なし)の値は確かめない。空のセルも確かめない
    candidate = use_ok & ~null_ontology
    rt = _checked_numbers(values["rt"], candidate & (values["rt"] != ""), names["rt"], integer=False)
    c = _checked_numbers(values["c"], candidate & (values["c"] != ""), names["c"], integer=True)
    db = _checked_numbers(values["db"], candidate & (values["db"] != ""), names["db"], integer=True)

    checks = [
        ("ontology_null", null_ontology & (not keep_ontology_null)),
        ("missing_value", missing),
        ("use_column", ~use_ok),
    ]
    reasons = [";".join(code for code, hit in checks if hit.iloc[i]) for i in range(len(frame))]

    df = pd.DataFrame({
        "MasterAlignmentID": np.arange(len(frame)),  # 同点のときは行番号の小さい方から外す
        "rt": rt,
        "c_s_total": c.round().astype("Int64"),
        "db_s_total": db.round().astype("Int64"),
        "Ontology": ontology,
        "base_ontology": ontology,
        "exclusion_reason": reasons,
    }, index=frame.index)
    df["inclusion"] = df["exclusion_reason"] == ""
    df["curation"] = df["inclusion"]
    return df


# 判定結果として右に足す列(rt_group は補正後の Ontology。入力にすでに Ontology 列があるので別の名前にする)
ADDED_COLUMNS = [
    "rt_group", "c_s_total", "db_s_total", "inclusion", "pair_evaluated", "curation", "exclusion_reason",
    "n_miss", "n_correct", "n_total", "miss%", "exclusion_round",
]


def check_added_columns(table: RawTable) -> None:
    """入力に、足す列と同じ名前の列があればエラー(上書きしない)。"""
    clash = [col for col in ADDED_COLUMNS if col in set(table.frame.columns)]
    if clash:
        raise TableError(f"入力に、判定結果の列と同じ名前の列があります: {', '.join(clash)}"
                         "(列名を変えてから実行してください)")


def write_table(table: RawTable, result: pd.DataFrame, path: Path, rt_order_only: bool = False) -> None:
    """入力の全列を元の文字列のまま左に置き、判定結果の列を右に足して書き出す。行は入力の順のまま。"""
    added = result.rename(columns={"Ontology": "rt_group"})[ADDED_COLUMNS]
    keep = (result["exclusion_reason"] == "rt_order") if rt_order_only else pd.Series(True, index=result.index)
    # 足す列の書式(NA は空欄など)は arf2 の TSV と同じ to_csv に任せる。入力の列は引用符を付けずにそのまま書く
    added_lines = added.loc[keep].to_csv(sep="\t", index=False, lineterminator="\n").split("\n")
    input_lines = ["\t".join(table.frame.columns)] + ["\t".join(row) for row in table.frame.loc[keep].itertuples(
        index=False, name=None)]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf_8_sig") as f:  # 改行は arf2 の TSV と同じく OS の既定
        f.write("".join(f"{a}\t{b}\n" for a, b in zip(input_lines, added_lines)))
