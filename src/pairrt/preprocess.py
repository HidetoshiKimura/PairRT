"""前処理: 除外フィルタ、炭素数・二重結合数、Ontology の補正。

4.4 の raw_data_for_order_2.py の移植。変更点は設計書 §5 を参照。
"""

import re
from dataclasses import dataclass

import pandas as pd

_C_DB = re.compile(r"(\d+):(\d+)")
_OXYGEN = re.compile(r";O\d*")  # MS-DIAL5 の表記 (;O, ;O2, ;O3 ...)
_ETHER_LINK = re.compile(r"(?<![A-Za-z])([OP])-\d")  # PE O-36:5 / PE P-36:4 の O- / P-


@dataclass
class FilterOptions:
    low_score: bool = True
    no_ms2: bool = True
    wo_ms2: bool = True
    riken: bool = True
    ontology_null: bool = True
    oxygen_tag: bool = True
    ether_split: bool = True


def calc_c_db(name: str) -> tuple[int | None, int | None]:
    """`|` より後ろ(なければ名前全体)の `数字:数字` を合計し、(炭素数, 二重結合数) を返す。"""
    name = str(name)
    if name in ("Others", "Unknown", ""):
        return None, None

    if "|" in name:
        name = name.split("|", 1)[1]

    matches = _C_DB.findall(name)
    if not matches:
        return None, None

    return sum(int(c) for c, _ in matches), sum(int(db) for _, db in matches)


def oxygen_tag(name: str) -> str:
    """種レベルの表記(`|` より前)から最初の `;O\\d*` を返す。無ければ空文字。"""
    m = _OXYGEN.search(str(name).split("|", 1)[0])
    return m.group() if m else ""


def ether_link(name: str) -> str:
    """種レベルの表記(`|` より前)の `O-` / `P-` から "O" か "P" を返す。無ければ空文字。"""
    m = _ETHER_LINK.search(str(name).split("|", 1)[0])
    return m.group(1) if m else ""


def is_null_ontology(ontology: str) -> bool:
    return str(ontology).strip().lower() in ("", "unknown", "null", "nan")


def _exclusion_reasons(name: str, ontology: str, opts: FilterOptions) -> str:
    checks = [
        ("ontology_null", opts.ontology_null and is_null_ontology(ontology)),
        ("low_score", opts.low_score and "low score" in name),
        ("no_MS2", opts.no_ms2 and "no MS2" in name),
        ("wo_MS2", opts.wo_ms2 and "w/o MS2" in name),
        ("RIKEN", opts.riken and "RIKEN" in name),
    ]
    return ";".join(code for code, hit in checks if hit)


def _corrected_ontology(name: str, ontology: str, opts: FilterOptions) -> str:
    if is_null_ontology(ontology):
        return ontology
    # P- は二重結合を 1 つ少なく書く(P-16:0 = O-16:1)ので、O- と同じグループで比べない
    if opts.ether_split and ontology.startswith("Ether"):
        link = ether_link(name)
        if link:
            ontology += f"_{link}"
    if opts.oxygen_tag:
        ontology += oxygen_tag(name)
    return ontology


def preprocess(spots: pd.DataFrame, opts: FilterOptions) -> pd.DataFrame:
    """load_spots の出力に、判定に使う列を足して返す。"""
    df = spots.copy()
    names = df["name"].fillna("").astype(str)
    base = df["Ontology"].fillna("").astype(str)

    df["base_ontology"] = base
    df["Ontology"] = [_corrected_ontology(n, o, opts) for n, o in zip(names, base)]

    c_db = [calc_c_db(n) for n in names]
    df["c_s_total"] = pd.array([c for c, _ in c_db], dtype="Int64")
    df["db_s_total"] = pd.array([db for _, db in c_db], dtype="Int64")

    df["exclusion_reason"] = [_exclusion_reasons(n, o, opts) for n, o in zip(names, base)]
    df["inclusion"] = df["exclusion_reason"] == ""
    df["curation"] = df["inclusion"]
    return df
