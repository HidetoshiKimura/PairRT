"""判定結果の TSV 書き出しと、件数の内訳表示。"""

from pathlib import Path

import pandas as pd

MAIN_COLUMNS = [
    "MasterAlignmentID", "name", "mz", "rt", "Ontology", "base_ontology", "Adduct_type",
    "c_s_total", "db_s_total", "inclusion", "pair_evaluated", "curation", "exclusion_reason",
    "n_miss", "n_correct", "n_total", "miss%", "exclusion_round",
]


def _sorted(df: pd.DataFrame) -> pd.DataFrame:
    """n_miss 降順(NA は最後)、同点は MasterAlignmentID 昇順。"""
    return df.sort_values(["n_miss", "MasterAlignmentID"], ascending=[False, True],
                          na_position="last", kind="stable")[MAIN_COLUMNS]


def _write(df: pd.DataFrame, path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    _sorted(df).to_csv(path, sep="\t", index=False, encoding="utf_8_sig")


def write_main(df: pd.DataFrame, path: Path) -> None:
    _write(df, path)


def write_excluded(df: pd.DataFrame, path: Path) -> None:
    """RT 順序で自動除外した行だけを書き出す。"""
    _write(df.loc[df["exclusion_reason"] == "rt_order"], path)


def summary(df: pd.DataFrame) -> str:
    pair_evaluated = df["pair_evaluated"].astype(bool)
    rows = [
        ("検出スポット数", len(df)),
        ("前処理フィルタ後", int(df["inclusion"].astype(bool).sum())),
        ("ペア評価対象", int(pair_evaluated.sum())),
        ("RT順序で自動除外", int((df["exclusion_reason"] == "rt_order").sum())),
        ("ペア評価で残ったもの", int((pair_evaluated & df["curation"].astype(bool)).sum())),
    ]
    return "\n".join(f"{label}: {value}" for label, value in rows)
