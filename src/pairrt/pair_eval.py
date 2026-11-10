"""ペア評価と反復除外。

4.4 の auto_pair_eval_3.py の移植。同じ Ontology の中で、炭素数・二重結合数から
予想される溶出順と実際の RT が矛盾するペアを数え、矛盾の多いものから 1 件ずつ外す。
"""

from collections.abc import Iterable

import numpy as np
import pandas as pd

from pairrt.preprocess import is_null_ontology

COUNT_COLUMNS = ["n_miss", "n_correct", "n_total", "miss%"]
_RESULT_DTYPES = {"n_miss": "Int64", "n_correct": "Int64", "n_total": "Int64",
                  "miss%": "Float64", "exclusion_round": "Int64"}


def count_pairs(group: pd.DataFrame) -> pd.DataFrame:
    """`curation=True` の行どうしの全ペアを判定し、行ごとの n_miss などを返す(index は元のまま)。"""
    active = group.loc[group["curation"].astype(bool)]
    c = active["c_s_total"].to_numpy(dtype=float)
    db = active["db_s_total"].to_numpy(dtype=float)
    rt = active["rt"].to_numpy(dtype=float)

    ci, cj = c[:, None], c[None, :]
    dbi, dbj = db[:, None], db[None, :]
    rti, rtj = rt[:, None], rt[None, :]

    fast = (ci <= cj) & (dbi >= dbj) & ((ci != cj) | (dbi != dbj))  # i が先に溶出するはず
    slow = fast.T
    miss = ((rti >= rtj) & fast) | ((rti <= rtj) & slow)
    correct = ((rti < rtj) & fast) | ((rti > rtj) & slow)

    n_miss = miss.sum(axis=1)
    n_correct = correct.sum(axis=1)
    n_total = n_miss + n_correct
    miss_pct = np.where(n_total > 0, 100.0 * n_miss / np.maximum(n_total, 1), 0.0)

    return pd.DataFrame(
        {"n_miss": n_miss, "n_correct": n_correct, "n_total": n_total, "miss%": miss_pct},
        index=active.index,
    )


def auto_curate_group(group: pd.DataFrame, n_miss_threshold: int,
                      miss_pct_threshold: float) -> pd.DataFrame:
    """1 つの Ontology グループで、候補が無くなるまで最悪の 1 行を外し続ける。"""
    g = group.copy()
    for col, dtype in _RESULT_DTYPES.items():
        g[col] = pd.Series(pd.NA, index=g.index, dtype=dtype)

    round_no = 0
    while True:
        counts = count_pairs(g)
        candidates = counts.loc[(counts["n_miss"] >= n_miss_threshold)
                                & (counts["miss%"] >= miss_pct_threshold)]
        if candidates.empty:
            g.loc[counts.index, COUNT_COLUMNS] = counts[COUNT_COLUMNS]
            return g

        round_no += 1
        ranked = candidates.assign(_id=g.loc[candidates.index, "MasterAlignmentID"]).sort_values(
            ["n_miss", "miss%", "_id"], ascending=[False, False, True])
        worst = ranked.index[0]
        g.loc[worst, COUNT_COLUMNS] = counts.loc[worst, COUNT_COLUMNS].to_numpy()
        g.loc[worst, "exclusion_round"] = round_no
        g.loc[worst, "curation"] = False
        g.loc[worst, "exclusion_reason"] = "rt_order"


def run(df: pd.DataFrame, n_miss_threshold: int = 4, miss_pct_threshold: float = 10.0,
        exclude_ontologies: Iterable[str] = ()) -> pd.DataFrame:
    """前処理済みの全スポットに対し、Ontology ごとにペア評価と反復除外を行う。"""
    out = df.copy()
    out["pair_evaluated"] = False
    for col, dtype in _RESULT_DTYPES.items():
        out[col] = pd.Series(pd.NA, index=out.index, dtype=dtype)

    excluded = set(exclude_ontologies)
    target = (
        out["inclusion"].astype(bool)
        & ~out["base_ontology"].map(is_null_ontology).astype(bool)
        & ~out["base_ontology"].isin(excluded)
        & ~out["Ontology"].isin(excluded)
        & out["c_s_total"].notna()
        & out["db_s_total"].notna()
        & out["rt"].notna()
    )
    out.loc[target, "pair_evaluated"] = True

    result_cols = [*COUNT_COLUMNS, "exclusion_round", "curation", "exclusion_reason"]
    for _, group in out.loc[target].groupby("Ontology", sort=True):
        curated = auto_curate_group(group, n_miss_threshold, miss_pct_threshold)
        # 複数列を一度に .loc 代入すると、pandas 3.1 では NA を含む拡張型の列で失敗するので 1 列ずつ書く
        for col in result_cols:
            out.loc[curated.index, col] = curated[col]
    return out
