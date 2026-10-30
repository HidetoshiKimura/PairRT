"""4.4 の auto_pair_eval_3 と、1 回分の数え方 (n_miss / n_correct) と反復除外の結果が一致することを確かめる。

反復除外の同点処理は設計書 §5-2 で意図的に変えたので、同点のラウンドがあるデータは比べない。
"""

import importlib
import sys

import numpy as np
import pandas as pd
import pytest

from conftest import ROOT, grp
from rtcurate.pair_eval import auto_curate_group, count_pairs

LEGACY_DIR = ROOT.parent / "RerativeRT_auto_checher4.4" / "script"

# 4.4 側の `drop(inplace=True)` が pandas 3 で出す非推奨警告(4.4 は変更しない)
pytestmark = pytest.mark.filterwarnings("ignore:The inplace keyword:DeprecationWarning",
                                        "ignore::pandas.errors.Pandas4Warning")


@pytest.fixture(scope="module")
def legacy():
    if not (LEGACY_DIR / "auto_pair_eval_3.py").exists():
        pytest.skip(f"4.4 のスクリプトがありません: {LEGACY_DIR}")
    sys.path.insert(0, str(LEGACY_DIR))
    try:
        yield importlib.import_module("auto_pair_eval_3")
    finally:
        sys.path.remove(str(LEGACY_DIR))


@pytest.mark.parametrize("seed", range(20))
def test_counts_match_legacy(legacy, seed):
    rng = np.random.default_rng(seed)
    c = rng.choice(np.arange(30, 45, 2), size=20)
    db = rng.integers(0, 7, size=20)
    rt = 0.3 * c - 0.8 * db + rng.normal(0, 1.5, size=20)
    g = grp([(i + 1, int(c[i]), int(db[i]), float(rt[i])) for i in range(20)])

    old_in = pd.DataFrame({
        "MasterAlignment_ID": g["MasterAlignmentID"], "c_s_total": c, "db_s_total": db,
        "rt": rt, "curation": True, "Ontology": "PC",
    })
    pairs, df_true, _ = legacy.pair_evaluation(old_in)
    old = legacy.miss_count(pairs, df_true).set_index("MasterAlignment_ID")

    new = count_pairs(g).set_axis(g["MasterAlignmentID"])
    for col in ["n_miss", "n_correct"]:
        assert new[col].to_dict() == old[col].astype(int).to_dict(), col


def _has_tie(g: pd.DataFrame, curated: pd.DataFrame, n_miss_threshold: int, pct: float) -> bool:
    """新実装の除外順をなぞり、どこかのラウンドで n_miss 最大の候補が同点なら True。"""
    rounds = curated["exclusion_round"]
    for r in range(1, int(rounds.max()) + 1 if rounds.notna().any() else 1):
        state = g.assign(curation=~(rounds < r).fillna(False).to_numpy())
        counts = count_pairs(state)
        cand = counts.loc[(counts["n_miss"] >= n_miss_threshold) & (counts["miss%"] >= pct)]
        if (cand["n_miss"] == cand["n_miss"].max()).sum() > 1:
            return True
    return False


@pytest.mark.parametrize("seed", range(30))
def test_auto_curation_matches_legacy(legacy, seed):
    """同点が起きないデータで、反復除外の結果(外す ID と残った行の数え方)が 4.4 と一致する。"""
    rng = np.random.default_rng(seed)
    n = 25
    c = rng.choice(np.arange(30, 45, 2), size=n)
    db = rng.integers(0, 7, size=n)
    rt = 0.3 * c - 0.8 * db + rng.normal(0, 0.3, size=n)
    outliers = rng.choice(n, size=3, replace=False)
    rt[outliers] += rng.choice([-1, 1], size=3) * rng.uniform(3, 6, size=3)
    g = grp([(i + 1, int(c[i]), int(db[i]), float(rt[i])) for i in range(n)])

    new = auto_curate_group(g, 4, 10.0)
    if _has_tie(g, new, 4, 10.0):
        pytest.skip("同点のラウンドがある(§5-2 で意図的に 4.4 と変えた部分)")

    old_in = pd.DataFrame({
        "MasterAlignment_ID": g["MasterAlignmentID"], "c_s_total": c, "db_s_total": db,
        "rt": rt, "curation": True, "Ontology": "PC",
    })
    old_true, _, old_excluded = legacy.auto_curation(old_in, 4, 10, [])
    old_ex_ids = [int(x["MasterAlignment_ID"].iloc[0]) for x in old_excluded]

    new_ex = new.loc[new["exclusion_reason"] == "rt_order"].sort_values("exclusion_round")
    assert new_ex["MasterAlignmentID"].tolist() == old_ex_ids  # 外す順番まで一致

    kept = new.loc[new["curation"].astype(bool)].set_index("MasterAlignmentID")
    old_kept = old_true.set_index("MasterAlignment_ID")
    for col in ["n_miss", "n_correct"]:
        assert kept[col].astype(int).to_dict() == old_kept[col].astype(int).to_dict(), col
