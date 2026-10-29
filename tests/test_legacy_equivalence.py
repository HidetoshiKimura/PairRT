"""4.4 の auto_pair_eval_3 と、1 回分の数え方 (n_miss / n_correct) が一致することを確かめる。

反復除外の同点処理は設計書 §5-2 で意図的に変えたので、ここでは比べない。
"""

import importlib
import sys

import numpy as np
import pandas as pd
import pytest

from conftest import ROOT, grp
from rtcurate.pair_eval import count_pairs

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
