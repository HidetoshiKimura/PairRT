from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
SAMPLE_ARF2 = (
    ROOT.parent / "test" / "landscape_of_mice" / "raw" / "kidney" / "Neg"
    / "AlignmentResult_2026_10_20_11_07_17.arf2"
)


@pytest.fixture
def sample_arf2() -> Path:
    if not SAMPLE_ARF2.exists():
        pytest.skip(f"サンプル arf2 がありません: {SAMPLE_ARF2}")
    return SAMPLE_ARF2


def grp(rows) -> pd.DataFrame:
    """(MasterAlignmentID, c, db, rt) から、ペア評価にかける 1 グループ分の DataFrame を作る。"""
    return pd.DataFrame(
        [
            {"MasterAlignmentID": i, "Ontology": "PC", "base_ontology": "PC",
             "c_s_total": c, "db_s_total": db, "rt": rt,
             "inclusion": True, "curation": True, "exclusion_reason": ""}
            for i, c, db, rt in rows
        ]
    ).astype({"c_s_total": "Int64", "db_s_total": "Int64"})


def spots_df(rows) -> pd.DataFrame:
    """(MasterAlignmentID, name, Ontology) から load_spots と同じ形の DataFrame を作る。"""
    return pd.DataFrame(
        [
            {"MasterAlignmentID": i, "name": name, "mz": 500.0 + i, "rt": float(i),
             "Ontology": onto, "Adduct_type": "[M-H]-"}
            for i, name, onto in rows
        ]
    )
