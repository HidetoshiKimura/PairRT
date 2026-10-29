import pandas as pd
import pytest

from conftest import spots_df
from rtcurate.output import MAIN_COLUMNS, summary, write_excluded, write_main
from rtcurate.pair_eval import run
from rtcurate.preprocess import FilterOptions, preprocess


@pytest.fixture
def evaluated_df() -> pd.DataFrame:
    # PC 系列 8 件 + 外れ値 1 件(rt_order)、前処理除外 1 件、c/db 無しの未評価 1 件
    rows = [(i, f"PC {30 + 2 * i}:0", "PC") for i in range(1, 9)]
    rows += [(9, "PC 48:0", "PC"), (10, "low score: PE 34:1", "PE"), (11, "TCA", "BileAcid")]
    df = preprocess(spots_df(rows), FilterOptions())
    df["rt"] = [3.0 + i for i in range(1, 9)] + [1.0, 5.0, 6.0]
    return run(df)


def test_write_main_order_and_bom(tmp_path, evaluated_df):
    p = tmp_path / "new_dir" / "out.tsv"
    write_main(evaluated_df, p)
    assert p.read_bytes().startswith(b"\xef\xbb\xbf")
    back = pd.read_csv(p, sep="\t", encoding="utf_8_sig")
    assert list(back.columns) == MAIN_COLUMNS and len(back) == len(evaluated_df)
    n = back["n_miss"]
    assert n.dropna().is_monotonic_decreasing
    assert n.isna().sum() == 0 or n.isna().idxmax() > n.last_valid_index()  # NA は最後にまとまる
    assert back.loc[0, "MasterAlignmentID"] == 9  # n_miss 最大の外れ値が先頭


def test_write_main_ties_sorted_by_id(tmp_path, evaluated_df):
    p = tmp_path / "out.tsv"
    write_main(evaluated_df, p)
    back = pd.read_csv(p, sep="\t", encoding="utf_8_sig")
    zero = back.loc[back["n_miss"] == 0, "MasterAlignmentID"]
    assert zero.is_monotonic_increasing


def test_write_excluded_only_rt_order(tmp_path, evaluated_df):
    p = tmp_path / "ex.tsv"
    write_excluded(evaluated_df, p)
    back = pd.read_csv(p, sep="\t", encoding="utf_8_sig")
    assert list(back.columns) == MAIN_COLUMNS
    assert back["MasterAlignmentID"].tolist() == [9]
    assert (back["exclusion_reason"] == "rt_order").all()


def test_summary_labels_and_values(evaluated_df):
    s = summary(evaluated_df)
    for label, value in [("検出スポット数", 11), ("前処理フィルタ後", 10), ("ペア評価対象", 9),
                         ("RT順序で自動除外", 1), ("ペア評価で残ったもの", 8)]:
        assert f"{label}: {value}" in s
