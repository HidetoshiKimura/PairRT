import pandas as pd

from conftest import grp, spots_df
from pairrt.pair_eval import auto_curate_group, count_pairs, run
from pairrt.preprocess import FilterOptions, preprocess


def by_id(g: pd.DataFrame, result: pd.DataFrame) -> pd.DataFrame:
    return result.set_axis(g.loc[result.index, "MasterAlignmentID"])


def test_counts_simple_series():
    # PC 32:0, 34:0, 36:0 は rt 昇順が正しい。36:0 だけ早すぎる。
    g = grp([(1, 32, 0, 5.0), (2, 34, 0, 6.0), (3, 36, 0, 4.0)])
    r = by_id(g, count_pairs(g))
    assert r.loc[3, "n_miss"] == 2 and r.loc[3, "n_correct"] == 0
    assert r.loc[1, "n_miss"] == 1 and r.loc[1, "n_correct"] == 1
    assert r.loc[1, "miss%"] == 50.0


def test_more_double_bonds_elute_earlier():
    # 同じ炭素数で DB が多い方が先に出るのが正しい
    r = count_pairs(grp([(1, 34, 2, 4.0), (2, 34, 1, 5.0)]))
    assert r["n_correct"].tolist() == [1, 1] and r["n_miss"].tolist() == [0, 0]


def test_equal_rt_counts_as_miss():
    r = count_pairs(grp([(1, 32, 0, 5.0), (2, 34, 0, 5.0)]))
    assert r["n_miss"].tolist() == [1, 1]


def test_unordered_pair_not_counted():
    # c 多い & db 多い → 順序不定
    r = count_pairs(grp([(1, 34, 1, 5.0), (2, 36, 2, 4.0)]))
    assert r["n_total"].tolist() == [0, 0] and r["miss%"].tolist() == [0.0, 0.0]


def test_identical_composition_not_counted():
    r = count_pairs(grp([(1, 34, 1, 5.0), (2, 34, 1, 9.0)]))
    assert r["n_total"].tolist() == [0, 0]


def test_count_pairs_ignores_curation_false_rows():
    g = grp([(1, 32, 0, 5.0), (2, 34, 0, 4.0), (3, 36, 0, 6.0)])
    g.loc[g["MasterAlignmentID"] == 2, "curation"] = False
    r = by_id(g, count_pairs(g))
    assert list(r.index) == [1, 3] and r["n_miss"].tolist() == [0, 0]


def test_auto_curation_removes_outlier():
    rows = [(i, 30 + 2 * i, 0, 3.0 + i) for i in range(1, 9)] + [(99, 48, 0, 1.0)]
    out = auto_curate_group(grp(rows), 4, 10.0).set_index("MasterAlignmentID")
    assert out.loc[99, "curation"] == False and out.loc[99, "exclusion_reason"] == "rt_order"
    assert out.loc[99, "exclusion_round"] == 1 and out.loc[99, "n_miss"] == 8
    rest = out.drop(index=99)
    assert rest["curation"].all() and (rest["n_miss"] == 0).all() and rest["exclusion_round"].isna().all()


def test_tie_breaks_by_miss_pct():
    # ID10 と ID5 はどちらも n_miss=1。ID5 は ID7 と正しい順序なので miss% が低い → ID10 が外れる
    g = grp([(5, 34, 0, 4.0), (7, 34, 1, 3.0), (10, 32, 0, 5.0)])
    out = auto_curate_group(g, 1, 0.0).set_index("MasterAlignmentID")
    assert out["curation"].to_dict() == {5: True, 7: True, 10: False}
    assert out.loc[10, "miss%"] == 100.0 and out.loc[10, "exclusion_round"] == 1


def test_tie_breaks_by_id_when_pct_equal():
    g = grp([(2, 34, 0, 4.0), (1, 32, 0, 5.0)])
    out = auto_curate_group(g, 1, 0.0).set_index("MasterAlignmentID")
    assert out["curation"].to_dict() == {1: False, 2: True}


def test_single_spot_group_and_empty_run():
    df = preprocess(spots_df([(1, "PC 34:1", "PC"), (2, "low score: PE 1:0", "PE")]), FilterOptions())
    out = run(df).set_index("MasterAlignmentID")
    assert out.loc[1, "pair_evaluated"] and out.loc[1, "n_total"] == 0 and out.loc[1, "curation"]
    assert not out.loc[2, "pair_evaluated"] and pd.isna(out.loc[2, "n_miss"])
    assert len(run(df.assign(inclusion=False, curation=False))) == 2


def test_no_c_db_untouched():
    df = preprocess(spots_df([(1, "TCA", "BileAcid")]), FilterOptions())
    out = run(df)
    assert not out["pair_evaluated"].iloc[0] and out["curation"].iloc[0]
    assert out["exclusion_reason"].iloc[0] == ""


def test_null_ontology_never_pair_evaluated():
    df = preprocess(spots_df([(1, "PC 34:1", "Unknown"), (2, "PC 36:1", "Unknown")]),
                    FilterOptions(ontology_null=False))
    assert not run(df)["pair_evaluated"].any()


def test_exclude_ontology_matches_base_and_corrected():
    df = preprocess(spots_df([(1, "OxPC 34:2;O2", "OxPC"), (2, "OxPC 36:2;O2", "OxPC")]), FilterOptions())
    assert not run(df, exclude_ontologies=["OxPC"])["pair_evaluated"].any()
    assert not run(df, exclude_ontologies=["OxPC;O2"])["pair_evaluated"].any()
    assert run(df)["pair_evaluated"].all()


def test_run_with_exclusion_in_non_first_group():
    # 実データでは Ontology グループの行が飛び飛びの index を持つ。そこで除外が起きても落ちないこと
    names = [("PE 34:0", "PE")] + [(f"PC {30 + 2 * i}:0", "PC") for i in range(1, 9)] + [("PC 48:0", "PC")]
    df = preprocess(spots_df([(i, n, o) for i, (n, o) in enumerate(names, start=1)]), FilterOptions())
    df["rt"] = [2.0] + [3.0 + i for i in range(1, 9)] + [1.0]
    out = run(df).set_index("MasterAlignmentID")
    assert out.loc[10, "exclusion_reason"] == "rt_order" and out.loc[10, "n_miss"] == 8
    assert out.loc[10, "exclusion_round"] == 1 and not out.loc[10, "curation"]
    assert out.loc[2, "curation"] and out.loc[2, "n_miss"] == 0
    assert str(out["n_miss"].dtype) == "Int64" and str(out["miss%"].dtype) == "Float64"


def test_groups_are_independent():
    # 同じ組成でも Ontology が違えば比較しない
    df = preprocess(spots_df([(1, "PC 32:0", "PC"), (2, "PE 34:0", "PE")]), FilterOptions())
    df["rt"] = [5.0, 1.0]
    out = run(df)
    assert (out["n_total"] == 0).all()
