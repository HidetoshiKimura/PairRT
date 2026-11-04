import pytest

from conftest import spots_df
from rtcurate.preprocess import FilterOptions, calc_c_db, is_null_ontology, oxygen_tag, preprocess


@pytest.mark.parametrize("name,expected", [
    ("PC 34:1|PC 16:0_18:1", (34, 1)),
    ("Cer 33:0;O3|Cer 19:0;O3/14:0", (33, 0)),
    ("TG 18:1_18:1_18:1", (54, 3)),
    ("TCA", (None, None)),
    ("Unknown", (None, None)),
    ("", (None, None)),
])
def test_calc_c_db(name, expected):
    assert calc_c_db(name) == expected


@pytest.mark.parametrize("name,tag", [
    ("OxPC 34:2;O2|OxPC 16:0_18:2;O2", ";O2"),
    ("SM 34:1;O2(d9)", ";O2"),
    ("FA 18:2;O", ";O"),
    ("Cer 33:0;O3|Cer 19:0;O2/14:0;O", ";O3"),
    ("PC 34:1|PC 16:0_18:1", ""),
])
def test_oxygen_tag(name, tag):
    assert oxygen_tag(name) == tag


@pytest.mark.parametrize("onto,expected", [
    ("", True), ("unknown", True), ("Unknown", True), (" null ", True), ("PC", False),
])
def test_is_null_ontology(onto, expected):
    assert is_null_ontology(onto) is expected


def test_filters_and_reasons():
    df = spots_df([
        (1, "PC 34:1|PC 16:0_18:1", "PC"),
        (2, "low score: no MS2: PC 30:0", "PC"),
        (3, "RIKEN N-VS1 ID-1", "Unknown"),
        (4, "OxPC 34:2;O2", "OxPC"),
        (5, "w/o MS2: PE 36:2", "PE"),
    ])
    out = preprocess(df, FilterOptions()).set_index("MasterAlignmentID")
    assert out.loc[1, "inclusion"] and out.loc[1, "exclusion_reason"] == ""
    assert out.loc[1, "c_s_total"] == 34 and out.loc[1, "db_s_total"] == 1
    assert out.loc[2, "exclusion_reason"] == "low_score;no_MS2"
    assert out.loc[3, "exclusion_reason"] == "ontology_null;RIKEN"
    assert out.loc[4, "Ontology"] == "OxPC;O2" and out.loc[4, "base_ontology"] == "OxPC"
    assert out.loc[5, "exclusion_reason"] == "wo_MS2"
    assert (out["curation"] == out["inclusion"]).all()


def test_no_c_db_is_na():
    out = preprocess(spots_df([(1, "TCA", "BileAcid")]), FilterOptions())
    assert out["c_s_total"].isna().iloc[0] and out["db_s_total"].isna().iloc[0]
    assert out["inclusion"].iloc[0]


def test_keep_flags_disable_filters():
    df = spots_df([(2, "low score: PC 30:0", "PC")])
    out = preprocess(df, FilterOptions(low_score=False))
    assert out["inclusion"].iloc[0] and out["exclusion_reason"].iloc[0] == ""


def test_no_oxygen_tag_option():
    out = preprocess(spots_df([(4, "OxPC 34:2;O2", "OxPC")]), FilterOptions(oxygen_tag=False))
    assert out["Ontology"].iloc[0] == "OxPC"


def test_p_and_o_tags():
    df = spots_df([(1, "PE P-36:2|PE P-18:0_18:2", "EtherPE"), (2, "PC O-34:1", "EtherPC")])
    out = preprocess(df, FilterOptions(p_tag=True, o_tag=True)).set_index("MasterAlignmentID")
    assert out.loc[1, "Ontology"] == "EtherPE"  # "p-" は小文字のみ一致(4.4 と同じ)
    assert out.loc[2, "Ontology"] == "EtherPCO-"


def test_null_ontology_not_tagged():
    out = preprocess(spots_df([(1, "FA 18:2;O", "")]), FilterOptions(ontology_null=False))
    assert out["Ontology"].iloc[0] == "" and out["inclusion"].iloc[0]


# 各フィルタにだけ当たる行(名前, Ontology)と、そのフィルタの FilterOptions の項目・理由コード
FILTER_CASES = [
    ("low_score", "low score: PC 30:0", "PC", "low_score"),
    ("no_ms2", "no MS2: PC 30:0", "PC", "no_MS2"),
    ("wo_ms2", "w/o MS2: PC 30:0", "PC", "wo_MS2"),
    ("riken", "RIKEN N-VS1 ID-1", "PC", "RIKEN"),
    ("ontology_null", "PC 30:0", "", "ontology_null"),
]


@pytest.mark.parametrize("field,name,onto,code", FILTER_CASES)
def test_each_filter_excludes_by_default(field, name, onto, code):
    out = preprocess(spots_df([(1, name, onto)]), FilterOptions())
    assert not out["inclusion"].iloc[0] and out["exclusion_reason"].iloc[0] == code


@pytest.mark.parametrize("field,name,onto,code", FILTER_CASES)
def test_each_filter_can_be_disabled(field, name, onto, code):
    out = preprocess(spots_df([(1, name, onto)]), FilterOptions(**{field: False}))
    assert out["inclusion"].iloc[0] and out["curation"].iloc[0]
    assert out["exclusion_reason"].iloc[0] == ""


@pytest.mark.parametrize("field,name,onto,code", FILTER_CASES)
def test_disabling_one_filter_keeps_the_others(field, name, onto, code):
    # 5 つのフィルタすべてに当たる行で、1 つだけ外すと、残りの 4 つの理由は残る
    out = preprocess(spots_df([(1, "low score: no MS2: w/o MS2: RIKEN PC 30:0", "")]),
                     FilterOptions(**{field: False}))
    reasons = out["exclusion_reason"].iloc[0].split(";")
    assert code not in reasons
    assert sorted(reasons) == sorted(c for *_, c in FILTER_CASES if c != code)
    assert not out["inclusion"].iloc[0]


def test_disable_all_filters():
    off = {f: False for f, *_ in FILTER_CASES}
    df = spots_df([(i, n, o) for i, (_, n, o, _) in enumerate(FILTER_CASES, start=1)])
    out = preprocess(df, FilterOptions(**off))
    assert out["inclusion"].all() and (out["exclusion_reason"] == "").all()


def test_kept_low_score_row_still_gets_c_db():
    # フィルタを外して残した行も、名前の接頭辞に関係なく c/db が取れてペア評価に回る
    out = preprocess(spots_df([(1, "low score: no MS2: PC 34:1|PC 16:0_18:1", "PC")]),
                     FilterOptions(low_score=False, no_ms2=False))
    assert out["c_s_total"].iloc[0] == 34 and out["db_s_total"].iloc[0] == 1
