from pathlib import Path

import pandas as pd
import pytest

from rtcurate.table import (
    TableError, apply_use_column, find_column, msdial_spots, read_table, use_mask,
)

# MS-DIAL5 の Height_1_….txt と同じ並び: メタ情報の行(1 列目が空)のあとに列名の行
HEADER = ["Alignment ID", "Average Rt(min)", "Average Mz", "Metabolite name", "Adduct type",
          "Ontology", "1", "1"]  # 末尾の "1" はクラス平均と標準偏差(同じ名前が 2 つ)
META4 = [["", "", "", "", "", "Class", "1", "1"], ["", "", "", "", "", "File type", "Sample", "Sample"],
         ["", "", "", "", "", "Injection order", "1", "2"], ["", "", "", "", "", "Batch ID", "1", "1"]]
META9 = META4 + [["", "", "", "", "", f"meta{i}", "", ""] for i in range(5)]
ROW_A = ["0", "13.793", "760.5851", "PC 34:1", "[M+H]+", "PC", "100", "5"]
ROW_B = ["1", "14.100", "788.6164", "PC 36:1", "[M+H]+", "PC", "200", "6"]


def write_tsv(path, lines: list[list[str]], encoding="utf-8") -> Path:
    path.write_text("\n".join("\t".join(r) for r in lines) + "\n", encoding=encoding)
    return path


def test_msdial_txt_detected_with_4_meta_rows(tmp_path):
    t = read_table(write_tsv(tmp_path / "Height_1.txt", META4 + [HEADER, ROW_A, ROW_B]))
    assert t.kind == "msdial_txt"
    assert list(t.frame.columns) == HEADER  # "1" が 2 つのまま
    assert len(t.frame) == 2 and t.frame.iloc[0, 1] == "13.793"
    assert list(t.frame.index) == [0, 1]


def test_msdial_txt_detected_with_9_meta_rows(tmp_path):
    t = read_table(write_tsv(tmp_path / "Height_1.txt", META9 + [HEADER, ROW_A, ROW_B]))
    assert t.kind == "msdial_txt" and list(t.frame.columns) == HEADER and len(t.frame) == 2


def test_header_after_row_20_is_generic(tmp_path):
    meta = [["", "", "", "", "", f"m{i}", "", ""] for i in range(20)]
    t = read_table(write_tsv(tmp_path / "late.txt", meta + [HEADER, ROW_A]))
    assert t.kind == "generic"


def test_alignment_id_row_without_required_columns_is_generic(tmp_path):
    header = [h if h != "Ontology" else "Class" for h in HEADER]
    t = read_table(write_tsv(tmp_path / "x.txt", META4 + [header, ROW_A]))
    assert t.kind == "generic"


def test_generic_uses_first_row_as_header(tmp_path):
    t = read_table(write_tsv(tmp_path / "x.tsv", [["RT", "C", "DB", "Ontology"], ["1.50", "34", "1", "PC"]]))
    assert t.kind == "generic"
    assert list(t.frame.columns) == ["RT", "C", "DB", "Ontology"]
    assert t.frame.iloc[0, 0] == "1.50"  # 文字列のまま


def test_cp932_fallback(tmp_path):
    t = read_table(write_tsv(tmp_path / "x.tsv", [["RT", "名前"], ["1.0", "ホスファチジルコリン"]],
                             encoding="cp932"))
    assert list(t.frame.columns) == ["RT", "名前"] and t.frame.iloc[0, 1] == "ホスファチジルコリン"


def test_utf8_bom(tmp_path):
    t = read_table(write_tsv(tmp_path / "x.tsv", [["RT", "C"], ["1", "2"]], encoding="utf-8-sig"))
    assert list(t.frame.columns) == ["RT", "C"]


def test_undecodable_raises(tmp_path):
    p = tmp_path / "x.tsv"
    p.write_bytes(b"RT\tC\n\x81\x7f\xff\xfe\n")
    with pytest.raises(TableError):
        read_table(p)


def test_row_with_wrong_column_count_raises(tmp_path):
    p = write_tsv(tmp_path / "x.tsv", [["RT", "C"], ["1", "2"], ["1"]])
    with pytest.raises(TableError, match="3"):  # ファイルの 3 行目
        read_table(p)


def test_row_with_too_many_columns_raises(tmp_path):
    p = write_tsv(tmp_path / "x.tsv", [["RT", "C"], ["1", "2", "3"]])
    with pytest.raises(TableError, match="2"):
        read_table(p)


def test_trailing_blank_lines_are_ignored(tmp_path):
    p = tmp_path / "x.tsv"
    p.write_text("RT\tC\n1\t2\n\n\n", encoding="utf-8")
    assert len(read_table(p).frame) == 1


def test_empty_file_raises(tmp_path):
    p = tmp_path / "x.tsv"
    p.write_text("", encoding="utf-8")
    with pytest.raises(TableError):
        read_table(p)


def test_find_column_case_insensitive():
    f = pd.DataFrame([["1", "2"]], columns=["rt", "Ontology"])
    assert find_column(f, "RT") == 0
    assert find_column(f, "ontology") == 1


def test_find_column_missing_or_ambiguous():
    with pytest.raises(TableError, match="RT"):
        find_column(pd.DataFrame(columns=["C"]), "RT")
    with pytest.raises(TableError, match="RT"):
        find_column(pd.DataFrame(columns=["RT", "rt"]), "RT")
    with pytest.raises(TableError, match="RT"):
        find_column(pd.DataFrame(columns=["RT", "RT"]), "RT")


def msdial(tmp_path, rows, header=HEADER) -> "RawTable":
    return read_table(write_tsv(tmp_path / "h.txt", META4 + [header] + rows))


def test_msdial_spots_matches_load_spots_shape(tmp_path):
    t = msdial(tmp_path, [["5", "13.793", "760.5851", "PC 34:1", "[M+H]+", "PC", "1", "2"],
                          ["7", "null", "null", "null", "null", "null", "1", "2"]])
    df = msdial_spots(t)
    assert list(df.columns) == ["MasterAlignmentID", "name", "mz", "rt", "Ontology", "Adduct_type"]
    assert df.loc[0, "MasterAlignmentID"] == 5 and df.loc[0, "rt"] == 13.793
    assert df.loc[0, "mz"] == 760.5851 and df.loc[0, "name"] == "PC 34:1"
    assert df.loc[0, "Adduct_type"] == "[M+H]+" and df.loc[0, "Ontology"] == "PC"
    assert pd.isna(df.loc[1, "rt"]) and pd.isna(df.loc[1, "mz"])
    assert df.loc[1, "name"] == "" and df.loc[1, "Ontology"] == "" and df.loc[1, "Adduct_type"] == ""
    assert df["MasterAlignmentID"].dtype.kind == "i" and df["rt"].dtype.kind == "f"


@pytest.mark.parametrize("rows,header,match", [
    ([ROW_A, ROW_A], HEADER, "重複"),
    ([["x"] + ROW_A[1:]], HEADER, "Alignment ID"),
    ([["0", "abc"] + ROW_A[2:]], HEADER, "Average Rt"),
    ([["0", "1.0", "abc"] + ROW_A[3:]], HEADER, "Average Mz"),
    ([ROW_A[:4] + ROW_A[5:]], [h for h in HEADER if h != "Adduct type"], "Adduct type"),
    ([ROW_A + ["PC"]], HEADER + ["Ontology"], "Ontology"),
])
def test_msdial_spots_errors(tmp_path, rows, header, match):
    with pytest.raises(TableError, match=match):
        msdial_spots(msdial(tmp_path, rows, header))


@pytest.mark.parametrize("value,expected", [("TRUE", True), ("true", True), (" 1 ", True), ("Yes", True),
                                            ("FALSE", False), ("0", False), ("no", False)])
def test_use_mask_values(value, expected):
    f = pd.DataFrame({"Keep": [value]})
    assert use_mask(f, "keep").tolist() == [expected]


def test_use_mask_rejects_blank_and_other():
    f = pd.DataFrame({"keep": ["TRUE", "", "maybe"]})
    with pytest.raises(TableError, match="maybe") as e:
        use_mask(f, "keep")
    assert "keep" in str(e.value)


def test_use_mask_missing_column():
    with pytest.raises(TableError, match="keep"):
        use_mask(pd.DataFrame({"x": ["1"]}), "keep")


def test_apply_use_column_appends_reason():
    df = pd.DataFrame({"exclusion_reason": ["", "low_score", ""], "inclusion": [True, False, True],
                       "curation": [True, False, True]})
    out = apply_use_column(df, pd.Series([True, False, False]))
    assert out["exclusion_reason"].tolist() == ["", "low_score;use_column", "use_column"]
    assert out["inclusion"].tolist() == [True, False, False]
    assert out["curation"].tolist() == [True, False, False]
    assert df["exclusion_reason"].tolist() == ["", "low_score", ""]  # 元の表は変えない
