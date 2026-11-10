from pathlib import Path

import pandas as pd
import pytest

from pairrt.pair_eval import run
from pairrt.preprocess import FilterOptions, preprocess
from pairrt.table import (
    ADDED_COLUMNS, GenericColumns, TableError, apply_use_column, find_column, generic_frame, msdial_spots, read_table,
    check_added_columns, use_mask, write_table,
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


def test_msdial_spots_unique_ids_are_used(tmp_path):
    df = msdial_spots(msdial(tmp_path, [["87"] + ROW_A[1:], ["38"] + ROW_B[1:]]))
    assert df["MasterAlignmentID"].tolist() == [87, 38] and not df.attrs.get("ids_from_rows")


def test_msdial_spots_duplicate_ids_use_row_number(tmp_path):
    # Pos と Neg をまとめた表などで Alignment ID が重複していたら、行番号(0 始まり)を ID の代わりに使う
    df = msdial_spots(msdial(tmp_path, [["87"] + ROW_A[1:], ["38"] + ROW_B[1:], ["87"] + ROW_B[1:]]))
    assert df["MasterAlignmentID"].tolist() == [0, 1, 2] and df.attrs["ids_from_rows"]


@pytest.mark.parametrize("rows,header,match", [
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


GEN_HEADER = ["ID", "RT", "C", "DB", "Ontology", "Name"]


def gen(tmp_path, rows, header=GEN_HEADER):
    return read_table(write_tsv(tmp_path / "g.tsv", [header] + rows))


def test_generic_basic(tmp_path):
    t = gen(tmp_path, [["a", "1.5", "34", "1", "PC", "low score: PC 34:1"],
                       ["b", "2.5", "36", "2", "PC", "PC 36:2"]])
    df = generic_frame(t, GenericColumns(), None, keep_ontology_null=False)
    assert df["MasterAlignmentID"].tolist() == [0, 1]
    assert df.loc[0, "rt"] == 1.5 and df.loc[0, "c_s_total"] == 34 and df.loc[0, "db_s_total"] == 1
    assert df.loc[0, "Ontology"] == "PC" and df.loc[0, "base_ontology"] == "PC"
    assert bool(df.loc[0, "inclusion"]) and df.loc[0, "exclusion_reason"] == ""  # 名前の "low score" は見ない
    assert bool(df.loc[0, "curation"])
    assert str(df["c_s_total"].dtype) == "Int64" and str(df["db_s_total"].dtype) == "Int64"


def test_generic_custom_and_case_insensitive_columns(tmp_path):
    t = gen(tmp_path, [["5.0", "36", "2", "PE"]], header=["rt (min)", "carbon", "dbl", "class"])
    df = generic_frame(t, GenericColumns("RT (min)", "Carbon", "DBL", "Class"), None, False)
    assert df.loc[0, "rt"] == 5.0 and df.loc[0, "db_s_total"] == 2 and df.loc[0, "Ontology"] == "PE"


def test_generic_missing_column_raises(tmp_path):
    t = gen(tmp_path, [["1.0", "34", "PC"]], header=["RT", "C", "Ontology"])
    with pytest.raises(TableError, match="DB"):
        generic_frame(t, GenericColumns(), None, False)


def test_generic_reasons_and_order(tmp_path):
    t = gen(tmp_path, [["a", "", "34", "1", "", "x"],
                       ["b", "1.0", "", "1", "PC", "x"],
                       ["c", "1.0", "34", "1", "unknown", "x"],
                       ["d", "1.0", "34", "1", "PC", "x"]])
    df = generic_frame(t, GenericColumns(), pd.Series([False, True, True, True]), keep_ontology_null=False)
    assert df["exclusion_reason"].tolist() == ["ontology_null;missing_value;use_column",
                                               "missing_value", "ontology_null", ""]
    assert df["inclusion"].tolist() == [False, False, False, True]
    assert df["curation"].tolist() == [False, False, False, True]


def test_generic_keep_ontology_null(tmp_path):
    t = gen(tmp_path, [["a", "1.0", "34", "1", "null", "x"]])
    df = generic_frame(t, GenericColumns(), None, keep_ontology_null=True)
    assert df.loc[0, "exclusion_reason"] == "" and bool(df.loc[0, "inclusion"])


@pytest.mark.parametrize("col,value", [("RT", "abc"), ("RT", "NaN"), ("RT", "inf"), ("C", "36.5"),
                                       ("DB", "x"), ("C", "nan")])
def test_generic_bad_value_raises(tmp_path, col, value):
    row = {"ID": "a", "RT": "1.0", "C": "34", "DB": "1", "Ontology": "PC", "Name": "x"}
    row[col] = value
    t = gen(tmp_path, [list(row.values())])
    with pytest.raises(TableError) as e:
        generic_frame(t, GenericColumns(), None, False)
    msg = str(e.value)
    assert f"'{col}'" in msg and value in msg and "1 件目" in msg
    if value.lower() in ("nan", "inf"):
        assert "空欄" in msg


def test_generic_integer_like_float_accepted(tmp_path):
    t = gen(tmp_path, [["a", "1.0", "36.0", "2.0", "PC", "x"]])
    df = generic_frame(t, GenericColumns(), None, False)
    assert df.loc[0, "c_s_total"] == 36 and df.loc[0, "db_s_total"] == 2


def test_generic_bad_values_ignored_on_excluded_rows(tmp_path):
    t = gen(tmp_path, [["a", "abc", "34", "1", "PC", "x"],   # use=False → チェックしない
                       ["b", "abc", "34", "1", "", "x"],     # Ontology が空 → チェックしない
                       ["c", "1.0", "34", "1", "PC", "x"]])
    df = generic_frame(t, GenericColumns(), pd.Series([False, True, True]), keep_ontology_null=True)
    assert pd.isna(df.loc[0, "rt"]) and pd.isna(df.loc[1, "rt"]) and df.loc[2, "rt"] == 1.0
    assert df["exclusion_reason"].tolist() == ["use_column", "", ""]


def test_write_table_keeps_input_and_appends(tmp_path):
    t = read_table(write_tsv(tmp_path / "h.txt", META4 + [HEADER, ROW_A, ROW_B]))
    result = run(preprocess(msdial_spots(t), FilterOptions()))
    out = tmp_path / "o" / "h_rtcurated.tsv"
    write_table(t, result, out)
    raw = out.read_bytes()
    assert raw.startswith(b"\xef\xbb\xbf")
    lines = raw.decode("utf-8-sig").splitlines()
    assert lines[0].split("\t") == HEADER + ADDED_COLUMNS  # "1" が 2 つのまま
    assert lines[1].split("\t")[:len(HEADER)] == ROW_A     # "13.793" などが元の文字列のまま
    assert lines[2].split("\t")[:len(HEADER)] == ROW_B
    assert len(lines) == 3


def test_added_columns_order():
    assert ADDED_COLUMNS == ["rt_group", "c_s_total", "db_s_total", "inclusion", "pair_evaluated", "curation",
                             "exclusion_reason", "n_miss", "n_correct", "n_total", "miss%", "exclusion_round"]


def test_write_table_rt_order_only(tmp_path):
    # PC 系列 6 件 + RT が大きく外れた 1 件(入力の 2 行目)
    rows = [[str(i), f"{rt:.3f}", "700", f"PC {c}:0", "[M+H]+", "PC", "1", "1"]
            for i, (c, rt) in enumerate([(30, 1), (48, 0.5), (32, 2), (34, 3), (36, 4), (38, 5), (40, 6)])]
    t = read_table(write_tsv(tmp_path / "h.txt", META4 + [HEADER] + rows))
    result = run(preprocess(msdial_spots(t), FilterOptions()))
    out = tmp_path / "ex.tsv"
    write_table(t, result, out, rt_order_only=True)
    lines = out.read_text(encoding="utf-8-sig").splitlines()
    assert len(lines) == 2 and lines[1].split("\t")[0] == "1"


def test_rt_group_is_corrected_ontology(tmp_path):
    t = read_table(write_tsv(tmp_path / "h.txt", META4 + [HEADER,
        ["0", "1.0", "800", "OxPC 34:2;O2", "[M+H]+", "OxPC", "1", "1"]]))
    result = run(preprocess(msdial_spots(t), FilterOptions()))
    out = tmp_path / "o.tsv"
    write_table(t, result, out)
    back = out.read_text(encoding="utf-8-sig").splitlines()
    header, row = back[0].split("\t"), back[1].split("\t")
    assert row[header.index("rt_group")] == "OxPC;O2"
    assert row[header.index("Ontology")] == "OxPC"  # 入力の列は元のまま


def test_check_added_columns_conflict(tmp_path):
    t = gen(tmp_path, [["a", "1.0", "34", "1", "PC", "x"]], header=["ID", "RT", "C", "DB", "Ontology", "curation"])
    with pytest.raises(TableError, match="curation"):
        check_added_columns(t)
    check_added_columns(gen(tmp_path, [["a", "1.0", "34", "1", "PC", "x"]]))  # ぶつからなければ何もしない


def test_quote_at_start_of_field_is_kept_verbatim(tmp_path):
    # MS-DIAL は値を引用符で囲まない。" で始まる値があっても、行をまとめずにそのまま読んで書き戻す
    rows = [ROW_A[:3] + ['"a'] + ROW_A[4:], ROW_B, ["2"] + ROW_B[1:3] + ['b"'] + ROW_B[4:], ROW_B[:0] + ["3"] + ROW_B[1:]]
    src = write_tsv(tmp_path / "h.txt", META4 + [HEADER] + rows)
    t = read_table(src)
    assert len(t.frame) == 4 and t.frame.iloc[0, 3] == '"a' and t.frame.iloc[2, 3] == 'b"'
    out = tmp_path / "o.tsv"
    write_table(t, run(preprocess(msdial_spots(t), FilterOptions())), out)
    lines = out.read_text(encoding="utf-8-sig").splitlines()
    assert len(lines) == 5
    assert [l.split("\t")[:len(HEADER)] for l in lines[1:]] == rows


def test_quotes_inside_field_are_written_verbatim(tmp_path):
    t = gen(tmp_path, [["a\"b", "1.0", "34", "1", "PC", '"quoted"']])
    out = tmp_path / "o.tsv"
    write_table(t, run(generic_frame(t, GenericColumns(), None, False)), out)
    row = out.read_text(encoding="utf-8-sig").splitlines()[1].split("\t")
    assert row[0] == 'a"b' and row[5] == '"quoted"'


def test_utf16_file_gets_clear_error(tmp_path):
    p = tmp_path / "x.txt"
    p.write_bytes("RT\tC\tDB\tOntology\r\n1\t34\t1\tPC\r\n".encode("utf-16"))  # BOM 付き UTF-16LE
    with pytest.raises(TableError, match="UTF-16"):
        read_table(p)


# ---- CSV と --header-row ----

def write_csv(path, text: str) -> Path:
    path.write_bytes(text.encode("utf-8"))  # 改行を OS の既定に変換せず、書いたとおりに保存する
    return path


def test_csv_msdial_detected(tmp_path):
    lines = META4 + [HEADER, ROW_A, ROW_B]
    t = read_table(write_csv(tmp_path / "m.csv", "\r\n".join(",".join(r) for r in lines) + "\r\n"))
    assert t.kind == "msdial_txt" and list(t.frame.columns) == HEADER and t.frame.iloc[1, 1] == "14.100"


def test_csv_quoted_field_keeps_comma(tmp_path):
    t = read_table(write_csv(tmp_path / "g.csv", 'RT,C,DB,Ontology,Note\n1.0,34,1,PC,"a, b"\n'))
    assert t.kind == "generic" and t.frame.iloc[0, 4] == "a, b"
    out = tmp_path / "o.tsv"
    write_table(t, run(generic_frame(t, GenericColumns(), None, False)), out)
    assert out.read_text(encoding="utf-8-sig").splitlines()[1].split("\t")[4] == "a, b"


@pytest.mark.parametrize("value", ['"x\ty"', '"x\ny"'])
def test_csv_value_with_tab_or_newline_raises(tmp_path, value):
    p = write_csv(tmp_path / "g.csv", f"RT,C,DB,Ontology,Note\n1.0,34,1,PC,{value}\n")
    with pytest.raises(TableError, match="Note"):
        read_table(p)


def test_csv_column_count_mismatch_raises(tmp_path):
    with pytest.raises(TableError, match="2 行目"):
        read_table(write_csv(tmp_path / "g.csv", "RT,C\n1,2,3\n"))


def test_tsv_is_not_split_on_commas(tmp_path):
    t = read_table(write_tsv(tmp_path / "g.tsv", [["RT", "Note"], ["1.0", "a, b"]]))
    assert list(t.frame.columns) == ["RT", "Note"] and t.frame.iloc[0, 1] == "a, b"


def test_header_row_selects_generic_header(tmp_path):
    p = write_tsv(tmp_path / "g.tsv", [["exported by tool X", ""], ["", ""], ["RT", "C"], ["1.0", "34"]])
    t = read_table(p, header_row=3)
    assert t.kind == "generic" and list(t.frame.columns) == ["RT", "C"] and len(t.frame) == 1


def test_header_row_selects_msdial_header_beyond_auto_scan(tmp_path):
    meta = [["", "", "", "", "", f"m{i}", "", ""] for i in range(25)]
    p = write_tsv(tmp_path / "h.txt", meta + [HEADER, ROW_A])
    assert read_table(p).kind == "generic"                    # 自動では 20 行目までしか探さない
    t = read_table(p, header_row=26)
    assert t.kind == "msdial_txt" and len(t.frame) == 1


def test_header_row_column_count_error_uses_file_line(tmp_path):
    p = write_tsv(tmp_path / "g.tsv", [["x"], ["RT", "C"], ["1.0"]])
    with pytest.raises(TableError, match="3 行目"):
        read_table(p, header_row=2)


@pytest.mark.parametrize("row", [0, -1, 4])
def test_header_row_out_of_range_raises(tmp_path, row):
    p = write_tsv(tmp_path / "g.tsv", [["RT", "C"], ["1.0", "34"], ["2.0", "36"]])
    with pytest.raises(TableError, match=str(row)):
        read_table(p, header_row=row)
