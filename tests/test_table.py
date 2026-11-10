from pathlib import Path

import pandas as pd
import pytest

from rtcurate.table import TableError, find_column, read_table

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
