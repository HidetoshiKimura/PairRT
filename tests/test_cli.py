import os
from pathlib import Path

import pandas as pd
import pytest

from conftest import spots_df
from test_table import HEADER, META4, write_tsv
from rtcurate.cli import main
from rtcurate.tags import read_tags

STEM = "AlignmentResult_2026_10_20_11_07_17"


def read(path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", encoding="utf_8_sig")


def all_tags(path) -> list[int]:
    return [t for ts in read_tags(path).peaks.values() for t in ts]


def test_missing_file_returns_2(capsys, tmp_path):
    assert main([str(tmp_path / "nope.arf2")]) == 2
    assert "nope.arf2" in capsys.readouterr().err


def test_wrong_extension_returns_2(capsys, tmp_path):
    f = tmp_path / "a.csv"
    f.write_text("x")
    assert main([str(f)]) == 2
    assert ".arf2" in capsys.readouterr().err


def test_broken_arf2_returns_2(capsys, tmp_path):
    f = tmp_path / "broken.arf2"
    f.write_bytes(b"not msgpack")
    assert main([str(f)]) == 2
    assert "broken.arf2" in capsys.readouterr().err


def test_default_writes_only_tags_xml(sample_arf2, tmp_path, capsys):
    assert main([str(sample_arf2), "-o", str(tmp_path)]) == 0
    assert [p.name for p in tmp_path.iterdir()] == [f"{STEM}_tags.xml"]
    tags = all_tags(tmp_path / f"{STEM}_tags.xml")
    assert tags.count(5) == 106 and 1 not in tags  # 既定ではタグ 1 を付けない
    stdout = capsys.readouterr().out
    assert "検出スポット数: 7955" in stdout and f"{STEM}_tags.xml" in stdout


def test_tag_confirmed_option(sample_arf2, tmp_path):
    assert main([str(sample_arf2), "-o", str(tmp_path), "--tag-confirmed"]) == 0
    tags = all_tags(tmp_path / f"{STEM}_tags.xml")
    assert tags.count(5) == 106 and tags.count(1) == 1448  # ペア評価で残ったもの


def test_tsv_and_excluded_options(sample_arf2, tmp_path):
    assert main([str(sample_arf2), "-o", str(tmp_path), "--tsv", "--write-excluded"]) == 0
    assert len(read(tmp_path / f"{STEM}_rtcurated.tsv")) == 7955
    ex = read(tmp_path / f"{STEM}_rt_excluded.tsv")
    assert len(ex) == 106 and (ex["exclusion_reason"] == "rt_order").all()


def test_existing_tags_are_backed_up_and_merged(sample_arf2, tmp_path, capsys):
    target = tmp_path / f"{STEM}_tags.xml"
    # ID 0 = 手作業の Misannotation、ID 1〜9 = 前回の rtcurate の印(タグ 5)
    peaks = '<Peak Id="0"><Tag>3</Tag></Peak>' + "".join(
        f'<Peak Id="{i}"><Tag>5</Tag></Peak>' for i in range(1, 10))
    target.write_text('<?xml version="1.0" encoding="utf-8"?><PeakSpotTags><Definitions/>'
                      f"<Peaks>{peaks}</Peaks></PeakSpotTags>", encoding="utf-8")
    old = target.read_bytes()
    assert main([str(sample_arf2), "-o", str(tmp_path), "--tsv"]) == 0

    backups = [p for p in tmp_path.iterdir() if p.is_dir()]
    assert len(backups) == 1 and backups[0].name.startswith("rtcurate_backup_")
    assert (backups[0] / f"old_{STEM}_tags.xml").read_bytes() == old
    assert (backups[0] / f"new_{STEM}_tags.xml").read_bytes() == target.read_bytes()
    assert str(backups[0]) in capsys.readouterr().out

    new = read_tags(target).peaks
    assert new[0] == [3]  # 手作業のタグは残る
    rt_order = set(read(tmp_path / f"{STEM}_rtcurated.tsv")
                   .query("exclusion_reason == 'rt_order'")["MasterAlignmentID"])
    for i in range(1, 10):  # 前回の印は外れ、今回 rt_order のものだけに付く
        assert (5 in new.get(i, [])) == (i in rt_order)
    assert all_tags(target).count(5) == 106


def test_rt_tag_option(sample_arf2, tmp_path):
    assert main([str(sample_arf2), "-o", str(tmp_path), "--rt-tag", "3"]) == 0
    tags = all_tags(tmp_path / f"{STEM}_tags.xml")
    assert tags.count(3) == 106 and 5 not in tags


def test_threshold_option_changes_result(sample_arf2, tmp_path):
    assert main([str(sample_arf2), "-o", str(tmp_path), "--n-miss-threshold", "1000"]) == 0
    assert 5 not in all_tags(tmp_path / f"{STEM}_tags.xml")


def test_keep_and_exclude_options(sample_arf2, tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    main([str(sample_arf2), "-o", str(a), "--tsv"])
    main([str(sample_arf2), "-o", str(b), "--tsv", "--keep-low-score", "--exclude-ontology", "CL"])
    da, db = read(a / f"{STEM}_rtcurated.tsv"), read(b / f"{STEM}_rtcurated.tsv")
    assert db["inclusion"].sum() > da["inclusion"].sum()
    assert not db.loc[db["base_ontology"] == "CL", "pair_evaluated"].any()


def test_unwritable_output_returns_2(sample_arf2, tmp_path, capsys):
    # 出力先に同名のフォルダがあると書き込めない(Excel で開いたままの TSV と同じく OSError になる)
    (tmp_path / f"{STEM}_rtcurated.tsv").mkdir()
    assert main([str(sample_arf2), "-o", str(tmp_path), "--tsv"]) == 2
    assert f"{STEM}_rtcurated.tsv" in capsys.readouterr().err


def test_broken_existing_tags_returns_2_and_keeps_file(sample_arf2, tmp_path, capsys):
    target = tmp_path / f"{STEM}_tags.xml"
    target.write_text("<PeakSpotTags><Peaks>", encoding="utf-8")
    assert main([str(sample_arf2), "-o", str(tmp_path)]) == 2
    assert target.read_text(encoding="utf-8") == "<PeakSpotTags><Peaks>"
    assert [p.name for p in tmp_path.iterdir()] == [f"{STEM}_tags.xml"]
    assert f"{STEM}_tags.xml" in capsys.readouterr().err


def _cp932_stream():
    import io
    return io.TextIOWrapper(io.BytesIO(), encoding="cp932")


def test_git_bash_output_is_utf8(monkeypatch):
    from rtcurate.cli import _utf8_for_git_bash
    out, err = _cp932_stream(), _cp932_stream()
    monkeypatch.setattr("sys.stdout", out)
    monkeypatch.setattr("sys.stderr", err)
    monkeypatch.setenv("MSYSTEM", "MINGW64")
    _utf8_for_git_bash()
    assert out.encoding == "utf-8" and err.encoding == "utf-8"


def test_other_shells_keep_encoding(monkeypatch):
    from rtcurate.cli import _utf8_for_git_bash
    out = _cp932_stream()
    monkeypatch.setattr("sys.stdout", out)
    monkeypatch.delenv("MSYSTEM", raising=False)
    _utf8_for_git_bash()
    assert out.encoding == "cp932"


# 各フィルタ / Ontology 補正に 1 行ずつ当たるデータ(ID, 名前, Ontology)
WIRING_ROWS = [
    (1, "low score: PC 30:0", "PC"),
    (2, "no MS2: PC 30:0", "PC"),
    (3, "w/o MS2: PC 30:0", "PC"),
    (4, "RIKEN N-VS1 ID-1", "PC"),
    (5, "PC 30:0", ""),
    (6, "OxPC 34:2;O2", "OxPC"),
    (7, "PE P-36:2", "EtherPE"),
    (8, "PC O-34:1", "EtherPC"),
]


def run_wiring(monkeypatch, tmp_path, *options) -> pd.DataFrame:
    """.arf2 の読み込みを WIRING_ROWS に差し替えて実行し、TSV を ID で引ける形で返す。"""
    monkeypatch.setattr("rtcurate.cli.load_spots", lambda _: spots_df(WIRING_ROWS))
    src = tmp_path / "dummy.arf2"
    src.write_bytes(b"")
    assert main([str(src), "-o", str(tmp_path), "--tsv", *options]) == 0
    return read(tmp_path / "dummy_rtcurated.tsv").set_index("MasterAlignmentID")


@pytest.mark.parametrize("option,row", [
    ("--keep-low-score", 1),
    ("--keep-no-ms2", 2),
    ("--keep-wo-ms2", 3),
    ("--keep-riken", 4),
    ("--keep-ontology-null", 5),
])
def test_keep_option_disables_only_its_filter(monkeypatch, tmp_path, option, row):
    df = run_wiring(monkeypatch, tmp_path, option)
    filtered = {1, 2, 3, 4, 5}
    kept = set(df.index[df["inclusion"]]) & filtered
    assert kept == {row}  # 指定したフィルタの行だけが残る


def test_default_applies_all_filters(monkeypatch, tmp_path):
    df = run_wiring(monkeypatch, tmp_path)
    assert not df.loc[[1, 2, 3, 4, 5], "inclusion"].any()
    assert df.loc[6, "Ontology"] == "OxPC;O2"  # 酸素数タグは既定で付く
    assert df.loc[7, "Ontology"] == "EtherPE_P" and df.loc[8, "Ontology"] == "EtherPC_O"


def test_all_keep_options_together(monkeypatch, tmp_path):
    df = run_wiring(monkeypatch, tmp_path, "--keep-low-score", "--keep-no-ms2", "--keep-wo-ms2",
                    "--keep-riken", "--keep-ontology-null")
    assert df["inclusion"].all()


def test_ontology_options(monkeypatch, tmp_path):
    df = run_wiring(monkeypatch, tmp_path, "--no-oxygen-tag", "--no-ether-split")
    assert df.loc[6, "Ontology"] == "OxPC"
    assert df.loc[7, "Ontology"] == "EtherPE"
    assert df.loc[8, "Ontology"] == "EtherPC"


@pytest.mark.parametrize("option", ["--p-tag", "--o-tag"])
def test_removed_tag_options_are_rejected(option, tmp_path):
    with pytest.raises(SystemExit) as e:
        main([str(tmp_path / "a.arf2"), option])
    assert e.value.code == 2


def test_sample_with_all_filters_off(sample_arf2, tmp_path):
    # 実データで前処理フィルタをすべて外しても最後まで動き、全スポットが残る
    assert main([str(sample_arf2), "-o", str(tmp_path), "--tsv", "--keep-low-score", "--keep-no-ms2",
                 "--keep-wo-ms2", "--keep-riken", "--keep-ontology-null"]) == 0
    df = read(tmp_path / f"{STEM}_rtcurated.tsv")
    assert len(df) == 7955 and df["inclusion"].all()
    assert df["pair_evaluated"].sum() > 1554  # 既定(フィルタあり)のペア評価対象より多い


# ---- MS-DIAL の txt と汎用 TSV ----

def msdial_txt(tmp_path) -> Path:
    """PC 系列 8 件 + RT が外れた 1 件(ID 8)+ low score 1 件(ID 9)の MS-DIAL の txt。"""
    rows = [[str(i), f"{3.0 + i:.3f}", "700", f"PC {30 + 2 * i}:0", "[M+H]+", "PC", "1", "1"] for i in range(8)]
    rows.append(["8", "1.000", "700", "PC 48:0", "[M+H]+", "PC", "1", "1"])
    rows.append(["9", "5.000", "700", "low score: PE 34:1", "[M+H]+", "PE", "1", "1"])
    return write_tsv(tmp_path / "Height_1.txt", META4 + [HEADER] + rows)


def generic_tsv(tmp_path, keep_last=False) -> Path:
    """PC 5 件は順序どおり。use=FALSE の 6 件目 (30:0, RT 20) を入れると矛盾が生じる。"""
    rows = [[f"s{i}", f"{rt}", str(c), "0", "PC", "TRUE"]
            for i, (c, rt) in enumerate([(32, 5), (34, 6), (36, 7), (38, 8), (40, 9)])]
    rows.append(["s5", "20", "30", "0", "PC", "TRUE" if keep_last else "FALSE"])
    return write_tsv(tmp_path / "other.tsv", [["ID", "RT", "C", "DB", "Ontology", "keep"]] + rows)


def test_msdial_txt_end_to_end(tmp_path, capsys):
    src = msdial_txt(tmp_path)
    assert main([str(src)]) == 0
    out = read(tmp_path / "Height_1_rtcurated.tsv")
    assert len(out) == 10 and out["Alignment ID"].tolist() == list(range(10))  # 行は入力の順
    reason = dict(zip(out["Alignment ID"], out["exclusion_reason"].fillna("")))
    assert reason[8] == "rt_order" and reason[9] == "low_score"
    assert [r for i, r in reason.items() if i < 8] == [""] * 8
    assert not list(tmp_path.glob("*_tags.xml"))
    assert "検出スポット数: 10" in capsys.readouterr().out


def test_msdial_txt_write_excluded(tmp_path):
    assert main([str(msdial_txt(tmp_path)), "-o", str(tmp_path / "o"), "--write-excluded"]) == 0
    ex = read(tmp_path / "o" / "Height_1_rt_excluded.tsv")
    assert ex["Alignment ID"].tolist() == [8]


def test_msdial_txt_use_column(tmp_path):
    src = tmp_path / "h.txt"
    lines = msdial_txt(tmp_path).read_text(encoding="utf-8").splitlines()
    head = len(META4)
    lines = [l + "\t" for l in lines[:head]] + [lines[head] + "\tuse"] + \
            [l + ("\tFALSE" if l.startswith("8\t") else "\tTRUE") for l in lines[head + 1:]]
    src.write_text("\n".join(lines) + "\n", encoding="utf-8")
    assert main([str(src), "--use-column", "USE"]) == 0
    out = read(tmp_path / "h_rtcurated.tsv").set_index("Alignment ID")
    assert out.loc[8, "exclusion_reason"] == "use_column" and "rt_order" not in set(out["exclusion_reason"].fillna(""))


def test_generic_end_to_end_with_use_column(tmp_path):
    src = generic_tsv(tmp_path)
    assert main([str(src)]) == 0  # use 列を指定しなければ 6 件目も比べる → 1 件 rt_order
    assert (read(tmp_path / "other_rtcurated.tsv")["exclusion_reason"] == "rt_order").sum() == 1
    assert main([str(src), "--use-column", "keep"]) == 0
    out = read(tmp_path / "other_rtcurated.tsv")
    assert (out["exclusion_reason"] == "rt_order").sum() == 0
    assert out.loc[5, "exclusion_reason"] == "use_column" and out.loc[5, "ID"] == "s5"


def test_generic_custom_columns(tmp_path):
    src = write_tsv(tmp_path / "x.tsv", [["rt (min)", "Carbon", "Double", "Class"], ["1.0", "34", "1", "PC"]])
    assert main([str(src), "--rt-column", "RT (min)", "--c-column", "carbon", "--db-column", "double",
                 "--ontology-column", "class"]) == 0
    assert read(tmp_path / "x_rtcurated.tsv").loc[0, "rt_group"] == "PC"


def test_generic_prints_row_label(tmp_path, capsys):
    assert main([str(generic_tsv(tmp_path))]) == 0
    assert "行数: 6" in capsys.readouterr().out


def test_generic_ignored_options_warn(tmp_path, capsys):
    assert main([str(generic_tsv(tmp_path)), "--keep-low-score", "--no-ether-split",
                 "--keep-ontology-null"]) == 0
    err = capsys.readouterr().err
    assert "効果がありません" in err and "--keep-low-score" in err and "--no-ether-split" in err
    assert "--keep-ontology-null" not in err


@pytest.mark.parametrize("kind,options", [
    ("txt", ["--rt-tag", "5"]), ("txt", ["--tag-confirmed"]),
    ("tsv", ["--rt-tag", "3"]), ("tsv", ["--tag-confirmed"]),
    ("txt", ["--rt-column", "RT"]), ("arf2", ["--rt-column", "RT"]),
    ("arf2", ["--ontology-column", "X"]), ("arf2", ["--use-column", "keep"]),
])
def test_option_input_mismatch_returns_2(kind, options, tmp_path, monkeypatch, capsys):
    if kind == "arf2":
        monkeypatch.setattr("rtcurate.cli.load_spots", lambda _: spots_df(WIRING_ROWS))
        src = tmp_path / "dummy.arf2"
        src.write_bytes(b"")
    else:
        src = msdial_txt(tmp_path) if kind == "txt" else generic_tsv(tmp_path)
    assert main([str(src), "-o", str(tmp_path / "o"), *options]) == 2
    assert options[0] in capsys.readouterr().err
    assert not (tmp_path / "o").exists()


def test_table_bad_value_returns_2(tmp_path, capsys):
    src = write_tsv(tmp_path / "x.tsv", [["RT", "C", "DB", "Ontology"], ["abc", "34", "1", "PC"]])
    assert main([str(src)]) == 2
    assert "'RT'" in capsys.readouterr().err
    assert not (tmp_path / "x_rtcurated.tsv").exists()


def test_table_added_column_clash_returns_2(tmp_path, capsys):
    src = write_tsv(tmp_path / "x.tsv", [["RT", "C", "DB", "Ontology", "curation"], ["1", "34", "1", "PC", "a"]])
    assert main([str(src)]) == 2
    assert "curation" in capsys.readouterr().err


def test_txt_matches_arf2_on_real_files(tmp_path):
    # 同じアライメントの txt と arf2 を、環境変数 RTCURATE_PAIR_TXT / RTCURATE_PAIR_ARF2 で渡す
    txt, arf2 = os.environ.get("RTCURATE_PAIR_TXT"), os.environ.get("RTCURATE_PAIR_ARF2")
    if not (txt and arf2 and Path(txt).is_file() and Path(arf2).is_file()):
        pytest.skip("環境変数 RTCURATE_PAIR_TXT と RTCURATE_PAIR_ARF2 で実ファイルを指定してください")
    assert main([txt, "-o", str(tmp_path)]) == 0
    assert main([arf2, "-o", str(tmp_path), "--tsv"]) == 0
    t = read(tmp_path / f"{Path(txt).stem}_rtcurated.tsv")
    a = read(tmp_path / f"{Path(arf2).stem}_rtcurated.tsv")
    assert set(t.loc[t["exclusion_reason"] == "rt_order", "Alignment ID"]) == \
        set(a.loc[a["exclusion_reason"] == "rt_order", "MasterAlignmentID"])
