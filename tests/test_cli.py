import pandas as pd

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
    f = tmp_path / "a.txt"
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
