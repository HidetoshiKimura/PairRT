import pandas as pd

from rtcurate.cli import main


def read(path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", encoding="utf_8_sig")


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


def test_sample_end_to_end(sample_arf2, tmp_path, capsys):
    out = tmp_path / "r.tsv"
    assert main([str(sample_arf2), "-o", str(out), "--write-excluded"]) == 0
    df = read(out)
    assert len(df) == 7955
    ex = read(tmp_path / "AlignmentResult_2026_10_20_11_07_17_rt_excluded.tsv")
    assert len(ex) > 0 and (ex["exclusion_reason"] == "rt_order").all()
    stdout = capsys.readouterr().out
    assert "検出スポット数: 7955" in stdout and str(out) in stdout


def test_excluded_not_written_by_default(sample_arf2, tmp_path):
    main([str(sample_arf2), "-o", str(tmp_path / "r.tsv")])
    assert [p.name for p in tmp_path.iterdir()] == ["r.tsv"]


def test_threshold_option_changes_result(sample_arf2, tmp_path):
    b = tmp_path / "b.tsv"
    assert main([str(sample_arf2), "-o", str(b), "--n-miss-threshold", "1000"]) == 0
    assert (read(b)["exclusion_reason"] == "rt_order").sum() == 0


def test_keep_and_exclude_options(sample_arf2, tmp_path):
    a, b = tmp_path / "a.tsv", tmp_path / "b.tsv"
    main([str(sample_arf2), "-o", str(a)])
    main([str(sample_arf2), "-o", str(b), "--keep-low-score", "--exclude-ontology", "CL"])
    da, db = read(a), read(b)
    assert db["inclusion"].sum() > da["inclusion"].sum()
    assert not db.loc[db["base_ontology"] == "CL", "pair_evaluated"].any()


def test_unwritable_output_returns_2(sample_arf2, tmp_path, capsys):
    # 既存フォルダを出力先にすると書き込めない(Excel で開いたままの TSV と同じく OSError になる)
    out = tmp_path / "is_a_dir"
    out.mkdir()
    assert main([str(sample_arf2), "-o", str(out)]) == 2
    assert str(out) in capsys.readouterr().err
