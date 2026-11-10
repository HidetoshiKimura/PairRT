import pytest

from pairrt.io import load_spots


def test_load_sample_columns_and_count(sample_arf2):
    df = load_spots(sample_arf2)
    assert list(df.columns) == ["MasterAlignmentID", "name", "mz", "rt", "Ontology", "Adduct_type"]
    assert len(df) == 7955
    assert df["MasterAlignmentID"].is_unique
    assert df["name"].str.startswith("low score:").sum() == 1849


def test_duplicate_ids_raise(monkeypatch):
    fake = [dict(MasterAlignmentID=1, Name="a", MassCenter=1.0, RT=1.0, Ontology="PC", AdductType="")] * 2
    monkeypatch.setattr("pairrt.io.read_arf2", lambda p, slim: fake)
    with pytest.raises(ValueError, match="duplicate"):
        load_spots("x.arf2")
