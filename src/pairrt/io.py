"""MS-DIAL5 の AlignmentResult (.arf2) を読み込む。"""

import os

import pandas as pd
from msdialio import read_arf2

# msdialio のキー → pairrt の列名(この順で出力する)
_COLUMNS = {
    "MasterAlignmentID": "MasterAlignmentID",
    "Name": "name",
    "MassCenter": "mz",
    "RT": "rt",
    "Ontology": "Ontology",
    "AdductType": "Adduct_type",
}


def load_spots(path: str | os.PathLike) -> pd.DataFrame:
    """`.arf2` のアライメントスポットを、1 行 = 1 スポットの DataFrame で返す。"""
    records = read_arf2(path, slim=False)
    df = pd.DataFrame(
        [{new: rec[old] for old, new in _COLUMNS.items()} for rec in records],
        columns=list(_COLUMNS.values()),
    )
    df["MasterAlignmentID"] = df["MasterAlignmentID"].astype(int)
    df["mz"] = df["mz"].astype(float)
    df["rt"] = df["rt"].astype(float)

    duplicated = df.loc[df["MasterAlignmentID"].duplicated(), "MasterAlignmentID"]
    if not duplicated.empty:
        raise ValueError(
            f"duplicate MasterAlignmentID in {os.fspath(path)}: {sorted(set(duplicated))[:10]} "
            "(LC-IMS のドリフトスポットを含む arf2 には対応していません)"
        )
    return df
