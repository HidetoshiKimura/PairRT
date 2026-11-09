import os
from datetime import datetime
from pathlib import Path

import pandas as pd
import pytest

from rtcurate.tags import (
    DEFAULT_DEFINITIONS, TagFile, merge_tags, read_tags, render_tags, save_with_backup,
)

EMPTY_XML = (
    "\ufeff<?xml version=\"1.0\" encoding=\"utf-8\"?>\r\n"
    "<PeakSpotTags>\r\n"
    "  <Definitions>\r\n"
    "    <Tag>\r\n      <Id>1</Id>\r\n      <Label>Confirmed</Label>\r\n    </Tag>\r\n"
    "    <Tag>\r\n      <Id>2</Id>\r\n      <Label>Low quality spectrum</Label>\r\n    </Tag>\r\n"
    "    <Tag>\r\n      <Id>3</Id>\r\n      <Label>Misannotation</Label>\r\n    </Tag>\r\n"
    "    <Tag>\r\n      <Id>4</Id>\r\n      <Label>Coelution (mixed spectra)</Label>\r\n    </Tag>\r\n"
    "    <Tag>\r\n      <Id>5</Id>\r\n      <Label>Overannotation</Label>\r\n    </Tag>\r\n"
    "  </Definitions>\r\n"
    "  <Peaks />\r\n"
    "</PeakSpotTags>"
)

# MS-DIAL 5.5 が実際に書いたタグ入りのファイル(ID 71 に Misannotation、ID 20 に Confirmed)。
# 手元のファイルのパスを環境変数 RTCURATE_TAGGED_XML で渡す。無ければこのテストはスキップする。
REAL_TAGGED = Path(os.environ.get("RTCURATE_TAGGED_XML", ""))


def result_df(rows) -> pd.DataFrame:
    """(MasterAlignmentID, pair_evaluated, curation, exclusion_reason) から判定結果を作る。"""
    return pd.DataFrame(rows, columns=["MasterAlignmentID", "pair_evaluated", "curation",
                                       "exclusion_reason"])


def test_render_empty_matches_msdial_bytes():
    assert render_tags(TagFile(DEFAULT_DEFINITIONS, {})) == EMPTY_XML.encode("utf-8")


def test_read_render_roundtrip_of_real_file(tmp_path):
    if not os.environ.get("RTCURATE_TAGGED_XML") or not REAL_TAGGED.is_file():
        pytest.skip("環境変数 RTCURATE_TAGGED_XML でタグ入りの実ファイルを指定してください")
    original = REAL_TAGGED.read_bytes()
    tf = read_tags(REAL_TAGGED)
    assert tf.peaks[71] == [3] and tf.peaks[20] == [1]
    assert render_tags(tf) == original


def test_render_peaks_sorted_with_multiple_tags():
    xml = render_tags(TagFile(DEFAULT_DEFINITIONS, {7: [1, 5], 3: [5]})).decode("utf-8-sig")
    assert xml.index('<Peak Id="3">') < xml.index('<Peak Id="7">')
    assert '<Peak Id="7">\r\n      <Tag>1</Tag>\r\n      <Tag>5</Tag>\r\n    </Peak>' in xml


def test_read_missing_file_gives_defaults(tmp_path):
    tf = read_tags(tmp_path / "none_tags.xml")
    assert tf.definitions == DEFAULT_DEFINITIONS and tf.peaks == {}


def test_merge_tags_kept_and_excluded():
    df = result_df([
        (1, True, True, ""),           # 評価して残った → 1
        (2, True, False, "rt_order"),  # RT 順序で外れた → 5
        (3, False, True, ""),          # 評価していない → タグなし
        (4, False, False, "low_score"),
    ])
    merged, both = merge_tags({}, df, rt_tag=5, confirm=True)
    assert merged == {1: [1], 2: [5]} and both == []


def test_merge_default_does_not_touch_tag1():
    old = {1: [1], 2: [1], 5: [3]}  # 1, 2 は手作業の Confirmed
    df = result_df([(1, True, True, ""), (2, True, False, "rt_order"), (3, True, True, "")])
    merged, both = merge_tags(old, df, rt_tag=5)
    assert merged == {1: [1], 2: [1, 5], 5: [3]}  # 3 には 1 を付けず、手作業の 1 も外さない
    assert both == [2]


def test_merge_resets_old_rt_tag_but_keeps_manual_tags():
    old = {2: [5], 9: [5], 10: [3], 11: [1]}  # 9 は前回の rtcurate の印、10/11 は手作業
    df = result_df([(2, True, False, "rt_order"), (9, True, True, "")])
    merged, _ = merge_tags(old, df, rt_tag=5, confirm=True)
    assert merged == {2: [5], 9: [1], 10: [3], 11: [1]}


def test_merge_never_removes_tag1_and_reports_both():
    old = {2: [1]}  # 前回 Confirmed、今回は外れた
    df = result_df([(2, True, False, "rt_order")])
    merged, both = merge_tags(old, df, rt_tag=5, confirm=True)
    assert merged == {2: [1, 5]} and both == [2]


def test_merge_with_other_rt_tag():
    merged, _ = merge_tags({3: [4], 5: [5]}, result_df([(3, True, False, "rt_order")]), rt_tag=4)
    assert merged == {3: [4], 5: [5]}  # 4 は付け直し、手作業の 5 は残る


def test_save_without_existing_file_makes_no_backup(tmp_path):
    target = tmp_path / "A_tags.xml"
    backup = save_with_backup(target, b"new", datetime(2026, 10, 30, 12, 0, 0))
    assert backup is None
    assert target.read_bytes() == b"new"
    assert [p.name for p in tmp_path.iterdir()] == ["A_tags.xml"]


def test_save_with_existing_file_backs_up_old_and_new(tmp_path):
    target = tmp_path / "A_tags.xml"
    target.write_bytes(b"old")
    backup = save_with_backup(target, b"new", datetime(2026, 10, 30, 12, 0, 0))
    assert backup == tmp_path / "rtcurate_backup_20261030_120000"
    assert (backup / "old_A_tags.xml").read_bytes() == b"old"
    assert (backup / "new_A_tags.xml").read_bytes() == b"new"
    assert target.read_bytes() == b"new"


def test_save_backup_folder_name_collision(tmp_path):
    target = tmp_path / "A_tags.xml"
    now = datetime(2026, 10, 30, 12, 0, 0)
    target.write_bytes(b"v1")
    first = save_with_backup(target, b"v2", now)
    second = save_with_backup(target, b"v3", now)
    assert first != second and (second / "old_A_tags.xml").read_bytes() == b"v2"
