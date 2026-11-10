"""MS-DIAL5 のタグファイル (`<アライメント名>_tags.xml`) の読み書き。

RT 順序の判定結果をタグとして書き込み、MS-DIAL の画面で確認できるようにする。
書き出す形は MS-DIAL 5.5 が保存するファイルと同じ(UTF-8 BOM 付き、CRLF、インデント 2)。
"""

import shutil
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

import pandas as pd

CONFIRMED_TAG = 1
DEFAULT_RT_TAG = 5  # Overannotation を「RT 順序で外れた」の印として使う
DEFAULT_DEFINITIONS = [
    (1, "Confirmed"),
    (2, "Low quality spectrum"),
    (3, "Misannotation"),
    (4, "Coelution (mixed spectra)"),
    (5, "Overannotation"),
]


@dataclass
class TagFile:
    definitions: list[tuple[int, str]]
    peaks: dict[int, list[int]] = field(default_factory=dict)  # MasterAlignmentID → タグ ID


def read_tags(path: Path) -> TagFile:
    """既存のタグファイルを読む。無ければ MS-DIAL の既定の定義だけを返す。"""
    path = Path(path)
    if not path.exists():
        return TagFile(list(DEFAULT_DEFINITIONS), {})
    root = ET.parse(path).getroot()
    definitions = [(int(t.findtext("Id")), t.findtext("Label") or "")
                   for t in root.iterfind("Definitions/Tag")]
    peaks = {int(p.get("Id")): [int(t.text) for t in p.iterfind("Tag")]
             for p in root.iterfind("Peaks/Peak")}
    return TagFile(definitions or list(DEFAULT_DEFINITIONS), peaks)


def render_tags(tf: TagFile) -> bytes:
    lines = ['<?xml version="1.0" encoding="utf-8"?>', "<PeakSpotTags>", "  <Definitions>"]
    for tag_id, label in tf.definitions:
        lines += ["    <Tag>", f"      <Id>{tag_id}</Id>",
                  f"      <Label>{_escape(label)}</Label>", "    </Tag>"]
    lines.append("  </Definitions>")
    peaks = {i: t for i, t in tf.peaks.items() if t}
    if peaks:
        lines.append("  <Peaks>")
        for peak_id in sorted(peaks):
            lines.append(f'    <Peak Id="{peak_id}">')
            lines += [f"      <Tag>{t}</Tag>" for t in peaks[peak_id]]
            lines.append("    </Peak>")
        lines.append("  </Peaks>")
    else:
        lines.append("  <Peaks />")
    lines.append("</PeakSpotTags>")
    return "\r\n".join(lines).encode("utf-8-sig")


def _escape(text: str) -> str:
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def merge_tags(existing: dict[int, list[int]], df: pd.DataFrame, rt_tag: int = DEFAULT_RT_TAG,
               confirm: bool = False) -> tuple[dict[int, list[int]], list[int]]:
    """判定結果を既存のタグに反映する。

    - rt_tag は pairrt 専用とみなし、前回の分をすべて外してから今回の結果で付け直す。
    - Confirmed (1) は手作業のキュレーション用なので既定では触れない。
      confirm=True のときだけ、ペア評価で残ったものに付け足す(外すことはしない)。
    - それ以外のタグはそのまま残す。
    戻り値の 2 つ目は、Confirmed と rt_tag の両方が付いたスポットの ID。
    """
    merged = {i: [t for t in tags if t != rt_tag] for i, tags in existing.items()}
    evaluated = df["pair_evaluated"].astype(bool)
    kept = df.loc[evaluated & df["curation"].astype(bool), "MasterAlignmentID"]
    excluded = df.loc[df["exclusion_reason"] == "rt_order", "MasterAlignmentID"]
    additions = [(excluded, rt_tag)] + ([(kept, CONFIRMED_TAG)] if confirm else [])
    for ids, tag in additions:
        for i in ids.astype(int):
            tags = merged.setdefault(i, [])
            if tag not in tags:
                tags.append(tag)
    merged = {i: sorted(t) for i, t in merged.items() if t}
    both = sorted(i for i, t in merged.items() if CONFIRMED_TAG in t and rt_tag in t)
    return merged, both


def save_with_backup(target: Path, content: bytes, now: datetime) -> Path | None:
    """タグファイルを書く。既にあれば、古いものと新しいものをバックアップフォルダにも置く。

    戻り値はバックアップフォルダ(作らなかったときは None)。
    """
    target = Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    backup = None
    if target.exists():
        stamp = now.strftime("%Y%m%d_%H%M%S")
        backup = target.parent / f"rtcurate_backup_{stamp}"
        n = 2
        while backup.exists():
            backup = target.parent / f"rtcurate_backup_{stamp}_{n}"
            n += 1
        backup.mkdir()
        shutil.copy2(target, backup / f"old_{target.name}")
        (backup / f"new_{target.name}").write_bytes(content)
    target.write_bytes(content)
    return backup
