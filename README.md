# rtcurate

MS-DIAL5 の AlignmentResult(`.arf2`)を読み、同じ脂質クラス(Ontology)の中で
「炭素数・二重結合数から予想される溶出順」と実際の RT が矛盾するアノテーションに印を付ける CLI。
判定は RT 順序だけで行う。`RerativeRT_auto_checher4.4` の前処理とペア評価を移植したもの。

設計: [docs/2026-10-22-rtcurate-design.md](docs/2026-10-22-rtcurate-design.md)

## セットアップ

```
uv sync
```

## 使い方

```
uv run rtcurate AlignmentResult_xxx.arf2 [options]
```

1 回の実行で 1 ファイル(1 極性)を処理する。Neg と Pos はそれぞれ実行する。

| オプション | 既定値 | 意味 |
|---|---|---|
| `-o, --output PATH` | 入力と同じフォルダの `<入力名>_rtcurated.tsv` | 出力 TSV のパス |
| `--write-excluded` | オフ | 自動で外したものの一覧を `<入力名>_rt_excluded.tsv` として追加出力する(出力 TSV と同じフォルダ) |
| `--n-miss-threshold INT` | 4 | 除外候補になる n_miss の下限(以上) |
| `--miss-percent-threshold FLOAT` | 10 | 除外候補になる miss% の下限(以上) |
| `--keep-low-score` / `--keep-no-ms2` / `--keep-wo-ms2` / `--keep-riken` / `--keep-ontology-null` | すべてオフ(=除外する) | それぞれの前処理フィルタを無効にする |
| `--no-oxygen-tag` | オフ(=付ける) | Ontology への酸素数タグ(`;O2` など)付与を無効にする |
| `--p-tag` / `--o-tag` | オフ | 4.4 の `p_tag` / `o_tag` に相当 |
| `--exclude-ontology NAME` | なし | 指定した Ontology をペア評価から外す(複数回指定可、`OxPC` で `OxPC;O2` なども外れる) |

## 出力の読み方

1 行 = 1 アライメントスポット。主な列:

- `curation`: `True` = RT 順序の点で問題なし(残す)/ `False` = 外す
- `exclusion_reason`: 外した理由。`rt_order`(RT 順序)、`low_score` / `no_MS2` / `wo_MS2` / `RIKEN` / `ontology_null`(前処理)
- `pair_evaluated`: RT 順序の比較にかけたか。c/db が取れない名前などは `False` で、`curation=True` のまま残る
- `n_miss` / `n_correct` / `miss%`: 同じ Ontology 内で順序が矛盾した/合っていたペア数。自動除外した行は除外したときの値
- `exclusion_round`: その Ontology の中で何回目の反復で外したか

## テスト

```
uv run pytest
```
