# rtcurate

同じ脂質クラス(Ontology)の中で「炭素数・二重結合数から予想される溶出順」と実際の RT が
矛盾するアノテーションに印を付ける CLI。判定は RT 順序だけで行う。
`RerativeRT_auto_checher4.4` の前処理とペア評価を移植したもの。

入力は 3 種類(拡張子と中身から自動で見分ける):

| 入力 | 見分け方 | 出力 |
|---|---|---|
| MS-DIAL5 の AlignmentResult(`.arf2`) | 拡張子 `.arf2` | MS-DIAL のタグファイル `_tags.xml` |
| MS-DIAL のアライメント結果の txt(`Height_1_….txt` など) | `.txt` / `.tsv` / `.csv` で、先頭 20 行以内に `Alignment ID` の列名の行がある | 判定結果の列を足した表 `_rtcurated.tsv` |
| 汎用 TSV(RT・C・DB・Ontology の列を持つ表。MS-DIAL 以外のツールの結果など) | 上のどちらでもない `.txt` / `.tsv` / `.csv` | 判定結果の列を足した表 `_rtcurated.tsv` |

`.txt` / `.tsv` はタブ区切り、`.csv` はカンマ区切りとして読む。

## セットアップ

```
uv sync
```

## 使い方

```
uv run rtcurate AlignmentResult_xxx.arf2 [options]
```

1 回の実行で 1 ファイル(1 極性)を処理する。Neg と Pos はそれぞれ実行する。

既定では MS-DIAL のタグファイル `<入力名>_tags.xml` を .arf2 と同じフォルダに 1 つだけ書き出す。
MS-DIAL でそのアライメントを開くと、判定結果がタグとして表示される。

- RT 順序で外れたスポット → **Overannotation(5)**(`--rt-tag` で変更可)
- それ以外 → タグなし。Confirmed(1)は手作業のキュレーション用として触れない
  (`--tag-confirmed` を付けると、ペア評価で残ったスポットに Confirmed(1)を付け足す)

既存のタグファイルがある場合:

- 手作業で付けたタグは残す。ただし rt タグ(既定 5)は rtcurate 専用として、前回の分を外して付け直す。
- Confirmed(1)は外さない。手作業で Confirmed にしたスポットが RT 順序で外れた場合は 1 と 5 の両方が付くので、件数と ID を表示する。
- 書き換える前に `rtcurate_backup_<日時>/` を作り、古いファイル(`old_…`)と新しいファイル(`new_…`)を置く。

| オプション | 既定値 | 意味 |
|---|---|---|
| `-o, --out-dir DIR` | 入力と同じフォルダ | 出力先フォルダ |
| `--tsv` | オフ | 判定結果の表 `<入力名>_rtcurated.tsv` も出力する |
| `--write-excluded` | オフ | 自動で外したものの一覧を `<入力名>_rt_excluded.tsv` も出力する |
| `--overwrite` | オフ | 出力の表(`_rtcurated.tsv` / `_rt_excluded.tsv`)がすでにあれば上書きする。既定では何も書かずにエラーで止める(`_tags.xml` は常に古いものを退避してから書くので対象外) |
| `--rt-tag {2,3,4,5}` | 5 | RT 順序で外れたスポットに付けるタグ ID |
| `--tag-confirmed` | オフ | ペア評価で残ったスポットに Confirmed(1)を付け足す |
| `--n-miss-threshold INT` | 4 | 除外候補になる n_miss の下限(以上) |
| `--miss-percent-threshold FLOAT` | 10 | 除外候補になる miss% の下限(以上) |
| `--keep-low-score` / `--keep-no-ms2` / `--keep-wo-ms2` / `--keep-riken` / `--keep-ontology-null` | すべてオフ(=除外する) | それぞれの前処理フィルタを無効にする |
| `--no-oxygen-tag` | オフ(=付ける) | Ontology への酸素数タグ(`;O2` など)付与を無効にする |
| `--no-ether-split` | オフ(=分ける) | Ether 系の O- / P- を分けずに判定する。既定では `EtherPE_O` / `EtherPE_P` のように別の Ontology にする(P- は二重結合を 1 つ少なく書くため) |
| `--exclude-ontology NAME` | なし | 指定した Ontology をペア評価から外す(複数回指定可、`OxPC` で `OxPC;O2` なども外れる) |
| `--rt-column` / `--c-column` / `--db-column` / `--ontology-column NAME` | `RT` / `C` / `DB` / `Ontology` | 汎用 TSV の列名(大文字小文字は区別しない) |
| `--use-column NAME` | なし | この列が `FALSE` / `0` / `no` の行を判定から完全に外す(txt / 汎用 TSV のみ) |
| `--header-row N` | 自動で探す | 列名の行をファイルの N 行目(1 始まり、Excel の行番号と同じ)に指定する(txt / TSV / CSV のみ) |

## MS-DIAL の txt と汎用 TSV

```
uv run rtcurate Height_1_2026_11_09_16_43_20.txt
uv run rtcurate other_tool.tsv --rt-column "RT (min)" --c-column Carbon --db-column DoubleBond --ontology-column Class
uv run rtcurate other_tool.tsv --use-column keep
```

- 出力は `<入力名>_rtcurated.tsv`(常に出す)。入力の全列を元の順番・元の文字列のまま左に置き、
  判定結果の列(`rt_group`、`c_s_total`、`db_s_total`、`inclusion`、`pair_evaluated`、`curation`、
  `exclusion_reason`、`n_miss`、`n_correct`、`n_total`、`miss%`、`exclusion_round`)を右に足す。
  行は入力と同じ順番・同じ件数。`rt_group` は実際に比べたグループ(補正後の Ontology)。
- タグファイルは作らない。`--rt-tag` / `--tag-confirmed` を付けるとエラー。
- MS-DIAL の txt は arf2 と同じ前処理(除外フィルタ、名前からの炭素数・二重結合数、Ontology の補正)をかける。
  txt の RT は小数 3 桁に丸められているので、RT が同じ値になるペアがあると arf2 と結果が変わることがある。
- 汎用 TSV では名前を使った処理は何もしない。C・DB・Ontology は列の値をそのまま使う。
  - RT・C・DB のどれかが空欄の行は `missing_value` として判定から外す。
  - 数値として読めない値(`NaN`、`inf` も含む)や、C・DB が整数でない値はエラー。
    ただし `--use-column` が偽の行や Ontology が空の行は確かめない。
  - n_miss と miss% が同じ候補は、行番号の小さい方から外す。
  - `--keep-low-score` などの名前を使うオプションは効果がない(付けると注意を表示する)。`--keep-ontology-null` は効く。
- MS-DIAL の txt で `Alignment ID` が重複している場合(Pos と Neg をまとめた表など)は、同点のときに外す順番に
  ID の代わりに行番号を使う(注意を表示する)。Pos と Neg の行は、同じ Ontology なら同じグループで比べられる。
- 列名の行は自動で探す(MS-DIAL の形なら先頭 20 行以内の `Alignment ID` の行、それ以外は 1 行目)。
  表の上に説明の行があるファイルは `--header-row N` で指定する。指定した行より上は出力に書かない。
- `.csv` では `"` で囲まれた値の中のカンマも値の一部として読む。`.txt` / `.tsv` では `"` も普通の文字として読む。
  出力はタブ区切りなので、タブや改行を含む値があるとエラー。
- 文字コードは UTF-8(BOM の有無どちらも)、読めなければ cp932。UTF-16 は非対応。

## TSV の読み方(`--tsv`)

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

実ファイルを使うテストは、次の環境変数でファイルを指定したときだけ動く(無ければスキップする):

- `RTCURATE_TAGGED_XML`: MS-DIAL がタグを書いた `_tags.xml`
- `RTCURATE_PAIR_TXT` / `RTCURATE_PAIR_ARF2`: 同じアライメントの txt と `.arf2`(自動除外が一致するかを確かめる)
