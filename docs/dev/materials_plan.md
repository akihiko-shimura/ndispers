# 材料追加計画

v0.10 で群クラスが汎用化されたので、既存の点群（3m, 32, 4̄2m, 4mm, mm2）に属する結晶は
「Sellmeier + 参照 d 値 + 出典」の 1 ファイルで追加できる。優先順位は需要 × 実装コスト。

## 原則

- 係数は手元の PDF か Crossref で照合できる文献からのみ転写する。記憶から書かない。
- 追加ごとに独立の再抽出（別エージェント or 別日の自分）で係数を突き合わせる。
  このリポジトリの典型バグは物理ではなく転写ミス。
- 温度項のない Sellmeier は `dndT = 0` とし、docstring に「T は受けるが無視」と明記
  （α-BBO, calcite と同じ扱い）。
- 非線形係数の Miller スケーリングの実験的裏付けは結晶ごとに `_d_note` に書く
  （支持: BBO, KDP, LiNbO₃ — Alford & Smith 2001; 粗い: KTP — Shoji 1997; 未検証: それ以外）。

## 第 1 段（実装済み, v0.11）

| 材料 | クラス | 群 | Sellmeier の出典 | d の出典 | 備考 |
|---|---|---|---|---|---|
| サファイア | `Sapphire` | 3̄m（中心対称） | Tropf 1995 Table 22 (Malitson & Dodge) | — | 0.2–5.5 µm, T 項なし |
| DKDP (KD\*P) | `DKDP` | 4̄2m | Ghosh 1992 Table I/II（KDP と同形） | d₃₆ = 0.37 (Roberts 1992) | dn/dT あり |
| α-quartz | `Quartz` | 32 | Tropf 1995 Table 22 (Radhakrishnan 1951) | d₁₁ = 0.30 (Shoji 1997) | 絶対スケールの基準、正単軸、旋光性は未モデル |
| 5%MgO:CLN 両波 | `MgOLN_Zelmon1997` | 3m | Zelmon 1997 Table 2 (21 °C) | d₃₃ 25.0, d₃₁ 4.4 (Shoji 1997), d₂₂ 2.1 (Roberts 1992, 逆符号) | LN の複屈折 d_eff が初めて計算可能 |

## 第 2 段（候補、着手順）

| 材料 | 群 | コスト | 理由・出典の当て |
|---|---|---|---|
| ZGP, AgGaS₂, AgGaSe₂ | 4̄2m | 低 | 中赤外 OPO。新規コードゼロ。Sellmeier: Bhar / Zondy; d₃₆ は Roberts 1992 表 |
| CBO (CsB₃O₅) | 222 | 中 | 355 nm THG。群クラス `Biax_222`（d₁₄ のみ）が要る。CBO フォルダの文献 |
| KTA, RTP | mm2 | 低 | KTP と同じ群。mid-IR OPO・EO |
| BiBO (BiB₃O₆) | 2 | 中〜高 | Ti:Sa SHG/OPA で人気。単斜晶で誘電主軸が結晶軸と一致せず、テンソル 8 成分。理論ノートに節を足してから |
| SBO (SrB₄O₇) | mm2 | 中 | DUV 透過材料として。複屈折が小さく複屈折位相整合は不可、と明記する前提 |
| YAG, MgF₂, BK7 | — | 低 | 分散のみ。RII が既に網羅しており優先度は低い |

## 機能（材料より優先度が高い可能性）

- **QPM 周期計算** (PPLN/PPKTP): Δk − 2π/Λ = 0 を解く。SLN/MgOLN の d₃₃ eee が v0.10 で評価可能になったので、実装は小さい。

## 除外

- 本計画は公開リポジトリの文書なので、未発表・検討中の候補材料は載せない。
