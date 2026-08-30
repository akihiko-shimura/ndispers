# SNLO 結晶リスト（ndispers 実装状況つき）

SNLO v80（2025-01-22）が内蔵する全結晶。出典は Web 版 Qmix
(<https://smithjj.github.io>) に埋め込まれた `INLINE_CRYSTAL_INFO` を
そのまま抽出したもの（2026-08-30 取得）。組成・点群・透過域・光学クラスは
SNLO のデータそのままで、加筆していない（点群だけ SNLO の ASCII 表記
`4(bar)2m` / `6m2` を国際記号 4̄2m / 6̄m2 に直した）。

- 実体は **85 結晶**（+ DIY 枠 `ZZ_U` / `ZZ_B` の 2 つ）。
- `ndispers` 列: ✅ 実装済み / △ 部分的（別組成で代用）/ — 未実装。
  現状 **19 エントリが ✅、1 つが △、65 が未実装**（DIY 枠 2 つを除く）。
- 「主な用途」は SNLO のデータには含まれていない。一般的な文献知識から
  書いた要約で、**出典照合はしていない**。材料選定の当たりをつける用途に留める。
- 透過域は SNLO の Sellmeier 有効範囲（µm）であって、結晶の実透過域とは限らない。

## 対応表の読み方

SNLO は同一組成で Sellmeier が複数あるとき別エントリにする（KTP_F/H/K、KTA_1/2/3、
LNB_C/M/S など）。ndispers も同じ方針（`AGS_Kato1996` と `AGS_Takaoka1999` など）なので、
両者は素直に 1 対 1 で対応づけられる。ただし文献の選び方は一致していない。


## リン酸・ヒ酸塩系（KDP ファミリー, 4̄2m）

| 略称 | 組成 | 点群 | 光学クラス | 透過域 (µm) | 主な用途 | ndispers |
|---|---|---|---|---|---|---|
| ADP | NH_4H_2PO_4 | 4̄2m | 負一軸 | 0.184–1.5 | 水溶液成長の古典材料。可視 SHG、ポッケルスセル | — |
| ADA | NH_4H_2AsO_4 | 4̄2m | 負一軸 | 0.22–1.2 | ADP の As 置換体。1064 nm 付近で 90° 位相整合 SHG | — |
| DADP | ND_4D_2PO_4 | 4̄2m | 負一軸 | 0.22–1.7 | 重水素化 ADP。IR 吸収低減、電気光学変調 | — |
| DADA | ND_4D_2AsO_4 | 4̄2m | 負一軸 | 0.22–1.2 | 重水素化 ADA。低温 NCPM SHG | — |
| KDP | KH_2PO_4 | 4̄2m | 負一軸 | 0.177–1.7 | 大口径成長が容易。慣性核融合レーザーの 2ω/3ω 変換、ポッケルスセル | ✅ `crystals.KDP` |
| KDA | KH_2AsO_4 | 4̄2m | 負一軸 | 0.216–1.7 | KDP の As 置換体。1064 nm SHG、電気光学 | — |
| DKDP | KD_2PO_4 | 4̄2m | 負一軸 | 0.2–2 | KD*P。1064 nm の 3ω/4ω、大型ポッケルスセル。KDP より IR 吸収が小さい | ✅ `crystals.DKDP` |
| RDP | RbH_2PO_4 | 4̄2m | 負一軸 | 0.22–1.5 | 電気光学変調器（低半波電圧） | — |
| RDA | RbH_2AsO_4 | 4̄2m | 負一軸 | 0.26–1.46 | 1064 nm で 90° 位相整合 SHG、Q スイッチ | — |
| DRDP | RhD_2PO_4 | 4̄2m | 負一軸 | 0.22–1.5 | 重水素化 RDP。電気光学変調 | — |
| DRDA | RhD_2AsO_4 | 4̄2m | 負一軸 | 0.26–1.7 | 重水素化 RDA。NCPM SHG | — |
| CDA | CsH_2AsO_4 | 4̄2m | 負一軸 | 0.26–1.43 | 1064 nm 90° NCPM SHG の定番（温度整合） | — |
| DCDA | CsD_2AsO_4 | 4̄2m | 負一軸 | 0.27–1.66 | 重水素化 CDA。1064 nm NCPM SHG | — |

## ボレート系（UV／深紫外の主力）

| 略称 | 組成 | 点群 | 光学クラス | 透過域 (µm) | 主な用途 | ndispers |
|---|---|---|---|---|---|---|
| BBO | beta-BaB_2O_4 | 3m | 負一軸 | 0.185–2.6 | β-BBO。UV SHG/THG/FHG、fs OPA、ポッケルスセル。最も汎用 | ✅ `crystals.BetaBBO_{Eimerl1987,Ghosh1995,KK2010,Tamosauskas2018}` |
| LBO | LiB_3O_5 | mm2 | 二軸 | 0.16–2.6 | 1064 nm SHG/THG、高平均出力の緑色発生、NCPM OPO。損傷閾値が高い | ✅ `crystals.LBO_{Castech,Ghosh1995,KK1994,KK2018,Newlight}_{xy,yz,zx}` |
| CLBO | CsLiB_6O_10 | 4̄2m | 負一軸 | 0.18–2.75 | 266/213 nm の深紫外発生（4ω/5ω）。潮解性あり | ✅ `crystals.CLBO` |
| CBO | CsB_3O_5 | 222 | 二軸 | 0.17–3 | 1064 nm THG、深紫外。LBO 類似で複屈折が大きい | — |
| LB4 | Li_2B_4O_7 | 4mm | 負一軸 | 0.16–3.5 | 深紫外 SHG（~240 nm）、高損傷閾値、非潮解 | ✅ `crystals.LB4` |
| KBBF | KBe_2BO_3F_2 | 32 | 負一軸 | 0.147–3.5 | 177.3 nm の直接 SHG（唯一の実用材料）。層状剥離のためプリズム結合が必要 | ✅ `crystals.KBBF` |
| RBBF | RBe_2BO_3F_2 | 32 | 負一軸 | 0.165–3.5 | KBBF 類似の深紫外材料。成長性が改善 | ✅ `crystals.RBBF` |
| CBBF | CsBe_2BO_3F_2 | 32 | 負一軸 | 0.15–3.5 | KBBF 系。深紫外 SHG | — |
| KABO | K_2Al_2B_2O_7 | 32 | 負一軸 | 0.18–3.6 | 深紫外 SHG（~200 nm）、非潮解 | — |
| KBO | KB_5O_8:4H_2O | mm2 | 二軸 | 0.165–1.4 | 深紫外 SHG（~217 nm）。潮解性が強い | — |
| BIBO | BiB_3O_6 | 2 | 二軸 | 0.286–2.5 | 1064 nm SHG で d_eff 大。fs OPA/OPO、緑色発生 | ✅ `crystals.BiBO_Miyata2009_{xy,yz,zx}` |
| LCB | La_2CaB_10O_19 | 2 | 二軸 | 0.18–3.3 | 紫外 SHG、非潮解 | — |
| NLBO | Na_3La_9O_3(BO_3)_8 | 6̄m2 | 負一軸 | 0.2–3.5 | 紫外 SHG | — |
| BBPO | BaBPO_5 | 32 | 負一軸 | 0.18–3.3 | 紫外 SHG | — |
| BPO | BPO_4 | 4̄ | 正一軸 | 0.135–3.5 | 真空紫外まで透過。深紫外 SHG | — |
| LBGO | LaBGeO_5 lanthanum borogermanate | 3 | 正一軸 | 0.195–2.6 | 紫外 SHG、周期分極反転（QPM）候補 | — |
| GCOB | GdCa_4O(BO_3)_3 | m | 二軸 | 0.32–2.7 | 自己周波数変換レーザー母材、1064 nm SHG | — |
| YCOB | YCa_4O(BO_3)_3 | m | 二軸 | 0.22–3.5 | 大口径成長可。高平均出力 SHG、自己周波数変換 | — |
| TCOB | TbCa_4O(BO_3)_3 | m | 二軸 | 0.52–1.6 | COB 系。SHG、自己周波数変換 | — |

## ヨウ素酸塩・有機・その他の酸化物

| 略称 | 組成 | 点群 | 光学クラス | 透過域 (µm) | 主な用途 | ndispers |
|---|---|---|---|---|---|---|
| LIIO3 | LiIO_3 | 6 | 負一軸 | 0.3–6 | LiIO3。可視〜近赤外 SHG/OPO の古典材料。潮解性・低損傷閾値 | ✅ `crystals.LiIO3` |
| DLAP | L-Arginine phosphate | 2 | 二軸 | 0.23–1.2 | 有機系。可視 SHG、高損傷閾値 | — |
| LFM | LiCOOH-H_2O | mm2 | 二軸 | 0.23–1.2 | 有機系。紫外〜可視 SHG | — |
| GEO | GeO_2 alpha-Germanium dioxide | 32 | 正一軸 | 0.4–5.1 | α-GeO2。石英類似で複屈折が大きい。紫外〜中赤外の複屈折素子 | — |
| LGN | La_3 Ga_5.5 Nb_0.5 O_14 | 32 | 正一軸 | 0.43–6.8 | ランガサイト系。中赤外 OPO/DFG、圧電 | — |
| LGTO | La_3Ga_5.5Ta_0.5O_14 | 32 | 正一軸 | 0.5–6 | ランガテイト。中赤外、圧電 | — |
| CTW | Cs_2TeW_3O_12 | 6 | 負一軸 | 0.41–5.31 | テルル酸タングステン。可視〜中赤外 SHG | — |
| NTW | Na_2TeW_2O_9 | m | 二軸 | 0.45–5 | テルル酸タングステン。可視〜中赤外 SHG | — |

## KTP ファミリー（mm2 二軸）

| 略称 | 組成 | 点群 | 光学クラス | 透過域 (µm) | 主な用途 | ndispers |
|---|---|---|---|---|---|---|
| KTP_F | KTiOPO_4 | mm2 | 二軸 | 0.35–4.5 | フラックス成長 KTP。1064 nm SHG、OPO、電気光学。PPKTP 母材 | ✅ `crystals.KTP_{xy,yz,zx}` |
| KTP_H | KTiOPO_4 | mm2 | 二軸 | 0.35–4.5 | 水熱合成 KTP。灰色化（gray tracking）耐性が高い | ✅ `crystals.KTP_{xy,yz,zx}` |
| KTP_K | KTiOPO_4 | mm2 | 二軸 | 0.35–4.5 | フラックス成長 KTP（Kato の Sellmeier） | ✅ `crystals.KTP_{xy,yz,zx}` |
| KTA_1 | KTiOAsO_4 | mm2 | 二軸 | 0.35–4 | KTP より中赤外透過が広い。3-4 µm OPO | — |
| KTA_2 | KTiOAsO_4 | mm2 | 二軸 | 0.35–4 | KTA（別 Sellmeier） | — |
| KTA_3 | KTiOAsO_4 | mm2 | 二軸 | 0.35–4 | KTA（別 Sellmeier） | — |
| RTP | RbTiOPO_4 | mm2 | 二軸 | 0.35–4.5 | 高繰り返し電気光学 Q スイッチ／ポッケルスセル。低イオン伝導 | — |
| RTA | RbTiOAsO_4 | mm2 | 二軸 | 0.35–4.5 | OPO、電気光学。KTA 類似 | — |
| CTA | CsTiOAsO_4 | mm2 | 二軸 | 0.35–5.3 | OPO（~2 µm アイセーフ）、中赤外まで透過 | — |

## ニオブ酸・タンタル酸塩（QPM 母材）

| 略称 | 組成 | 点群 | 光学クラス | 透過域 (µm) | 主な用途 | ndispers |
|---|---|---|---|---|---|---|
| LNB_C | LiNbO_3 congruent | 3m | 負一軸 | 0.33–5.5 | コングルエント LiNbO3。導波路、電気光学変調、PPLN 母材。光損傷あり | — |
| LNB_M | LiNbO_3 5% MgO doped congruent lithium niobate | 3m | 負一軸 | 0.33–5.5 | 5%MgO:CLN。光損傷耐性が高い。PPMgLN による OPO/DFG/中赤外発生 | ✅ `crystals.MgOLN_Zelmon1997` |
| LNB_S | stoichiometric LiNbO_3 | 3m | 負一軸 | 0.33–5.5 | ストイキオメトリック LiNbO3。分極反転電圧が低く厚い PPLN が作れる | △ `crystals.SLN (1%MgO 添加体のみ。SNLO の LNB_S は無添加)` |
| LITA_C | LiTaO_3 congruent, | 3m | 正一軸 | 0.28–5.5 | コングルエント LiTaO3。PPLT、SAW、紫外側が LN より広い | — |
| LITA_S | LiTaO_3 stoichiometric, | 3m | 正一軸 | 0.28–5.5 | ストイキオメトリック LiTaO3。低電圧分極反転 | — |
| LITA_M | 0.5% MgO doped LiTaO_3 stoichiometric, | 3m | 負一軸 | 0.28–5.5 | 0.5%MgO:SLT。緑色〜UV 発生の PPSLT、光損傷耐性 | ✅ `crystals.SLT (0.5%MgO 添加ストイキオメトリック)` |
| KNBO3 | KNbO_3 | mm2 | 二軸 | 0.4–4.5 | d_eff が大きい。青色 SHG（860→430 nm）、NCPM。機械的に脆い | — |

## 中赤外カルコゲナイド

| 略称 | 組成 | 点群 | 光学クラス | 透過域 (µm) | 主な用途 | ndispers |
|---|---|---|---|---|---|---|
| AGS | AgGaS_2 Silver thiogallate | 4̄2m | 負一軸 | 0.5–13 | 中赤外 DFG/OPO（〜13 µm）。CO2 レーザー波長変換 | ✅ `crystals.AGS_Kato1996 / AGS_Takaoka1999` |
| AGSE | AgGaSe_2 Silver gallium selenide | 4̄2m | 負一軸 | 0.71–18 | 中赤外 OPO/DFG（〜18 µm）。AGS より長波長側 | ✅ `crystals.AGSe_Kato2021` |
| AGGS | AgGaGeS_4 | mm2 | 二軸 | 0.7–13 | AGS 系四元。中赤外 OPO、複屈折を調整可能 | — |
| AGGSE | AgGaGe5_Se_12 | mm2 | 二軸 | 0.6–16 | 中赤外 OPO/DFG（〜16 µm） | — |
| HGS | HgGa_2S_4 | 4̄ | 負一軸 | 0.5–13 | 中赤外 DFG。1 µm 励起が可能（AGS より熱伝導良） | — |
| GS | GaSe | 6̄m2 | 負一軸 | 0.65–18 | THz 発生・検出、中赤外 DFG。層状で劈開のみ、任意角カット不可 | — |
| LIS | LiInS_2 | mm2 | 二軸 | 0.4–12 | 1 µm 直接励起の中赤外 OPO（〜12 µm） | — |
| LISE | LiInSe_2 | mm2 | 二軸 | 0.5–12 | 1 µm 励起中赤外 OPO。LIS より長波長 | — |
| LGS | LiGaS_2 | mm2 | 二軸 | 0.32–11.6 | 1 µm 励起中赤外 DFG/OPO、THz 発生 | — |
| LGSE | LiGaSe_2 | mm2 | 二軸 | 0.37–13.2 | 中赤外 OPO/DFG | — |
| LGT | LiGaTe_2 | 4̄2m | 正一軸 | 0.5–15 | 中赤外（〜15 µm）。屈折率・非線形係数が大きい | — |
| BGS | BaGa_4S_7 Barium thiogallate | mm2 | 二軸 | 0.62–9.2 | 1 µm 励起中赤外 OPO（〜9 µm）。近年の主力 | — |
| BGSE | BaGa_4Se_7 Barium gallium selenate | m | 二軸 | 0.47–17 | 1 µm 励起中赤外 OPO（〜17 µm）。長波長中赤外の本命 | — |
| BGGSE | BaGa_2GeSe_6 | 3 | 正一軸 | 0.58–16 | 中赤外 OPO（〜16 µm） | — |
| B2GGS | Ba_2Ga_8GeS_16 | 6mm | 負一軸 | 0.42–12 | 中赤外 OPO/DFG | — |
| AAS | Ag_3AsS_3 | 3m | 負一軸 | 0.6–13 | プルースタイト。中赤外 DFG、CO2 波長変換。損傷閾値が低い | — |
| ASS | Ag_3SbS_3 | 3m | 負一軸 | 0.7–14 | ピラルギライト。中赤外 DFG | — |
| TAS | Tl_3AsSe_3 | 3m | 負一軸 | 1.25–20 | 長波長中赤外（〜20 µm）DFG、CO2 レーザー変換 | — |
| HS | HgS | 32 | 正一軸 | 0.63–13.5 | シナバー。中赤外、旋光性が非常に強い | — |
| CDSE | CdSe | 6mm | 正一軸 | 0.75–16.5 | 中赤外 OPO/DFG（〜16 µm）。CO2 励起 | — |

## ポニクタイド・半導体

| 略称 | 組成 | 点群 | 光学クラス | 透過域 (µm) | 主な用途 | ndispers |
|---|---|---|---|---|---|---|
| ZGP | ZnGeP_2 | 4̄2m | 正一軸 | 0.74–12 | 2 µm 励起の中赤外 OPO（3-5 µm）の標準材料。熱伝導率が高い | ✅ `crystals.ZGP_Das2003 / ZGP_Zelmon2001` |
| CGA | CdGeAs_2 | 4̄2m | 正一軸 | 2.4–18 | 最長波長級（〜18 µm）。d が非常に大きいが 2.4 µm 以下は吸収 | — |
| CSP | CdSiP_2 | 4̄2m | 負一軸 | 0.95–9.5 | 1 µm 励起で 6.5 µm 帯 OPO。ZGP の短波長励起代替 | — |
| GAAS | GaAs | 4̄3m | 等方 | 0.997–17 | 等方性のため OP-GaAs（配向パターン QPM）で中赤外 OPO | — |
| GAP | GaP | 4̄3m | 等方 | 0.7–12.5 | OP-GaP で QPM、THz 発生 | — |
| ZNSE | ZnSe | 4̄3m | 等方 | 0.45–18 | 等方性。中赤外窓・レンズ材料（QPM 用の配向パターン化も） | ✅ `glasses.ZnSe (線形のみ、d 係数なし)` |
| TE | Te | 32 | 正一軸 | 3.8–32 | 超長波長（〜32 µm）DFG。屈折率 ~4.8、フレネル損失大 | — |
| SC4H | SiC-4H | 6mm | 正一軸 | 0.37–5.8 | SiC 4H。QPM 候補、広帯域透過・高熱伝導 | — |
| SC6H | SiC-6H | 6mm | 正一軸 | 0.4–5.8 | SiC 6H。QPM 候補 | — |

## DIY 枠（係数はユーザーが入力）

| 略称 | 組成 | 点群 | 光学クラス | 透過域 (µm) | 主な用途 | ndispers |
|---|---|---|---|---|---|---|
| ZZ_U | — | — | — | — | ユーザー定義の一軸結晶（Sellmeier と d を自分で入力する枠） | — |
| ZZ_B | — | — | — | — | ユーザー定義の二軸結晶（同上） | — |

## ndispers にあって SNLO にないもの

SNLO は非線形結晶専用なので、線形（中心対称）材料を持たない。ndispers 側の以下は
SNLO と重ならない。

| ndispers | 組成 | 点群 | 用途 |
|---|---|---|---|
| `crystals.AlphaBBO` | α-BaB₂O₄ | 3̄m | 複屈折素子、偏光子、パルス整形 |
| `crystals.Calcite` | CaCO₃ | 3̄m | 偏光子（Glan 系）、複屈折素子 |
| `crystals.Quartz` | α-SiO₂ | 32 | 波長板、旋光子。d の絶対スケール基準 |
| `crystals.Sapphire` | Al₂O₃ | 3̄m | 窓材、Ti:Sa 母材 |
| `crystals.MgF2` | MgF₂ | 4/mmm | 真空紫外窓、波長板 |
| `crystals.YVO4` | YVO₄ | 4/mmm | 偏光子、ウォークオフ補償、レーザー母材 |
| `glasses.*` | 溶融石英, CaF₂, BaF₂, LiF, YAG, N-BK7, SF10/11/57, Si, Ge, ZnS, Diamond | 等方 | 窓材・レンズ・分散管理 |

`glasses.ZnSe` だけは SNLO 側にも `ZNSE` として存在する（SNLO は QPM 用途で
非線形係数も持つが、ndispers は屈折率のみ）。

## 未実装のうち優先度が高そうなもの

判断材料は「実際に使われている頻度 × 点群が既存クラスで賄えるか」。

| 材料 | 点群 | 既存クラス | 理由 |
|---|---|---|---|
| KTA | mm2 | `Biax_mm2` | KTP の中赤外版。3–4 µm OPO の定番 |
| RTP | mm2 | `Biax_mm2` | 高繰り返しポッケルスセル。問い合わせが多い |
| KNbO₃ | mm2 | `Biax_mm2` | d_eff が大きく青色 SHG の教科書例 |
| LNB_C（無添加 CLN） | 3m | `Uniax_3m` | PPLN の母材。MgO:LN があるのに素の CLN がない |
| LiTaO₃（congruent） | 3m | `Uniax_3m` | PPLT。LN と同形で実装コストが低い |
| BaGa₄Se₇ / BaGa₄S₇ | mm2 | `Biax_mm2` | 現在の中赤外 OPO の主力 |
| CdSiP₂ | 4̄2m | `Uniax_42m` | ZGP と同形。1 µm 励起の 6.5 µm 帯 |
| GaSe | 6̄m2 | 新規（6̄m2 が必要） | THz 用途で需要。クラス追加が要る |

いずれも「Sellmeier + d + 出典」の 1 ファイルで追加できる（`materials_plan.md` の方針）。

## 再取得の手順

このファイルは手で書かず、`tools/gen_snlo_crystal_list.py` が生成する。

```bash
uv run python tools/gen_snlo_crystal_list.py
```

Web 版 Qmix (<https://smithjj.github.io>) のページソースに埋め込まれた
`INLINE_CRYSTAL_INFO` を取ってきて、スクリプト内の `ND`（SNLO 名 → ndispers クラス）、
`USE`（用途）、`FAM`（ファミリー分け）と突き合わせる。SNLO が更新されたとき、
または ndispers に材料を追加したときに回す。SNLO 側に未知の結晶が増えると
`FAM` / `USE` の網羅チェックで assert が落ちるので、取りこぼしはしない。

抽出される各エントリは `crystal_description`, `wavelength_range`, `iso_uni_or_bi`,
`crystal_class`, `ref_ind_source`, `thermo_optic_source`, `d_eff1`, `d_eff2`,
`d_string`（d テンソル）, `d_source`, `thermal_conductivity`, `thermal_expansion`,
`specific_heat`, `density` を持つ。Sellmeier 係数そのものは含まれないので、
実装時は `ref_ind_source` の文献に当たること。

