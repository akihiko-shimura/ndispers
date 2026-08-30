"""
Regenerate docs/dev/snlo_crystals.md, the SNLO crystal list with ndispers
implementation status.

SNLO ships its crystal database inside the app, but the web version of Qmix
embeds the same table as a JSON literal (`INLINE_CRYSTAL_INFO`) in its page
source. This scrapes that, joins it with the ND/USE/FAM tables below, and
writes the doc.

Run after an SNLO release, or after adding a medium to ndispers:

    uv run python tools/gen_snlo_crystal_list.py

Composition, point group, transparency range and optical class come from SNLO
verbatim. ND maps SNLO names to ndispers classes and USE holds the
application blurbs - both are hand-maintained here, not from SNLO.
"""
import json
import urllib.request
from datetime import date
from pathlib import Path

SRC = 'https://smithjj.github.io'
ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'docs' / 'dev' / 'snlo_crystals.md'


def fetch():
    """Pull INLINE_CRYSTAL_INFO out of the web-Qmix page by brace matching."""
    with urllib.request.urlopen(SRC) as r:
        page = r.read().decode('utf8', 'replace')
    key = 'const INLINE_CRYSTAL_INFO = '
    i = page.index(key) + len(key)
    depth = 0
    for j in range(i, len(page)):
        if page[j] == '{':
            depth += 1
        elif page[j] == '}':
            depth -= 1
            if depth == 0:
                break
    return json.loads(page[i:j + 1])


J = fetch()

# SNLO name -> ndispers class path ('' = not implemented, '~' prefix = partial)
ND = {
 'AGS':'crystals.AGS_Kato1996 / AGS_Takaoka1999',
 'AGSE':'crystals.AGSe_Kato2021',
 'BBO':'crystals.BetaBBO_{Eimerl1987,Ghosh1995,KK2010,Tamosauskas2018}',
 'BIBO':'crystals.BiBO_Miyata2009_{xy,yz,zx}',
 'CLBO':'crystals.CLBO',
 'DKDP':'crystals.DKDP',
 'KBBF':'crystals.KBBF',
 'KDP':'crystals.KDP',
 'KTP_F':'crystals.KTP_{xy,yz,zx}',
 'KTP_H':'crystals.KTP_{xy,yz,zx}',
 'KTP_K':'crystals.KTP_{xy,yz,zx}',
 'LB4':'crystals.LB4',
 'LBO':'crystals.LBO_{Castech,Ghosh1995,KK1994,KK2018,Newlight}_{xy,yz,zx}',
 'LIIO3':'crystals.LiIO3',
 'LNB_M':'crystals.MgOLN_Zelmon1997',
 'LNB_S':'~crystals.SLN (1%MgO 添加体のみ。SNLO の LNB_S は無添加)',
 'LITA_M':'crystals.SLT (0.5%MgO 添加ストイキオメトリック)',
 'RBBF':'crystals.RBBF',
 'ZGP':'crystals.ZGP_Das2003 / ZGP_Zelmon2001',
 'ZNSE':'glasses.ZnSe (線形のみ、d 係数なし)',
}

USE = {
 # KDP family
 'ADP':'水溶液成長の古典材料。可視 SHG、ポッケルスセル',
 'ADA':'ADP の As 置換体。1064 nm 付近で 90° 位相整合 SHG',
 'DADP':'重水素化 ADP。IR 吸収低減、電気光学変調',
 'DADA':'重水素化 ADA。低温 NCPM SHG',
 'KDP':'大口径成長が容易。慣性核融合レーザーの 2ω/3ω 変換、ポッケルスセル',
 'KDA':'KDP の As 置換体。1064 nm SHG、電気光学',
 'DKDP':'KD*P。1064 nm の 3ω/4ω、大型ポッケルスセル。KDP より IR 吸収が小さい',
 'RDP':'電気光学変調器（低半波電圧）',
 'RDA':'1064 nm で 90° 位相整合 SHG、Q スイッチ',
 'DRDP':'重水素化 RDP。電気光学変調',
 'DRDA':'重水素化 RDA。NCPM SHG',
 'CDA':'1064 nm 90° NCPM SHG の定番（温度整合）',
 'DCDA':'重水素化 CDA。1064 nm NCPM SHG',
 # borates
 'BBO':'β-BBO。UV SHG/THG/FHG、fs OPA、ポッケルスセル。最も汎用',
 'LBO':'1064 nm SHG/THG、高平均出力の緑色発生、NCPM OPO。損傷閾値が高い',
 'CLBO':'266/213 nm の深紫外発生（4ω/5ω）。潮解性あり',
 'CBO':'1064 nm THG、深紫外。LBO 類似で複屈折が大きい',
 'LB4':'深紫外 SHG（~240 nm）、高損傷閾値、非潮解',
 'KBBF':'177.3 nm の直接 SHG（唯一の実用材料）。層状剥離のためプリズム結合が必要',
 'RBBF':'KBBF 類似の深紫外材料。成長性が改善',
 'CBBF':'KBBF 系。深紫外 SHG',
 'KABO':'深紫外 SHG（~200 nm）、非潮解',
 'KBO':'深紫外 SHG（~217 nm）。潮解性が強い',
 'BIBO':'1064 nm SHG で d_eff 大。fs OPA/OPO、緑色発生',
 'LCB':'紫外 SHG、非潮解',
 'NLBO':'紫外 SHG',
 'BBPO':'紫外 SHG',
 'BPO':'真空紫外まで透過。深紫外 SHG',
 'LBGO':'紫外 SHG、周期分極反転（QPM）候補',
 'GCOB':'自己周波数変換レーザー母材、1064 nm SHG',
 'YCOB':'大口径成長可。高平均出力 SHG、自己周波数変換',
 'TCOB':'COB 系。SHG、自己周波数変換',
 'LIIO3':'LiIO3。可視〜近赤外 SHG/OPO の古典材料。潮解性・低損傷閾値',
 'DLAP':'有機系。可視 SHG、高損傷閾値',
 'LFM':'有機系。紫外〜可視 SHG',
 # KTP family
 'KTP_F':'フラックス成長 KTP。1064 nm SHG、OPO、電気光学。PPKTP 母材',
 'KTP_H':'水熱合成 KTP。灰色化（gray tracking）耐性が高い',
 'KTP_K':'フラックス成長 KTP（Kato の Sellmeier）',
 'KTA_1':'KTP より中赤外透過が広い。3-4 µm OPO',
 'KTA_2':'KTA（別 Sellmeier）',
 'KTA_3':'KTA（別 Sellmeier）',
 'RTP':'高繰り返し電気光学 Q スイッチ／ポッケルスセル。低イオン伝導',
 'RTA':'OPO、電気光学。KTA 類似',
 'CTA':'OPO（~2 µm アイセーフ）、中赤外まで透過',
 # niobates / tantalates
 'LNB_C':'コングルエント LiNbO3。導波路、電気光学変調、PPLN 母材。光損傷あり',
 'LNB_M':'5%MgO:CLN。光損傷耐性が高い。PPMgLN による OPO/DFG/中赤外発生',
 'LNB_S':'ストイキオメトリック LiNbO3。分極反転電圧が低く厚い PPLN が作れる',
 'LITA_C':'コングルエント LiTaO3。PPLT、SAW、紫外側が LN より広い',
 'LITA_S':'ストイキオメトリック LiTaO3。低電圧分極反転',
 'LITA_M':'0.5%MgO:SLT。緑色〜UV 発生の PPSLT、光損傷耐性',
 'KNBO3':'d_eff が大きい。青色 SHG（860→430 nm）、NCPM。機械的に脆い',
 # mid-IR chalcogenides
 'AGS':'中赤外 DFG/OPO（〜13 µm）。CO2 レーザー波長変換',
 'AGSE':'中赤外 OPO/DFG（〜18 µm）。AGS より長波長側',
 'AGGS':'AGS 系四元。中赤外 OPO、複屈折を調整可能',
 'AGGSE':'中赤外 OPO/DFG（〜16 µm）',
 'HGS':'中赤外 DFG。1 µm 励起が可能（AGS より熱伝導良）',
 'GS':'THz 発生・検出、中赤外 DFG。層状で劈開のみ、任意角カット不可',
 'LIS':'1 µm 直接励起の中赤外 OPO（〜12 µm）',
 'LISE':'1 µm 励起中赤外 OPO。LIS より長波長',
 'LGS':'1 µm 励起中赤外 DFG/OPO、THz 発生',
 'LGSE':'中赤外 OPO/DFG',
 'LGT':'中赤外（〜15 µm）。屈折率・非線形係数が大きい',
 'BGS':'1 µm 励起中赤外 OPO（〜9 µm）。近年の主力',
 'BGSE':'1 µm 励起中赤外 OPO（〜17 µm）。長波長中赤外の本命',
 'BGGSE':'中赤外 OPO（〜16 µm）',
 'B2GGS':'中赤外 OPO/DFG',
 'AAS':'プルースタイト。中赤外 DFG、CO2 波長変換。損傷閾値が低い',
 'ASS':'ピラルギライト。中赤外 DFG',
 'TAS':'長波長中赤外（〜20 µm）DFG、CO2 レーザー変換',
 'HS':'シナバー。中赤外、旋光性が非常に強い',
 'CDSE':'中赤外 OPO/DFG（〜16 µm）。CO2 励起',
 # pnictides / semiconductors
 'ZGP':'2 µm 励起の中赤外 OPO（3-5 µm）の標準材料。熱伝導率が高い',
 'CGA':'最長波長級（〜18 µm）。d が非常に大きいが 2.4 µm 以下は吸収',
 'CSP':'1 µm 励起で 6.5 µm 帯 OPO。ZGP の短波長励起代替',
 'GAAS':'等方性のため OP-GaAs（配向パターン QPM）で中赤外 OPO',
 'GAP':'OP-GaP で QPM、THz 発生',
 'ZNSE':'等方性。中赤外窓・レンズ材料（QPM 用の配向パターン化も）',
 'TE':'超長波長（〜32 µm）DFG。屈折率 ~4.8、フレネル損失大',
 'SC4H':'SiC 4H。QPM 候補、広帯域透過・高熱伝導',
 'SC6H':'SiC 6H。QPM 候補',
 # other oxides / misc
 'GEO':'α-GeO2。石英類似で複屈折が大きい。紫外〜中赤外の複屈折素子',
 'LGN':'ランガサイト系。中赤外 OPO/DFG、圧電',
 'LGTO':'ランガテイト。中赤外、圧電',
 'CTW':'テルル酸タングステン。可視〜中赤外 SHG',
 'NTW':'テルル酸タングステン。可視〜中赤外 SHG',
 # DIY
 'ZZ_U':'ユーザー定義の一軸結晶（Sellmeier と d を自分で入力する枠）',
 'ZZ_B':'ユーザー定義の二軸結晶（同上）',
}

FAM = [
 ('リン酸・ヒ酸塩系（KDP ファミリー, 4̄2m）',
  ['ADP','ADA','DADP','DADA','KDP','KDA','DKDP','RDP','RDA','DRDP','DRDA','CDA','DCDA']),
 ('ボレート系（UV／深紫外の主力）',
  ['BBO','LBO','CLBO','CBO','LB4','KBBF','RBBF','CBBF','KABO','KBO','BIBO','LCB','NLBO',
   'BBPO','BPO','LBGO','GCOB','YCOB','TCOB']),
 ('ヨウ素酸塩・有機・その他の酸化物',
  ['LIIO3','DLAP','LFM','GEO','LGN','LGTO','CTW','NTW']),
 ('KTP ファミリー（mm2 二軸）',
  ['KTP_F','KTP_H','KTP_K','KTA_1','KTA_2','KTA_3','RTP','RTA','CTA']),
 ('ニオブ酸・タンタル酸塩（QPM 母材）',
  ['LNB_C','LNB_M','LNB_S','LITA_C','LITA_S','LITA_M','KNBO3']),
 ('中赤外カルコゲナイド',
  ['AGS','AGSE','AGGS','AGGSE','HGS','GS','LIS','LISE','LGS','LGSE','LGT',
   'BGS','BGSE','BGGSE','B2GGS','AAS','ASS','TAS','HS','CDSE']),
 ('ポニクタイド・半導体',
  ['ZGP','CGA','CSP','GAAS','GAP','ZNSE','TE','SC4H','SC6H']),
 ('DIY 枠（係数はユーザーが入力）',
  ['ZZ_U','ZZ_B']),
]

covered = [k for _, ks in FAM for k in ks]
assert len(covered) == len(set(covered)), 'dup in FAM'
assert set(covered) == set(J), f'missing={set(J)-set(covered)} extra={set(covered)-set(J)}'
assert set(USE) == set(J), f'USE mismatch: {set(J)^set(USE)}'

OPT = {'neg. uniaxial':'負一軸','pos. uniaxial':'正一軸','positive uniaxial':'正一軸',
       'biaxial':'二軸','isotropic':'等方','DIY':'—'}

# SNLO の ASCII 表記 -> 国際記号。'6m2' は 6mm/622 と別物で 6-bar-m-2 のこと。
PG = {'': '—', '4(bar)': '4̄', '4(bar)2m': '4̄2m', '4(bar)3m': '4̄3m', '6m2': '6̄m2'}

def row(k):
    v = J[k]
    desc = v.get('crystal_description','').split(' or ')[0].strip() or '—'
    pg = PG.get(v.get('crystal_class',''), v.get('crystal_class') or '—')
    opt = OPT.get(v.get('iso_uni_or_bi',''), v.get('iso_uni_or_bi','—'))
    lo, hi = v.get('wavelength_range') or (0,0)
    rng = f'{lo/1000:g}–{hi/1000:g}' if lo else '—'
    nd = ND.get(k, '')
    if nd.startswith('~'):
        mark = f'△ `{nd[1:]}`'
    elif nd:
        mark = f'✅ `{nd}`'
    else:
        mark = '—'
    return f'| {k} | {desc} | {pg} | {opt} | {rng} | {USE[k]} | {mark} |'

out = []
HDR = '''# SNLO 結晶リスト（ndispers 実装状況つき）

SNLO v80（2025-01-22）が内蔵する全結晶。出典は Web 版 Qmix
(<https://smithjj.github.io>) に埋め込まれた `INLINE_CRYSTAL_INFO` を
そのまま抽出したもの（{DATE} 取得）。組成・点群・透過域・光学クラスは
SNLO のデータそのままで、加筆していない（点群だけ SNLO の ASCII 表記
`4(bar)2m` / `6m2` を国際記号 4̄2m / 6̄m2 に直した）。

- 実体は **85 結晶**（+ DIY 枠 `ZZ_U` / `ZZ_B` の 2 つ）。
- `ndispers` 列: ✅ 実装済み / △ 部分的（別組成で代用）/ — 未実装。
  現状 **{NOK} エントリが ✅、{NPART} つが △、{NMISS} が未実装**（DIY 枠 2 つを除く）。
- 「主な用途」は SNLO のデータには含まれていない。一般的な文献知識から
  書いた要約で、**出典照合はしていない**。材料選定の当たりをつける用途に留める。
- 透過域は SNLO の Sellmeier 有効範囲（µm）であって、結晶の実透過域とは限らない。

## 対応表の読み方

SNLO は同一組成で Sellmeier が複数あるとき別エントリにする（KTP_F/H/K、KTA_1/2/3、
LNB_C/M/S など）。ndispers も同じ方針（`AGS_Kato1996` と `AGS_Takaoka1999` など）なので、
両者は素直に 1 対 1 で対応づけられる。ただし文献の選び方は一致していない。
'''

NOK = sum(1 for k, v in ND.items() if not v.startswith('~'))
NPART = sum(1 for v in ND.values() if v.startswith('~'))
NMISS = len(J) - len(ND) - 2
out.append(HDR.format(NOK=NOK, NPART=NPART, NMISS=NMISS, DATE=date.today()))

for title, keys in FAM:
    out.append(f'\n## {title}\n')
    out.append('| 略称 | 組成 | 点群 | 光学クラス | 透過域 (µm) | 主な用途 | ndispers |')
    out.append('|---|---|---|---|---|---|---|')
    out += [row(k) for k in keys]

nd_only = '''
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
'''
out.append(nd_only)
OUT.write_text('\n'.join(out) + '\n')
print(f'wrote {OUT.relative_to(ROOT)} ({len(J)} crystals, {NOK} implemented)')
