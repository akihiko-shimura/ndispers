# Crystals

Each class carries, in its docstring, the Sellmeier equation it implements,
its validity range and the literature the coefficients come from — the same
text you get from `help(x)` or `x?`.

| Material | Classes |
|---|---|
| [β-BBO](#beta-bbo) | `BetaBBO_Eimerl1987`, `BetaBBO_Ghosh1995`, `BetaBBO_KK2010`, `BetaBBO_Tamosauskas2018` |
| [LBO](#lbo) | `LBO_Castech_xy`, `LBO_Castech_yz`, `LBO_Castech_zx`, `LBO_Ghosh1995_xy`, `LBO_Ghosh1995_yz`, `LBO_Ghosh1995_zx`, `LBO_KK1994_xy`, `LBO_KK1994_yz`, `LBO_KK1994_zx`, `LBO_KK2018_xy`, `LBO_KK2018_yz`, `LBO_KK2018_zx`, `LBO_Newlight_xy`, `LBO_Newlight_yz`, `LBO_Newlight_zx` |
| [KTP](#ktp) | `KTP_xy`, `KTP_yz`, `KTP_zx` |
| [BiBO](#bibo) | `BiBO_Miyata2009_xy`, `BiBO_Miyata2009_yz`, `BiBO_Miyata2009_zx` |
| [CLBO](#clbo) | `CLBO` |
| [KDP](#kdp) | `KDP`, `DKDP` |
| [KBBF](#kbbf) | `KBBF` |
| [RBBF](#rbbf) | `RBBF` |
| [LB4](#lb4) | `LB4` |
| [ZGP](#zgp) | `ZGP_Zelmon2001`, `ZGP_Das2003` |
| [AgGaS₂](#ags) | `AGS_Kato1996`, `AGS_Takaoka1999` |
| [AgGaSe₂](#agse) | `AGSe_Kato2021` |
| [α-BBO](#alpha-bbo) | `AlphaBBO` |
| [Calcite](#calcite) | `Calcite` |
| [Sapphire](#sapphire) | `Sapphire` |
| [α-quartz](#quartz) | `Quartz` |
| [LiIO₃](#liio3) | `LiIO3` |
| [YVO₄](#yvo4) | `YVO4` |
| [MgF₂](#mgf2) | `MgF2` |
| [SLN](#sln) | `SLN` |
| [MgO:LN](#mgoln) | `MgOLN_Zelmon1997` |
| [SLT](#slt) | `SLT` |

## β-BBO { #beta-bbo }

Barium borate. Four independent parameterisations.

::: ndispers.media.crystals.BetaBBO_Eimerl1987

::: ndispers.media.crystals.BetaBBO_Ghosh1995

::: ndispers.media.crystals.BetaBBO_KK2010

::: ndispers.media.crystals.BetaBBO_Tamosauskas2018

## LBO

Lithium triborate, biaxial: one class per principal plane, five sources.

::: ndispers.media.crystals.LBO_Castech_xy

::: ndispers.media.crystals.LBO_Castech_yz

::: ndispers.media.crystals.LBO_Castech_zx

::: ndispers.media.crystals.LBO_Ghosh1995_xy

::: ndispers.media.crystals.LBO_Ghosh1995_yz

::: ndispers.media.crystals.LBO_Ghosh1995_zx

::: ndispers.media.crystals.LBO_KK1994_xy

::: ndispers.media.crystals.LBO_KK1994_yz

::: ndispers.media.crystals.LBO_KK1994_zx

::: ndispers.media.crystals.LBO_KK2018_xy

::: ndispers.media.crystals.LBO_KK2018_yz

::: ndispers.media.crystals.LBO_KK2018_zx

::: ndispers.media.crystals.LBO_Newlight_xy

::: ndispers.media.crystals.LBO_Newlight_yz

::: ndispers.media.crystals.LBO_Newlight_zx

## KTP

Potassium titanyl phosphate, biaxial: one class per principal plane.

::: ndispers.media.crystals.KTP_xy

::: ndispers.media.crystals.KTP_yz

::: ndispers.media.crystals.KTP_zx

## BiBO

Bismuth triborate, monoclinic (point group 2), biaxial: one class per principal plane.

::: ndispers.media.crystals.BiBO_Miyata2009_xy

::: ndispers.media.crystals.BiBO_Miyata2009_yz

::: ndispers.media.crystals.BiBO_Miyata2009_zx

## CLBO

Caesium lithium borate.

::: ndispers.media.crystals.CLBO

## KDP

Potassium dihydrogen phosphate.

::: ndispers.media.crystals.KDP

::: ndispers.media.crystals.DKDP

## KBBF

Potassium beryllium fluoroborate, for deep-UV generation.

::: ndispers.media.crystals.KBBF

## RBBF

Rubidium beryllium fluoroborate.

::: ndispers.media.crystals.RBBF

## LB4

Lithium tetraborate.

::: ndispers.media.crystals.LB4

## ZGP

Zinc germanium phosphide, mid-infrared; two parameterisations.

::: ndispers.media.crystals.ZGP_Zelmon2001

::: ndispers.media.crystals.ZGP_Das2003

## AgGaS₂ { #ags }

Silver thiogallate, mid-infrared; two parameterisations.

::: ndispers.media.crystals.AGS_Kato1996

::: ndispers.media.crystals.AGS_Takaoka1999

## AgGaSe₂ { #agse }

Silver gallium selenide, mid-infrared.

::: ndispers.media.crystals.AGSe_Kato2021

## α-BBO { #alpha-bbo }

The other barium borate phase, used for birefringent optics.

::: ndispers.media.crystals.AlphaBBO

## Calcite

Calcium carbonate.

::: ndispers.media.crystals.Calcite

## Sapphire

α-Al₂O₃; centrosymmetric, birefringent optic.

::: ndispers.media.crystals.Sapphire

## α-quartz { #quartz }

Crystalline SiO₂; the absolute-scale reference of d coefficients.

::: ndispers.media.crystals.Quartz

## LiIO₃ { #liio3 }

Lithium iodate, point group 6.

::: ndispers.media.crystals.LiIO3

## YVO₄ { #yvo4 }

Yttrium orthovanadate; centrosymmetric, strongly positive uniaxial - polarizers and displacers.

::: ndispers.media.crystals.YVO4

## MgF₂ { #mgf2 }

Magnesium fluoride; centrosymmetric, positive uniaxial - VUV windows and wave plates.

::: ndispers.media.crystals.MgF2

## SLN

1% MgO-doped stoichiometric lithium niobate (extraordinary ray only).

::: ndispers.media.crystals.SLN

## MgO:LN { #mgoln }

5 mol% MgO-doped congruent LiNbO₃, both rays, at 21 °C.

::: ndispers.media.crystals.MgOLN_Zelmon1997

## SLT

0.5% MgO-doped stoichiometric lithium tantalate.

::: ndispers.media.crystals.SLT
