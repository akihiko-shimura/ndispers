# Validation

Every dispersion formula in this package is a transcription of somebody else's
measurement. This page records what has been checked against those sources and
what has not, so that a number taken from `ndispers` can be traced back to the
paper it came from.

Two kinds of check are reported here. A **fidelity check** compares the package
against the same paper whose coefficients it implements: it tests the
transcription, not the physics. An **independent check** compares it against a
different source — another group's measurement, a handbook table, or a second
parameterisation of the same crystal — and so tests both.

Everything below was recomputed at the version stated at the foot of this page.
Where a value is a phase-matching angle it comes from `pmAngles_sfg`, where it
is an index from `n`, and so on; nothing is quoted from memory. Where no
citable source could be found the row says so rather than guessing — an
omission is more useful than a fabricated reference.

## Refractive indices

| Medium | Quantity | ndispers | Literature | Source | |
|---|---|---|---|---|---|
| Sapphire | n_o, n_e at 589 nm | 1.76808, 1.76000 | 1.7681, 1.7600 | Malitson & Dodge 1972 | ✓ |
| α-Quartz | n_o, n_e at 589 nm | 1.54424, 1.55387 | 1.5443, 1.5534 | handbook value | ✓ |
| 5% MgO:LiNbO₃ | n_o, n_e at 1064 nm | 2.22884, 2.14737 | negative uniaxial, n_o > n_e | Zelmon *et al.* 1997, Fig. 2 | ✓ (see caveat 1) |
| MgF₂ | n_e − n_o at 1064 nm | +0.01159 | ≈ +0.0116 | Dodge 1984 | ✓ |
| YVO₄ | n_e − n_o at 1064 nm | +0.2082 | 0.204–0.21 | vendor and handbook data | ✓ |
| BiBO | n_x, n_y, n_z at 632.8 nm | 1.77668, 1.80641, 1.94582 | 1.77668, 1.80641, 1.94582 | Petrov *et al.* 2010, Table 3 | ✓ exact |
| BiBO | n_x, n_z at 1064 nm | 1.75752, 1.91711 | 1.75752, 1.91711 | Petrov *et al.* 2010, Table 3 | ✓ exact |
| ZnSe | n at 10.6 µm | 2.40266 | 2.403 | Connolly *et al.* 1979 | ✓ |
| ZnS | n at 1.064 µm | 2.28828 | 2.289 | Debenham 1984 | ✓ |
| Si | n at 5 µm | 3.42207 | 3.4255 | handbook value | ✓ 0.1 % |
| Ge | n at 10 µm | 4.00431 | 4.0043 | Barnes & Piltch 1979 | ✓ |
| Diamond | n at 587.6 nm | 2.41748 | 2.4175 | Peter 1923 | ✓ |
| YAG | n at 1064 nm | 1.81465 | 1.8145 | Zelmon *et al.* 1998 | ✓ |
| N-BK7, SF10, SF11, SF57 | n_d at 587.6 nm | 1.51680, 1.72825, 1.78471, 1.84666 | 1.5168, 1.72825, 1.78472, 1.84666 | SCHOTT catalogue headers | ✓ |
| LiF | n at 589.3 nm | 1.39211 | 1.3921 | Li 1976 | ✓ |
| BaF₂ | n at 589.3 nm | 1.47443 | 1.4744 | Malitson 1964 | ✓ |

## Phase-matching angles

| Crystal | Interaction | ndispers | Literature | Source | |
|---|---|---|---|---|---|
| β-BBO (Eimerl) | ooe, 1064 nm SHG | 22.884° | ≈ 22.8° | standard value | ✓ |
| β-BBO (Ghosh, KK2010) | ooe, 1064 nm SHG | 22.771°, 22.781° | ≈ 22.8° | standard value | ✓ |
| β-BBO (Tamosauskas) | ooe, 800 nm SHG | 29.216° | ≈ 29.2° | standard value | ✓ |
| LBO (KK2018) | ooe, 1064 nm SHG, xy | φ = 11.30° | φ ≈ 11.4° | Opt. Express 34, 6337 (2026) | ✓ |
| KTP | eoe, 1064 nm SHG, xy | φ = 23.70° | ≈ 23.5° | standard value | ✓ |
| KDP | ooe / oee, 1064 nm SHG | 41.07° / 58.91° | ≈ 41° / 59° | standard value | ✓ |
| DKDP | ooe / oee, 1064 nm SHG | 35.895° / 53.187° | 35.85° / 53.10° (model) | Ghosh 1992, Table III | ✓ fidelity, < 0.1° |
| CLBO | ooe / eoe, 1.0642 µm SHG | 29.240° / 42.395° | 29.2° / 42.4° calc, 29.5° / 42.4° meas | Umemura & Kato, ASSL 1999, Table 1 | ✓ ≤ 0.3° |
| CLBO | ooe / eoe, 355 nm SFG | 38.730° / 48.506° | 38.9° / 48.7° calc | same | ✓ ≤ 0.4° |
| CLBO | ooe, 266 nm SHG | 61.410° | 61.4° calc, 61.8° meas | same | ✓ |
| LiIO₃ | ooe, 1064 nm SHG | 30.205° | 30.2° | Eckardt *et al.* 1990, Table II | ✓ exact |
| Li₂B₄O₇ | ooe, 697.3 nm SHG | 43.43° | 45° ± 5° | Komatsu *et al.* 1997, Table I | ✓ |
| Li₂B₄O₇ | ooe, 532 nm SHG | 64.97° | 65.0° calc, 66° ± 5° meas | Sugawara 1998; Komatsu 1997 | ✓ |
| Li₂B₄O₇ | ooe, 266 + 1064 nm | 73.89° | 73° ± 5° | Komatsu *et al.* 1997, Table I | ✓ |
| AgGaS₂ (Kato 1996) | ooe, 10.591 µm SHG | 68.503° | 68.5° | Kato *et al.* 2019, Table 1 | ✓ exact |
| AgGaS₂ (Kato 1996) | ooe / eoe, 5.2955 µm SHG | 32.659° / 50.388° | 32.6° / 50.1° | same | ✓ ≤ 0.3° |
| AgGaSe₂ | ooe, 10.6 µm SHG | 55.558° | 55.5° calculated | Ionin *et al.* arXiv:1710.09601 | ✓ (see caveat 5) |
| ZnGeP₂ (Das 2003) | eeo, 10.6 µm SHG | 81.395° | — | see caveat 4 | ⚠ |
| BiBO (yz plane) | eeo, 1064 nm SHG | 168.811°, d_eff 3.62 pm/V | ≈ 168.5°, ≈ 3.5 pm/V | Petrov *et al.* 2010 | ✓ (see caveat 3) |
| LBO (yz, zx planes) | 1064 nm SHG | 20.35°, 32.26° | — | **not found in the available literature** | — |

## Noncritical phase matching

These sit exactly on the 0–90° endpoint and were reported as *no solution* before v0.13.0.

| Crystal | Quantity | ndispers | Literature | Source | |
|---|---|---|---|---|---|
| LBO (KK2018) | 90° temperature, 1064 nm SHG, xy | 147.0 °C | 149 °C calc, 148–151 °C observed | Kato, Grechin & Umemura 2018, Table 1 | ✓ |
| Li₂B₄O₇ | shortest type-I SHG | 487.63 → 243.82 nm | 487.6 → 243.8 nm at θ = 90° | Komatsu *et al.* 1997, Table I | ✓ exact |
| CLBO | shortest type-I SHG | 473.43 → 236.72 nm | 236.8 nm | Umemura & Kato 1997 | ✓ 0.08 nm |
| BiBO (yz) | type-I locus at θ = 90° | 545.6 nm | 545.7 nm | Petrov *et al.* 2010, p. 61 | ✓ |
| ZnGeP₂ (Das 2003) | type-I long-wave limit | 10.74 µm | 10.78 µm | Kato 1997 | ✓ |
| KTP (x-cut: θ = 90°, φ = 0) | type-II OPO signal / idler, 1064 nm pump (`tuning_dfg`) | 1570.8 / 3297.6 nm | ≈ 1.57 / 3.3 µm | standard value for the noncritical KTP OPO | ✓ |

## Derived quantities

| Medium | Quantity | ndispers | Literature | Source | |
|---|---|---|---|---|---|
| β-BBO | walk-off of the 532 nm e-ray at 22.88° | 3.192° | ≈ 3.2° | standard value | ✓ |
| β-BBO | d_eff, ooe, φ = 30° | 1.958 pm/V | ≈ 2.0 pm/V | Shoji *et al.* 1999 via Miller scaling | ✓ |
| LiIO₃ | walk-off at the SHG angle | 4.260° | 4.26° | Eckardt *et al.* 1990, Table II | ✓ exact |
| Li₂B₄O₇ | walk-off, 532 nm SHG | 1.660° | 1.66° | Sugawara 1998, §3.2 | ✓ exact |
| KTP | d_eff, type II, xy | 3.456 pm/V | 3.2–3.5 pm/V | Shoji *et al.* 1997 via Miller scaling | ✓ |
| AgGaS₂ | d₃₆ Miller-scaled 1.064 → 10.6 µm | 11.09 pm/V | 11.2 pm/V measured at 10.6 µm | Roberts 1992, Table V | ✓ — an independent test of Miller's rule over a tenfold span |
| 5% MgO:LiNbO₃ | first-order QPM period, 1064 nm SHG | 6.995 µm | 6.9–7.0 µm | standard PPLN value | ✓ |
| Fused silica | GVD at 800 nm / 1550 nm | +36.16 / −27.95 fs²/mm | ≈ +36 / ≈ −27 | standard values | ✓ |
| SF10 | GVD at 800 nm | +159.46 fs²/mm | ≈ +160 | standard value | ✓ |
| Sapphire | GVD at 800 nm, o-ray | +58.04 fs²/mm | ≈ +58 | standard value | ✓ |
| SF11, SF57 | GVD at 800 nm | +189.63, +223.58 fs²/mm | **no published value found** | — | derived only |

---

# Module status

Several media come in more than one parameterisation. They are not
interchangeable: each is faithful to its own source, and some sources are
better than others. This table says which to reach for.

| Status | Meaning |
|---|---|
| **recommended** | the default choice for that material |
| **faithful** | correct transcription of its source, but the source has a stated limitation; keep it to reproduce that source, not to get the best number |
| **deprecated** | superseded with no remaining advantage; slated for removal |

| Module | Status | Why |
|---|---|---|
| `LBO_KK2018_*` | **recommended** | the only LBO set that reproduces the observed 90° temperature (147.0 °C against 148–151 °C observed) |
| `LBO_Ghosh1995_*`, `LBO_Newlight_*` | faithful | 130 °C and 158 °C for the same point — usable for temperature work only as a rough guide |
| `LBO_Castech_*` | faithful | its dn/dT contradicts the 148 °C printed on the same data sheet, giving 240 °C. Room-temperature indices agree with every other set to 7 × 10⁻⁵; use it to reproduce the vendor sheet, not to tune temperature |
| `LBO_KK1994_*` | faithful | reaches no 90° solution at any temperature. Kato, Grechin & Umemura published the 2018 formula to correct exactly this |
| `AGS_Kato1996` | **recommended** | reproduces Kato *et al.* 2019's table to ≤ 0.3° |
| `AGS_Takaoka1999` | **deprecated** | deviates up to 1.35° on the same table, and the thermo-optic formula that is the paper's actual contribution is not implemented here, so nothing is gained over `AGS_Kato1996` |
| `ZGP_Das2003` | **recommended** | its type-I limit, 10.74 µm, matches Kato 1997's 10.78 µm, so it covers the CO₂ line |
| `ZGP_Zelmon2001` | faithful | valid over its stated 2–9 µm; do not use it near 10.6 µm, where its limit falls 0.2 µm short |
| `MgOLN_Zelmon1997` | **recommended** for birefringent phase matching | both rays, so ooe/oee d_eff can be evaluated |
| `SLN` | **recommended** for quasi-phase matching | *not* superseded by `MgOLN_Zelmon1997`: it is the only LiNbO₃ here with a temperature-dependent index, and the QPM period it gives moves from 7.03 µm at 25 °C to 6.74 µm at 200 °C, where `MgOLN_Zelmon1997` returns 6.995 µm at every temperature. Its *absolute* dn/dT is not physical (see caveat 10); the index *differences* that QPM depends on are what Gayer's coefficients were fitted to |
| `BetaBBO_*` (four) | all usable | they agree to 0.1° on the 1064 nm SHG angle; `Tamosauskas2018` is closest at 800 nm |

`LBO_EKSMA` was removed in 0.14.0. It had never been exported from
`ndispers.media.crystals`, so no code could reach it.

---

# Caveats

These came out of the same audit. None is a defect in the package's arithmetic;
each is a property of a source, a convention, or a genuine physical subtlety
that is easy to trip over. Each is also stated in the docstring of the class it
affects, so it reaches anyone who runs `help()`.

## 1. Zelmon 1997's Table 2 has its column headings interchanged

The printed table of Sellmeier coefficients for 5 mol% MgO:LiNbO₃ heads its
columns n_e, n_o in an order that, taken literally, makes LiNbO₃ a *positive*
uniaxial crystal. It is negative, the paper's own Fig. 2 shows n_o above n_e,
and the undoped Table 1 has the larger leading coefficient in the n_o column.
`MgOLN_Zelmon1997` assigns the coefficients by the physics, giving n_o = 2.229
and n_e = 2.147 at 1064 nm. **If you have typed that table in by hand, check
which way round you have it.** A test now requires every uniaxial crystal's
Sellmeier sets to agree with the optic sign stated in its docstring.

## 2. LBO: two of the five parameterisations cannot do temperature tuning

`LBO_Castech`'s thermo-optic coefficients contradict the noncritical
temperature printed in Table 3 of the same data sheet: they put 90° phase
matching for 1064 nm SHG at about 240 °C against the sheet's own 148 °C,
because the wavelength-independent dn_y/dT makes d(Δn)/dT roughly 1.8 times too
small. `LBO_KK1994` never reaches 90° at any temperature — which is precisely
why Kato, Grechin and Umemura published the 2018 formula. **Use `LBO_KK2018`
for temperature work** (147.0 °C). The other classes remain faithful to their
sources, which is the point of having them.

## 3. BiBO: two plane-naming conventions are in circulation

This package uses the dielectric frame with the two-fold axis along x, in which
its `zx` plane is the x–z plane of Petrov *et al.* 2010. The cut usually quoted
for 1064 nm SHG — "xz, θ = 168.5°, d_eff ≈ 3.5 pm/V" — is written in the other
convention (two-fold axis along y) and is this package's **`yz`** plane at
θ = 168.8°, where d_eff comes out 3.62 pm/V.

Separately, because point group 2 has no mirror plane through z, **θ and
180° − θ phase-match identically but give different d_eff**: for 800 nm type-I
SHG in the yz plane, 1.1 pm/V at 28.8° against 3.9 pm/V at 151.2°.
`pmAngles_sfg` returns the 0–90° root, so evaluate `deff_sfg` at π − θ as well.

## 4. ZnGeP₂: the two Sellmeier sets straddle the 10.6 µm boundary

`ZGP_Das2003`'s type-I limit is 10.74 µm and it phase-matches the CO₂ line at
81.4°; `ZGP_Zelmon2001`'s limit is 10.55 µm and it finds no solution at all,
falling short by 5 × 10⁻⁴ in index. Published ZGP sets differ by 10–20° near
the range boundaries, where dθ/dn is enormous. Treat any ZGP angle within about
0.2 µm of the limit as indicative only.

## 5. Where the literature itself disagrees

AgGaSe₂'s 10.6 µm SHG angle is 55.5° by one calculation and 57.5° in a
handbook table — a two-degree spread between sources, not a package error.
Similarly, β-BBO's d₂₂ has been quoted between 1.6 and 2.3 pm/V depending on
the calibration, and the absolute scale of second-order coefficients shifted
when Shoji *et al.* 1997 corrected for multiple reflection in plane-parallel
samples.

## 6. Miller's rule is an estimate, and its standing varies by crystal

`d_sfg` and `deff_sfg` scale the tensor components from one reference
measurement using Miller's rule. That rule is supported for BBO, KDP and
LiNbO₃ (Alford & Smith 2001), rough for KTP — Shoji *et al.* 1997 found
Miller's Δ to change by 35 % between 1.06 and 1.31 µm — and untested for LBO
and most of the rest. Each crystal's `_d_note` says which case it is. The
AgGaS₂ row in the table above is the one place where the rule could be tested
here directly, and it held to 1 %.

## 7. d_eff can vanish for an interaction that phase-matches

KTP `zx` `eeo` and LBO `zx` `ooe` have real phase-matching solutions with
d_eff identically zero: the tensor has no component coupling those
polarizations in that plane. This is correct physics, not a failure, and the
phase-matching app now says so rather than printing a bare "0 pm/V".

## 8. `qpm_period_sfg` at θ = 0 with an extraordinary wave

On the optic axis the extraordinary index degenerates to the ordinary one, so
asking for a QPM period at `angle_rad = 0` with `'e'` silently computes it from
n_o: 5.91 µm instead of 7.00 µm for 1064 nm SHG in PPLN. Use `angle_rad = π/2`,
which is the physical geometry for propagation perpendicular to the polar axis.

## 9. Transcription traps found in the sources themselves

The Handbook of Optics prints MgF₂'s extraordinary infrared pole as 12.771995;
Dodge 1984 gives 23.771995, which is what reproduces the measured n_e and
birefringence, and is what this package uses. LiIO₃'s ordinary UV pole is
0.0350823 in the Handbook and 0.0350832 in refractiveindex.info's transcription
of the same set (a difference of 3 × 10⁻⁷ in n).

## 10. DFG / OPA / OPO (0.17) is validated by construction, not yet by tuning curves

Every `_dfg` method is the SFG one evaluated at (signal, idler), and the test
suite checks that equality exactly on eleven crystals; `tuning_dfg` is checked
by round trip against `pmAngles_dfg` (angle → wavelength → the same angle).
The only literature tuning point in the tables is the noncritical KTP OPO
above. Published angle- and temperature-tuning curves (BBO at 355 nm, LBO,
PPLN after Myers *et al.* 1995) have not been compared because the sources
were not at hand; when they are, the rows go in the tables above.

## 11. Coverage limits

- **SLN** holds only d₃₃ and has no ordinary-ray Sellmeier equation, so only
  eee (quasi-phase-matched) interactions can be evaluated; every combination
  involving an o-wave correctly reports no solution.
- Media whose Sellmeier equation carries **no temperature term** accept the
  temperature argument and ignore it; their `dndT` returns 0. The docstring of
  each says so, and quotes a literature dn/dT where one is available.
- Biaxial crystals are limited to the three **principal dielectric planes**.
- The KDP-family **fifth-harmonic** noncritical temperature from Ghosh 1992's
  dispersion sits 15 °C below the measured value (Cui *et al.* 2022), although
  the fourth-harmonic cut-off wavelengths for KDP and DKDP match that paper
  exactly. The deviation is confined to the 215 nm end.

---

*Recomputed against ndispers 0.17.0. To reproduce any row, call the method
named in it; the package's own test suite pins the values marked "exact".*
