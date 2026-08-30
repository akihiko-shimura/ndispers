# Conventions and definitions

The facts on this page hold across the whole package. Everything here is
stated as implemented in the code; per-crystal specifics (coefficients,
validity ranges, references) live in each class docstring and in the
[media catalog](api/crystals.md). `ndispers.catalog()` returns the same
census programmatically — name, kind, point group, plane, validity range and
nonlinear-tensor components per medium.

A wavelength more than 1% outside a medium's Sellmeier validity range raises
a `ndispers.ValidityWarning` (the value is still returned — it is an
extrapolation, not an error). Silence it with
`warnings.simplefilter('ignore', ndispers.ValidityWarning)` when the
extrapolation is deliberate; `tuning_dfg` already does so internally.

## Units

| Quantity | Symbol / method | Unit |
|---|---|---|
| Wavelength (in vacuum) | `wl_um` | µm |
| Angle of the wavevector | `theta_rad`, `phi_rad` | rad |
| Temperature | `T_degC` | °C |
| Refractive index | `n` | — |
| Group delay | `GD` | fs/mm |
| Group velocity | `GV` | µm/fs |
| Group index | `ng` | — |
| Group-velocity dispersion | `GVD` | fs²/mm |
| Third-order dispersion | `TOD` | fs³/mm |
| Thermo-optic coefficient | `dndT` | 1/K |
| Walk-off angle | `woa_theta`, `woa_phi` | rad |
| Wavevector mismatch | `dk_sfg` | rad/µm |
| Crystal length | `L_mm` (in `pmFactor_sfg`) | mm |

The speed of light is taken as c = 2.99792458 × 10⁸ m/s
(0.299792458 µm/fs).

## Method arguments

Every dispersion method of a crystal takes the same four arguments:

```
x.n(wl_um, angle_rad, T_degC, pol='o')
```

Media in `ndispers.media.glasses` take only `(wl_um, T_degC)`, with no angle
and no polarization.

The two packages are split by **optical isotropy**, which is what decides the
call signature — not by whether a medium is crystalline, and not by
centrosymmetry. An isotropic medium has one refractive index with no direction
or polarization dependence, so there is no angle to pass. `CaF₂` is a crystal
(cubic, m3̄m) and sits with the glasses for that reason. α-BBO and calcite are
centrosymmetric as well, but they are uniaxial and therefore birefringent, so
they are crystals here.

- Scalar input returns a Python float. numpy arrays are accepted for any
  argument and broadcast together; plain Python lists and tuples are
  converted to arrays.
- `pol` defaults to `'o'` (ordinary ray), with two exceptions: the
  walk-off methods `woa_theta`/`woa_phi` default to `'e'` (the ordinary
  ray of a uniaxial crystal has no walk-off), and every method of `SLN`
  defaults to `'e'` because its Sellmeier equation exists only for the
  extraordinary ray.
- All derivative quantities (`GD`, `GV`, `ng`, `GVD`, `TOD`, `dndT`,
  `dndT2`, walk-off angles) are obtained by differentiating the medium's
  symbolic Sellmeier expression (`x.n_expr(pol)`), not by numerical
  finite differences.

## Which angle is the angle argument?

θ is the polar angle of the wavevector from the z principal axis; φ is
the azimuthal angle from the x axis.

- **Uniaxial crystals** (BBO, CLBO, KDP, ...): the angle argument is θ.
  Nothing depends on φ.
- **Biaxial crystals** (KTP, LBO): one class per principal dielectric
  plane; one angle is fixed by the plane and the other is the argument:

| Plane class | Fixed | Angle argument | o-wave polarized along | e-wave polarized in |
|---|---|---|---|---|
| `*_xy` | θ = π/2 | φ | z axis | xy plane |
| `*_yz` | φ = π/2 | θ | x axis | yz plane |
| `*_zx` | φ = π/2 | θ | y axis | zx plane |

The instance attributes tell you which case you have: `x.theta_rad` and
`x.phi_rad` are each either a fixed value in radians, the string
`'var'` (this is the angle you pass), or `'arb'` (the result does not
depend on it). The fixed angle is supplied internally; you never pass it.

## Temperature

The reference temperature — the T at which the temperature correction
vanishes — is a property of each literature source, **not** of the
package. As implemented:

| Media | Temperature term vanishes at |
|---|---|
| β-BBO (all four sources), CLBO, LBO (all sources), KTP | 20 °C |
| RBBF, fused silica, CaF₂ | 24 °C |
| SLN, SLT (Gayer 2008 form) | 24.5 °C |
| LB4 | 25 °C |
| KDP | 24.8 °C |
| KBBF | **no thermo-optic coefficients have been reported in the literature**, so the temperature dependence cannot be computed: the implemented coefficients are zero, n is T-independent and `dndT` is 0 |
| α-BBO, Calcite | no temperature term; `T_degC` is accepted and ignored, `dndT` is 0 |

A medium whose thermo-optic coefficients are all zero emits a one-time
`ndispers.TemperatureWarning` on first use: its methods accept `T_degC` but
ignore it, so a temperature sweep returns constant results.

Not every material property is known for every crystal. A method that
exists on all media is no guarantee that the underlying coefficients
have ever been measured — when they have not, the class implements
zeros. If a temperature dependence matters to your calculation, check
that `x.dndT(...)` is nonzero (or read the class docstring) before
relying on it.

## Dispersion quantities

With λ the vacuum wavelength and c the speed of light, the methods
implement:

| Method | Definition |
|---|---|
| `dn_wl` | dn/dλ |
| `ng` | n − λ·(dn/dλ) |
| `GD` | ng/c per unit length |
| `GV` | c/ng |
| `GVD` | (λ³/2πc²)·(d²n/dλ²) per unit length |
| `TOD` | −(λ⁴/4π²c³)·(3·d²n/dλ² + λ·d³n/dλ³) per unit length |
| `GVM` | GD(wl1) − GD(wl2), fs/mm; positive = wl1 arrives later |

Positive `GVD` is normal dispersion (longer wavelengths travel faster).

## Walk-off

`woa_theta` returns arctan(−(1/n)·∂n/∂θ), the angle between the
Poynting vector (energy-flow direction) and the wavevector, in radians;
`woa_phi` is the azimuthal analogue. For a uniaxial crystal the
ordinary ray has zero walk-off.

## Phase matching

For sum-frequency generation with input wavelengths λ₁, λ₂ and
1/λ₃ = 1/λ₁ + 1/λ₂, the collinear wavevector mismatch is

```
Δk = k₃ − k₂ − k₁ = 2π (n₃/λ₃ − n₂/λ₂ − n₁/λ₁)      [rad/µm]
```

so Δk = 0 is phase matching. Second-harmonic generation is the special
case λ₁ = λ₂. The conventions are stated here for SFG; the parametric
processes below mirror them.

### Sum-frequency generation and SHG

- The three-letter keys of `pmAngles_sfg` give the polarizations of
  (wave 1, wave 2, sum wave): `'ooe'` and `'eeo'` are Type I,
  `'oee'`/`'eoe'` and `'eoo'`/`'oeo'` are Type II. (`'ooe'`-type
  matching occurs in negative uniaxial crystals, `'eeo'`-type in
  positive ones.)
- `pmAngles_sfg` searches 0–90° and returns all solutions; an empty
  list means no phase-matching angle exists for that combination.
  `tol_deg` is the absolute tolerance of the returned angle in degrees.
- `pmFactor_sfg` returns sinc²(Δk·L/2), the relative conversion
  efficiency of a crystal of length `L_mm`.
- `acceptance_sfg(..., L_mm, param)` returns the acceptance width — the FWHM
  of sinc²(Δk·L/2) along one variable with the rest fixed: `param='wl'`
  (degenerate-SHG fundamental, wl1 = wl2 swept together), `'wl1'`/`'wl2'`
  (mix acceptance, one input wave alone), `'theta'` (internal angle) or
  `'T'`. Units: µm, rad, K respectively; `numpy.inf` along a noncritical
  direction. Evaluate at a phase-matched point.

- `qpm_period_sfg(wl1, wl2, angle_rad, T_degC, pol1, pol2, pol3, order=1)`
  returns the quasi-phase-matching period Λ = 2πm/|Δk| in µm for the same
  arguments; for PPLN/PPKTP use `angle_rad = pi/2` and `('e', 'e', 'e')`
  (all waves along the polar axis, d33).

### Difference-frequency generation, OPA and OPO

The parametric processes are the same interaction entered from the pump
side: k_p = k_s + k_i with 1/λ_p = 1/λ_s + 1/λ_i. Every `_dfg` method is the
`_sfg` one evaluated at (signal, idler), and the **three-letter keys and
polarization arguments read (signal, idler, pump)** — the order SNLO uses
(red1, red2, blue) — so `'ooe'` is a Type I OPA and `'oee'`/`'eoe'` Type II.

- `wl_idler(wl_p, wl_s)` is the idler; the signal must be longer than the pump.
- `pmAngles_dfg(wl_p, wl_s, T_degC)` returns the `pmAngles_sfg` dictionary
  with `'wl_i'` in place of `'wl3'`; `dk_dfg`, `pmFactor_dfg`,
  `qpm_period_dfg`, `d_dfg`, `deff_dfg` mirror their SFG counterparts.
- `tuning_dfg(wl_p, angle_rad, T_degC, pol_s, pol_i, pol_p)`
  solves Δk(λ_s) = 0 at a fixed angle and temperature and returns every
  (λ_s, λ_i) pair with λ_s ≤ λ_i — one point of an OPO tuning curve. The
  Type II branch with the other wave as signal is the triple with pol_s
  and pol_i exchanged. The search stops at the medium's Sellmeier validity
  edge by default; pass `wl_i_max=` to extrapolate further on purpose. With `qpm_period=Λ` (and
  `qpm_order`) it solves |Δk| = 2πm/Λ instead — the temperature tuning of a
  periodically poled crystal: `SLN().tuning_dfg(1.064, pi/2, T, 'e', 'e', 'e',
  qpm_period=30.85)` moves the signal from 1550 nm at 25 °C to 1640 nm at 200 °C.
- Bandwidths differ from SFG: for an OPA the pump is fixed and the idler
  follows the signal, so the gain bandwidth is set to first order by the
  signal–idler group-velocity mismatch. The phase-matching app computes it
  that way in its DFG mode.

## Second-order nonlinearity

Crystals whose point group allows a second-order nonlinearity inherit a
class from `ndispers.groups` (3m, 32, 4̄2m, 4mm or mm2) and carry one
reference measurement per independent tensor component in `_d_ref`:

```
_d_ref = {"d22": (2.2, 1.064, 1.064)}     # pm/V, for SFG of wl1 and wl2 in µm
```

- `d_sfg("d22", wl1, wl2, T_degC)` (also spelled `d22_sfg(wl1, wl2, T_degC)`)
  returns the component scaled to other wavelengths by **Miller's rule**,
  d_ijk ∝ χ_ii(ω₃) χ_jj(ω₁) χ_kk(ω₂) with χ_ii = n_i² − 1 the principal
  susceptibilities of the crystal frame. It reproduces the reference value
  exactly at the reference wavelengths. Miller's rule is an empirical,
  single-parameter estimate: good for BBO, KDP and LiNbO₃, rough for KTP,
  untested for several others — each crystal's `_d_note` says what is known.
- `deff_sfg(wl1, wl2, theta_rad, phi_rad, T_degC, pol1, pol2, pol3)` is
  the effective coefficient, d_eff = ê₃·d:(ê₁ê₂), contracted with the
  unit **E**-field vectors of the three waves. Walk-off is included (the
  e-wave field is perpendicular to the Poynting vector, not to k). The
  overall sign is a convention; compare magnitudes. For a principal plane of
  a biaxial crystal pass `None` for the angle the plane fixes.
- Kleinman symmetry is assumed throughout (transparency range). Under it
  d14 of class 32 vanishes, and type-II d_eff of class 4mm vanishes.
- Component names (d31, d32, ...) follow each crystal's literature. For LBO
  the dielectric axes are x, y, z = a, −c, b, so the polar axis c is
  dielectric y and "d31" means d_{c a a}, not d_{z x x}; `Biax_mm2` handles
  the mapping through the crystal's `_mm2_axes`.

The derivation and the closed forms per point group are in
`docs/dev/deff_theory.tex` (and checked numerically in
`ndispers/tests/test_nonlinear.py`).
