# Conventions and definitions

The facts on this page hold across the whole package. Everything here is
stated as implemented in the code; per-crystal specifics (coefficients,
validity ranges, references) live in each class docstring and in the
[media catalog](api.md).

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

Glasses are isotropic and take only `(wl_um, T_degC)`, with no
polarization.

- Scalar input returns a Python float. numpy arrays are accepted for any
  argument and broadcast together (plain Python lists are not accepted).
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
| KDP | 24.8 °C in the implemented form, but the thermo-optic coefficients are set to zero — n is T-independent, `dndT` is 0 |
| KBBF | **no thermo-optic coefficients have been reported in the literature**, so the temperature dependence cannot be computed: the implemented coefficients are zero, n is T-independent and `dndT` is 0 |
| α-BBO, Calcite | no temperature term; `T_degC` is accepted and ignored, `dndT` is 0 |

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
case λ₁ = λ₂.

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
