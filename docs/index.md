# ndispers

*ndispers* is a Python package for calculating refractive index dispersion of
various crystals and glasses used in nonlinear/ultrafast optics. It is based on
Sellmeier equations $n(\lambda)$ and thermo-optic coefficients (*dn/dT*)
reported in literature.

## Installation

*ndispers* works on Python 3.10 or higher:

```
pip install ndispers
```

## Quick start

For a guided tour with plots, see the
[tutorial notebook](https://github.com/akihiko-shimura/ndispers/blob/main/examples/basic_usage.ipynb),
readable directly on GitHub. The essentials:

```python
>>> import ndispers as nd
>>> bbo = nd.media.crystals.BetaBBO_Eimerl1987()
>>> bbo.n(0.532, 0, 25, pol='o')
1.674884049110459
```

The four arguments are wavelength (µm), theta angle (rad), temperature (°C) and
polarization (`pol='o'` or `'e'`, default `'o'`). Scalar input returns a float;
numpy-array input returns an array of the same shape:

```python
>>> import numpy as np
>>> bbo.n(np.arange(0.2, 1.5, 0.2), 0, 25, pol='o')
array([1.89001202, 1.69328828, 1.66985875, 1.6612891 , 1.65633946,
       1.65252664, 1.64903624])
```

The other dispersion methods — `GD`, `GV`, `ng`, `GVD`, `TOD`, `woa_theta`,
`woa_phi`, `dndT`, `dndT2` — take the same arguments. Phase-matching angles for
sum-frequency generation are solved directly:

```python
>>> bbo.pmAngles_sfg(1.064, 1.064, 25, deg=True)['ooe']
{'theta': [22.884169498625802], 'phi': None}
```

Use `help(bbo)` — or simply `bbo?` in IPython/Jupyter — to see any medium's
Sellmeier equation, validity range and the literature the coefficients come
from; the `constants` property returns the coefficient values themselves.

Medium objects are picklable, including after use, so they can be passed
directly to `multiprocessing` or `joblib` workers.

## Phase matching

`pmAngles_sfg` solves for the phase-matching angles of sum-frequency
generation. For SHG of 1064 nm in β-BBO at 25 °C:

```python
>>> bbo.pmAngles_sfg(1.064, 1.064, 25, deg=True)
{'wl3': 0.532,
 'ooe': {'theta': [22.884169498625802], 'phi': None},   # Type I
 'eeo': {'theta': [], 'phi': None},
 'oee': {'theta': [32.56045545648089], 'phi': None},    # Type II
 'eoe': {'theta': [32.56045545648089], 'phi': None},
 'eoo': {'theta': [], 'phi': None},
 'oeo': {'theta': [], 'phi': None}}
```

`dk_sfg` and `pmFactor_sfg` give the phase mismatch and the sinc² factor at
arbitrary angles; `acceptance_sfg` gives the FWHM of the sinc² curve along
wavelength, angle or temperature; `qpm_period_sfg` gives the
quasi-phase-matching period. For 1064 nm SHG in 5% MgO:LiNbO₃ along x with
all waves extraordinary (d₃₃),
`MgOLN_Zelmon1997().qpm_period_sfg(1.064, 1.064, np.pi/2, 20, 'e', 'e', 'e')`
returns 7.0 µm.

Difference-frequency generation and optical parametric processes enter the
same interaction from the pump side: `pmAngles_dfg(wl_p, wl_s, T)` gives the
angles and the idler wavelength; `tuning_dfg(wl_p, angle, T, pol_s, pol_i,
pol_p)` gives the signal/idler pairs that phase-match at a fixed angle — one
point of an OPO tuning curve — and accepts `qpm_period=` for the temperature
tuning of a periodically poled crystal. `dk_dfg`, `pmFactor_dfg`,
`qpm_period_dfg` and `deff_dfg` mirror their SFG counterparts with waves read
as (signal, idler, pump). For x-cut KTP pumped at 1064 nm,
`KTP_zx().tuning_dfg(1.064, np.pi/2, 25, 'o', 'e', 'o')` returns the
noncritical pair 1571/3298 nm.

The effective nonlinear coefficient at a phase-matching angle, for the
φ = 90° cut:

```python
>>> bbo.deff_sfg(1.064, 1.064, np.radians(22.88), np.radians(90), 25, 'o', 'o', 'e')
-1.9937...  # pm/V; the overall sign is a convention. d22 = 2.2 pm/V at
            # 1.064 µm SHG (Shoji et al. 1999), walk-off included
>>> bbo.d_sfg("d22", 0.8, 0.8, 25)
2.3300...   # d22 for 800 nm SHG by Miller's rule
```

Sign and angle conventions are defined on the
[Conventions](conventions.md) page.
Each crystal's `_d_note` states the source of its d coefficients and the
extent to which Miller scaling has been tested for it.

## Available media

Every medium is a class in `nd.media.crystals` or `nd.media.glasses`. The
[media catalog](api/crystals.md)
lists each class with its Sellmeier equation, validity range and references;
`nd.catalog()` returns the same census as data. Where a crystal has several
parameterisations, each is a separate class named after its source.

### Nonlinear optical crystals

Non-centrosymmetric; phase matching, acceptance widths and d_eff are available.

| Material | Abbreviation | Formula | Point group | Optical class |
|---|---|---|---|---|
| β-Barium borate | β-BBO | β-BaB₂O₄ | 3m | negative uniaxial |
| Lithium triborate | LBO | LiB₃O₅ | mm2 | biaxial, three principal planes |
| Potassium titanyl phosphate | KTP | KTiOPO₄ | mm2 | biaxial, three principal planes |
| Bismuth triborate | BiBO | BiB₃O₆ | 2 | biaxial (monoclinic), three principal planes |
| Cesium lithium borate | CLBO | CsLiB₆O₁₀ | 4̄2m | negative uniaxial |
| Potassium dihydrogen phosphate | KDP | KH₂PO₄ | 4̄2m | negative uniaxial |
| Deuterated potassium dihydrogen phosphate | DKDP, KD*P | KD₂PO₄ | 4̄2m | negative uniaxial |
| Potassium beryllium fluoroborate | KBBF | KBe₂BO₃F₂ | 32 | negative uniaxial |
| Rubidium beryllium fluoroborate | RBBF | RbBe₂BO₃F₂ | 32 | negative uniaxial |
| Lithium tetraborate | LB4 (also LBT) | Li₂B₄O₇ | 4mm | negative uniaxial |
| Lithium iodate | — | LiIO₃ | 6 | negative uniaxial |
| Zinc germanium phosphide | ZGP | ZnGeP₂ | 4̄2m | positive uniaxial, mid-infrared |
| Silver thiogallate | AGS | AgGaS₂ | 4̄2m | negative uniaxial, mid-infrared |
| Silver gallium selenide | AGSe | AgGaSe₂ | 4̄2m | negative uniaxial, mid-infrared |
| α-Quartz | — | SiO₂ | 32 | positive uniaxial |
| Lithium niobate, 5% MgO-doped congruent | MgO:LN | MgO:LiNbO₃ | 3m | negative uniaxial, both rays |
| Lithium niobate, 1% MgO-doped stoichiometric | MgO:SLN | MgO:LiNbO₃ | 3m | negative uniaxial, e-ray only |
| Lithium tantalate, 1% MgO-doped stoichiometric | MgO:SLT | MgO:LiTaO₃ | 3m | negative uniaxial |

### Birefringent optical crystals

Centrosymmetric, hence no second-order nonlinearity; dispersion, walk-off and
thermo-optics for windows, polarizers and compensators.

| Material | Abbreviation | Formula | Point group | Optical class |
|---|---|---|---|---|
| α-Barium borate | α-BBO | α-BaB₂O₄ | 3̄m | negative uniaxial |
| Calcite | — | CaCO₃ | 3̄m | negative uniaxial |
| Sapphire | — | α-Al₂O₃ | 3̄m | negative uniaxial |
| Magnesium fluoride | — | MgF₂ | 4/mmm | positive uniaxial |
| Yttrium orthovanadate | YVO₄ | YVO₄ | 4/m | positive uniaxial |

### Optically isotropic media

One refractive index; methods take `(wl_um, T_degC)`.

| Material | Formula |
|---|---|
| Fused silica | SiO₂ (amorphous) |
| Calcium fluoride | CaF₂ (cubic, m3̄m) |
| Lithium fluoride | LiF (cubic, m3̄m) |
| Barium fluoride | BaF₂ (cubic, m3̄m) |
| Yttrium aluminium garnet (YAG) | Y₃Al₅O₁₂ (cubic, m3̄m) |
| N-BK7, SF10, SF11, SF57 (SCHOTT) | borosilicate crown and dense flint glasses |
| Zinc selenide, zinc sulfide | ZnSe, ZnS (cubic, 4̄3m; CVD grades) |
| Silicon, germanium | Si, Ge (cubic, m3̄m) |
| Diamond | C (cubic, m3̄m) |

Parameterisations of the same material are not interchangeable: each is
faithful to its own source, and the sources differ in range and accuracy. The
[validation page](validation.md)
records the comparisons.

For biaxial crystals one class per principal dielectric plane is provided.
The angle argument is the one that varies in that plane (φ in xy, θ in yz and
zx); an instance's `theta_rad` / `phi_rad` attributes state which. Media
whose Sellmeier equation carries no temperature term accept the temperature
argument for a uniform signature and ignore it; their `dndT` returns 0 and
the first use emits a `TemperatureWarning`. A wavelength outside a medium's
Sellmeier validity range emits a `ValidityWarning`; the extrapolated value is
still returned.

## Parallel processing

Medium objects are picklable, including after use, and can be passed to
`multiprocessing`, `concurrent.futures` or `joblib` workers, or stored with
`pickle`/`shelve`. Evaluated functions are cached per class, so constructing
many instances of the same crystal is cheap.

## Where things are

- **[Conventions](conventions.md)** — units, argument order, angle meaning per
  medium, reference temperatures, and the phase-matching sign conventions.
- **[Media catalog](api/crystals.md)** — every crystal and glass with its Sellmeier
  equation, validity range and references.
- **[Tutorial notebook](https://github.com/akihiko-shimura/ndispers/blob/main/examples/basic_usage.ipynb)**
  — worked examples with plots.

## Why ndispers

Existing tools for these calculations are applications: SNLO
([as-photonics.com](https://as-photonics.com/products/snlo/), Windows),
[refractiveindex.info](https://refractiveindex.info/) and
[toolbox.lightcon.com](http://toolbox.lightcon.com/) (web), and
[iPhasematch](https://apps.apple.com/app/iphasematch/id492370060) (iOS). An
application answers the query typed into it; it cannot be called from a
program.

ndispers is a Python library. It computes the quantities SNLO's Qmix function
reports — refractive indices, group velocities and dispersion, walk-off,
phase-matching conditions, acceptance widths, effective nonlinearity — as
functions that accept and return numpy arrays, so they can be evaluated over
wavelength, angle and temperature grids inside a larger calculation. Every
coefficient is traceable: each class docstring cites the paper, the validity
range, and separate sources for the Sellmeier equation, the thermo-optic
coefficients and the d tensor. Naming the class and the package version
specifies a calculation completely.

Because the library is introspectable and its data sources are stated, an AI
assistant can operate it: media are enumerable with `dir()` or `catalog()`,
docstrings carry the required metadata, and
[llms.txt](https://ndispers.readthedocs.io/en/latest/llms.txt) states the API
on one page.

The scope is wider than SNLO's: linear birefringent crystals (α-BBO, calcite,
quartz, MgF₂, YVO₄, sapphire) and isotropic media (glasses, CaF₂, Si, Ge, …)
are included under the same interface.

**Features**

- All methods accept numpy arrays and broadcast over them; the evaluated
  functions are pre-generated, so numpy is the only runtime dependency.
  Medium objects are picklable and work with `multiprocessing`.
- Temperature derivatives dn/dT and d²n/dT² are obtained by differentiating
  the symbolic index expressions, not by finite differences.
- Phase matching: Δk, phase-matching angles, sinc² factor, acceptance widths,
  QPM period, OPO tuning curves; d_eff with Miller wavelength scaling for
  every non-centrosymmetric crystal.
- Several parameterizations of the same crystal are separate classes (LBO has
  five); disagreement between sources is visible rather than hidden behind
  one default.
- Coefficients are transcribed from the cited papers and cross-checked by
  independent re-extraction; regression tests pin the computed values
  ([validation record](validation.md)).
- A new crystal is one file: a point-group base class plus Sellmeier and d
  coefficients ([AGENTS.md](https://github.com/akihiko-shimura/ndispers/blob/main/AGENTS.md)). Requests and contributions:
  [GitHub](https://github.com/akihiko-shimura/ndispers/issues).
- Installation is `pip install ndispers` on any platform, Python ≥ 3.10.
