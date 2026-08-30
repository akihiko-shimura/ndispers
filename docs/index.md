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
  ([validation record](https://ndispers.readthedocs.io/en/latest/validation/)).
- A new crystal is one file: a point-group base class plus Sellmeier and d
  coefficients ([AGENTS.md](https://github.com/akihiko-shimura/ndispers/blob/main/AGENTS.md)). Requests and contributions:
  [GitHub](https://github.com/akihiko-shimura/ndispers/issues).
- Installation is `pip install ndispers` on any platform, Python ≥ 3.10.
