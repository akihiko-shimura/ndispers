<img src="https://raw.githubusercontent.com/akihiko-shimura/ndispers/main/docs/assets/logo.svg" alt="" width="132" align="right">

# ndispers

**English** | [日本語](https://github.com/akihiko-shimura/ndispers/blob/main/README.ja.md)

[![PyPI](https://img.shields.io/pypi/v/ndispers)](https://pypi.org/project/ndispers/)
[![test](https://github.com/akihiko-shimura/ndispers/actions/workflows/test.yml/badge.svg)](https://github.com/akihiko-shimura/ndispers/actions/workflows/test.yml)
[![apps](https://github.com/akihiko-shimura/ndispers/actions/workflows/pages.yml/badge.svg)](https://akihiko-shimura.github.io/ndispers/)

*ndispers* is a Python package for refractive-index dispersion of crystals and
glasses used in nonlinear and ultrafast optics. It implements Sellmeier
equations and thermo-optic coefficients reported in the literature.

It computes

- refractive index
- group delay, group velocity, group index
- group-velocity dispersion, third-order dispersion
- group-velocity mismatch between two waves
- walk-off angles
- dn/dT, d²n/dT²

as functions of

1. wavelength
2. polar (θ) or azimuthal (φ) angle of the wavevector with respect to the dielectric principal axes
3. temperature
4. polarization (ordinary or extraordinary ray)

For non-centrosymmetric crystals it further computes phase mismatch Δk,
phase-matching angles, the sinc²(Δk·L/2) factor, acceptance widths
(spectral, angular, thermal), quasi-phase-matching periods, OPO tuning
curves, and the effective nonlinear coefficient d_eff with tensor components
scaled to the working wavelengths by Miller's rule.

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
  coefficients ([AGENTS.md](AGENTS.md)). Requests and contributions:
  [GitHub](https://github.com/akihiko-shimura/ndispers/issues).
- Installation is `pip install ndispers` on any platform, Python ≥ 3.10.

## Installation

```
pip install ndispers
```

Requires Python 3.10 or later; the only dependency is numpy. The dispersion
formulas are written as sympy expressions, but the evaluated functions are
generated ahead of time and shipped. To work with the expressions themselves
(`n_expr` and related methods), or to evaluate your own subclass of a medium,
install sympy with `pip install ndispers[sym]`.

## Quick start

The [tutorial notebook](examples/basic_usage.ipynb) covers the workflow from
inspecting a Sellmeier equation to phase-matching curves, with plots. The
essentials:

```python
>>> import ndispers as nd
>>> bbo = nd.media.crystals.BetaBBO_Eimerl1987()
>>> bbo.n(0.532, 0, 25, pol='o')
1.674884049110459
>>> bbo.n(0.532, 3.1416/2, 25, pol='e')
1.5554658787539917
```

The four arguments are

1. wavelength (µm),
2. θ angle (rad),
3. temperature (°C),
4. polarization (`pol='o'` or `'e'`; default `'o'`).

Scalar input returns a float; numpy-array input broadcasts and returns an
array. The other dispersion methods (`GD`, `GV`, `ng`, `GVD`, `TOD`, `GVM`,
walk-off, `dndT`) and the phase-matching methods (`pmAngles_sfg`,
`acceptance_sfg`, `qpm_period_sfg`, `tuning_dfg`, `deff_sfg`, …) take the
same style of arguments. `help(bbo)` prints the material data and
references.

Worked examples — dispersion sweeps, phase-matching angles, OPO tuning,
d_eff — are in the [documentation](https://ndispers.readthedocs.io/en/latest/)
and the [tutorial notebook](examples/basic_usage.ipynb).

## Available media

Nonlinear crystals: β-BBO, LBO, KTP, BiBO, CLBO, KDP, DKDP, KBBF, RBBF, LB4,
LiIO₃, ZGP, AgGaS₂, AgGaSe₂, quartz, MgO-doped LiNbO₃ (congruent and
stoichiometric) and MgO:LiTaO₃. Linear birefringent crystals: α-BBO, calcite,
sapphire, MgF₂, YVO₄. Isotropic media: fused silica, CaF₂, LiF, BaF₂, YAG,
N-BK7, SF10/11/57, ZnSe, ZnS, Si, Ge, diamond.

Each medium is a class in `nd.media.crystals` or `nd.media.glasses`; several
crystals have more than one parameterisation, one class per source. The
[media catalog](https://ndispers.readthedocs.io/en/latest/api/crystals/)
lists every class with its Sellmeier equation, validity range and references;
`nd.catalog()` returns the same census as data.

## Documentation

- [Browser apps](https://akihiko-shimura.github.io/ndispers/) — a refractive-index explorer and a phase-matching calculator, running client-side.
- [Tutorial notebook](examples/basic_usage.ipynb) — the basic workflow, viewable on GitHub.
- [Validation](https://ndispers.readthedocs.io/en/latest/validation/) — comparisons against the literature, with numbers, sources and caveats.
- [ndispers.readthedocs.io](https://ndispers.readthedocs.io/en/latest/) — conventions (units, angles, signs) and the [media catalog](https://ndispers.readthedocs.io/en/latest/api/crystals/).
- [llms.txt](https://ndispers.readthedocs.io/en/latest/llms.txt) / [llms-full.txt](https://ndispers.readthedocs.io/en/latest/llms-full.txt) — the API reference for AI assistants.

## Development

```
git clone https://github.com/akihiko-shimura/ndispers.git
cd ndispers
uv sync
uv run pytest
```

To run the tutorial notebook, add the `notebook` group (matplotlib, IPython,
JupyterLab, marimo):

```
uv sync --group notebook
uv run jupyter lab examples/basic_usage.ipynb
```

Releases are published to PyPI by pushing a version tag; see
[docs/RELEASING.md](docs/RELEASING.md). [AGENTS.md](AGENTS.md) describes how
to add a crystal.

## License

MIT — see [LICENSE](LICENSE).
