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

Several tools already exist for this kind of calculation. The best-known and
most extensive is *SNLO* by Arlee V. Smith
([as-photonics.com](https://as-photonics.com/products/snlo/), a Windows GUI).
There are also web apps for refractive indices
([refractiveindex.info](https://refractiveindex.info/)) and phase matching
([toolbox.lightcon.com](http://toolbox.lightcon.com/)), and the iOS app
[iPhasematch](https://apps.apple.com/app/iphasematch/id492370060). Calculators
exist for desktop, web and mobile — but none of them can become a component of
*your* program.

*ndispers* provides the core of SNLO-style calculation — refractive index,
dispersion, phase matching and d_eff — **as a Python library**. Every method is
a plain function that takes and returns numpy arrays, so you can sweep
wavelength, angle and temperature at once and build the results directly into
your own numerical simulations (pulse propagation, OPO design, thermal
analysis) or Jupyter notebooks. And every coefficient is published with its
source: type `bbo?` and you see which table of which paper the numbers come
from, so a single class name in your Methods section makes the calculation
reproducible.

Being an ordinary, introspectable Python library has one more consequence:
**ndispers is agent-ready**. Everything an AI assistant needs is
self-describing — media enumerate with `dir()`, every docstring carries its
literature source and validity range, and
[llms.txt](https://ndispers.readthedocs.io/en/latest/llms.txt) condenses the
whole API, units and pitfalls into a single page. Hand that link to your
assistant and ask in plain language — *"find the phase-matching angle for
1064 nm type-I SHG in BBO and plot its temperature dependence"* — and it can
install, compute and plot without you writing a line of code. A GUI or web
calculator cannot be driven this way.

Nor is it limited to nonlinear crystals. Under the same interface it also
covers linear birefringent crystals that SNLO does not treat (α-BBO, calcite,
quartz, YVO₄, …) and optical glasses (fused silica, CaF₂, SF10, …), aiming to
handle ultrafast-optics calculations — waveplates, polarizers and dispersion
management included — in one package.

**Strengths at a glance**

- **Embeddable** — a library, not a GUI. All methods accept numpy arrays and
  run as pre-compiled numpy functions, so the only runtime dependency is numpy
  (no sympy, no first-call code-generation pause). Medium objects are
  picklable and can be passed straight to `multiprocessing` / `joblib`.
- **Agent-ready** — everything a GUI hides is exposed as introspectable
  Python, and [llms.txt](https://ndispers.readthedocs.io/en/latest/llms.txt)
  hands an AI assistant the whole API in one page: give it the link and ask in
  plain language. The calculator becomes conversational.
- **Transparent sourcing** — each medium's docstring states the Sellmeier
  reference and its wavelength/temperature validity range, and the `constants`
  property returns the coefficient values themselves. You can always trace
  where a number came from.
- **Multiple Sellmeier equations per crystal** — e.g. LBO ships Kato 1994,
  Kato–Kuroda 2018, Ghosh 1995 and manufacturer fits (Castech, Newlight) as
  separate classes, so you can compare disagreeing literature yourself.
- **Materials beyond SNLO** — linear birefringent crystals (α-BBO, calcite,
  quartz, sapphire, MgF₂, YVO₄) and isotropic media (fused silica, CaF₂,
  BaF₂, YAG, N-BK7, SF glasses, Si, Ge, …) with the same API.
- **Analytic temperature derivatives** — dn/dT and d²n/dT² are computed from
  differentiated expressions, not finite differences; useful for temperature
  tuning and thermal-lens analysis.
- **Nonlinear-optics toolkit** — phase mismatch Δk, direct phase-matching-angle
  solutions, the sinc² phase-matching factor, and d_eff wavelength-scaled by
  Miller's rule for every non-centrosymmetric crystal.
- **Easy to extend** — subclass the point-group base class and write one file
  of Sellmeier and d coefficients to add a crystal. Requests and contributions
  are welcome on [GitHub](https://github.com/akihiko-shimura/ndispers/issues).
- **Cross-platform** — `pip install ndispers` on Linux, macOS or Windows; runs
  on clusters and Colab alike.
- **Verified transcription** — coefficients are cross-checked against the
  literature by independent re-extraction and pinned by regression tests, so
  copy errors are caught by machinery, not luck.
