# ndispers

*ndispers* is a Python package for calculating refractive index dispersion of
various crystals and glasses used in nonlinear/ultrafast optics. It is based on
Sellmeier equations $n(\lambda)$ and thermo-optic coefficients (*dn/dT*)
reported in literature.

## Installation

*ndispers* works on Python 3.9 or higher:

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

## General overview

There are some softwares available for nonlinear optics. Probably the most
famous and extensive one is *SNLO*, by Arlee V. Smith
([as-photonics.com](https://as-photonics.com/products/snlo/), Windows only).
Other web applications exist to calculate refractive indices simply as a
function of wavelength ([refractiveindex.info](https://refractiveindex.info/))
or phase-matching conditions ([toolbox.lightcon.com](http://toolbox.lightcon.com/),
uniaxial crystals only). These apps are easy and quick to use, but most of them
are not open source; users cannot look into how a value was calculated or what
paper it is based on, and cannot extend the software for their own purposes.

This open-source Python project was created for those researchers, engineers
and students who want to study and employ *in depth* nonlinear/anisotropic
crystals and dispersive media, and is intended to be built into their numerical
simulation programs and Jupyter notebooks. You can request or contribute on
[GitHub](https://github.com/akihiko-shimura/ndispers/issues) to add new
crystals and methods.
