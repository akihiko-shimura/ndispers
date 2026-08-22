# ndispers

[![PyPI](https://img.shields.io/pypi/v/ndispers)](https://pypi.org/project/ndispers/)
[![test](https://github.com/akihiko-shimura/ndispers/actions/workflows/test.yml/badge.svg)](https://github.com/akihiko-shimura/ndispers/actions/workflows/test.yml)

*ndispers* is a Python package for calculating refractive index dispersion of various crystals and glasses used in the field of nonlinear/ultrafast optics. It is based on Sellmeier equations and thermo-optic coefficients (dn/dT) reported in literature.

You can easily compute

- Refractive index
- Group delay
- Group velocity
- Group index
- Group velocity dispersion
- Third-order dispersion
- Walk-off angles
- dn/dT
- d²n/dT²

as a function of

1. Wavelength of light
2. Polar (theta) or azimuthal (phi) angle of the wavevector with respect to the dielectric principal axes of an anisotropic crystal
3. Temperature of the medium
4. Polarization of light (ordinary or extraordinary ray)

The crystals also have nonlinear-optics methods:

- Phase mismatch, dk
- Phase-matching angles
- Phase-matching factor, sinc²(dk·L/2)
- Effective nonlinear coefficient, deff

## Installation

Requires Python 3.9 or later.

```
pip install ndispers
```

## Quick start

For a guided tour with plots — from inspecting a crystal's Sellmeier equation to phase-matching curves — see the [tutorial notebook](examples/basic_usage.ipynb), readable directly on GitHub. The essentials:

Make an object of a β-BBO crystal:

```python
>>> import ndispers as nd
>>> bbo = nd.media.crystals.BetaBBO_Eimerl1987()
```

To compute a refractive index:

```python
>>> bbo.n(0.532, 0, 25, pol='o')
1.674884049110459
>>> bbo.n(0.532, 3.1416/2, 25, pol='e')
1.5554658787539917
```

where the four arguments are, respectively,

1. wavelength (in micrometer),
2. theta angle (in radian),
3. temperature (in degree Celsius),
4. polarization (`pol='o'` or `'e'`, ordinary or extraordinary ray; default `'o'`).

At theta = 0 (propagation along the optic axis) the o-ray and e-ray indices coincide.
Scalar input returns a float. Each argument also accepts a numpy array, returning an array of the same shape:

```python
>>> import numpy as np
>>> bbo.n(np.arange(0.2, 1.5, 0.2), 0, 25, pol='o')
array([1.89001202, 1.69328828, 1.66985875, 1.6612891 , 1.65633946,
       1.65252664, 1.64903624])
```

The other dispersion methods (`GD`, `GV`, `ng`, `GVD`, `TOD`, `woa_theta`, `woa_phi`, `dndT`, `dndT2`) take the same arguments:

```python
>>> bbo.GVD(0.8, 0, 25, pol='o')   # fs^2/mm
71.86403019943364
```

Optically isotropic media — the glasses and the cubic crystal CaF₂ — take only wavelength and temperature, since there is no angle or polarization to specify:

```python
>>> silica = nd.media.glasses.FusedSilica()
>>> silica.n(1.064, 20)
1.4495857898590634
```

To look into the material information, its Sellmeier equation and the literature the coefficients come from, use `help(bbo)` — or simply `bbo?` in IPython/Jupyter. `bbo.constants` returns the coefficient values themselves.

## Phase matching

Phase-matching angles for sum-frequency generation are solved directly. For Type-I SHG of 1064 nm in β-BBO at 25 °C:

```python
>>> bbo.pmAngles_sfg(1.064, 1.064, 25, deg=True)
{'wl3': 0.532,
 'ooe': {'theta': [22.884169498625802], 'phi': None},
 'eeo': {'theta': [], 'phi': None},
 'oee': {'theta': [32.56045545648089], 'phi': None},
 'eoe': {'theta': [32.56045545648089], 'phi': None},
 'eoo': {'theta': [], 'phi': None},
 'oeo': {'theta': [], 'phi': None}}
```

`dk_sfg` and `pmFactor_sfg` give the phase mismatch and the sinc² phase-matching factor for arbitrary angles.

## Available media

Several crystals are provided in more than one parameterisation, named after the literature (or vendor) the coefficients come from — pick the source you want to rely on.

| Medium | Classes (`nd.media.crystals.*`) |
|---|---|
| β-BBO | `BetaBBO_Eimerl1987`, `BetaBBO_Ghosh1995`, `BetaBBO_KK2010`, `BetaBBO_Tamosauskas2018` |
| LBO (xy / yz / zx principal planes) | `LBO_Castech_*`, `LBO_Ghosh1995_*`, `LBO_KK1994_*`, `LBO_KK2018_*`, `LBO_Newlight_*` |
| KTP (xy / yz / zx principal planes) | `KTP_xy`, `KTP_yz`, `KTP_zx` |
| CLBO | `CLBO` |
| KDP | `KDP` |
| KBBF | `KBBF` |
| RBBF | `RBBF` |
| LB4 | `LB4` |
| α-BBO | `AlphaBBO` |
| Calcite | `Calcite` |
| 1% MgO-doped stoichiometric LiNbO₃ (e-ray) | `SLN` |
| 1% MgO-doped stoichiometric LiTaO₃ | `SLT` |
| Fused silica, CaF₂ (optically isotropic) | `nd.media.glasses.FusedSilica`, `nd.media.glasses.CaF2` |

For biaxial crystals (LBO, KTP), one class per principal dielectric plane is provided; the angle argument is the one that varies in that plane (phi in xy, theta in yz and zx). The `theta_rad` and `phi_rad` attributes of an instance tell which one it is.

Media whose Sellmeier equation carries no temperature term (α-BBO, Calcite) still take the temperature argument for a uniform signature and simply ignore it — their `dndT` returns 0.

## Parallel processing

Medium objects are picklable, including after use, so they can be passed directly to `multiprocessing`, `concurrent.futures` or `joblib` workers, and stored with `pickle`/`shelve`. Lambdified dispersion functions are cached per class, so constructing many instances of the same crystal is cheap.

## Documentation

- [Tutorial notebook](examples/basic_usage.ipynb) — the basic workflow, from inspecting a crystal's Sellmeier equation to phase-matching curves, viewable directly on GitHub.
- [ndispers.readthedocs.io](https://ndispers.readthedocs.io/en/latest/) — conventions (units, angles, sign conventions) and the [media catalog](https://ndispers.readthedocs.io/en/latest/api/crystals/) of every crystal and glass with its formula, validity range and references.

## Development

```
git clone https://github.com/akihiko-shimura/ndispers.git
cd ndispers
uv sync
uv run pytest
```

Releases are published to PyPI by pushing a version tag; see [docs/RELEASING.md](docs/RELEASING.md).

## License

MIT — see [LICENSE](LICENSE).
