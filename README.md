# ndispers

**English** | [日本語](https://github.com/akihiko-shimura/ndispers/blob/main/README.ja.md)

[![PyPI](https://img.shields.io/pypi/v/ndispers)](https://pypi.org/project/ndispers/)
[![test](https://github.com/akihiko-shimura/ndispers/actions/workflows/test.yml/badge.svg)](https://github.com/akihiko-shimura/ndispers/actions/workflows/test.yml)
[![apps](https://github.com/akihiko-shimura/ndispers/actions/workflows/pages.yml/badge.svg)](https://akihiko-shimura.github.io/ndispers/)

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
- Effective nonlinear coefficient, d_eff, with the tensor components scaled to the working wavelengths by Miller's rule (every non-centrosymmetric crystal in the package)

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

The effective nonlinear coefficient at that angle, for the φ = 90° cut:

```python
>>> bbo.deff_sfg(1.064, 1.064, np.radians(22.88), np.radians(90), 25, 'o', 'o', 'e')
-1.9937...  # pm/V (the overall sign is a convention); d22 = 2.2 pm/V at 1.064 µm SHG
            # (Shoji et al. 1999), walk-off included
>>> bbo.d_sfg("d22", 0.8, 0.8, 25)
2.3300...   # d22 for 800 nm SHG by Miller's rule
```

Every non-centrosymmetric crystal has this (KDP, CLBO, LBO, KTP, ...); `help(bbo.deff_sfg)` and the *Conventions* page describe the sign and angle conventions, and each crystal's `_d_note` says where its coefficients come from and how far Miller scaling has been tested for it.

## Available media

Every medium is a class in `nd.media.crystals` or `nd.media.glasses`; the
[media catalog](https://ndispers.readthedocs.io/en/latest/api/crystals/) lists
the class names with each one's Sellmeier equation, validity range and
references. Several crystals come in more than one parameterisation, named
after the literature or vendor the coefficients come from — pick the source
you want to rely on.

### Nonlinear optical crystals

Non-centrosymmetric: phase matching, acceptance bandwidths and d_eff are available.

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
| Zinc germanium phosphide | ZGP | ZnGeP₂ | 4̄2m | positive uniaxial, mid-infrared |
| Silver thiogallate | AGS | AgGaS₂ | 4̄2m | negative uniaxial, mid-infrared |
| Silver gallium selenide | AGSe | AgGaSe₂ | 4̄2m | negative uniaxial, mid-infrared |
| α-Quartz | — | SiO₂ | 32 | positive uniaxial |
| Lithium niobate, 5% MgO-doped congruent | MgO:LN | MgO:LiNbO₃ | 3m | negative uniaxial, both rays |
| Lithium niobate, 1% MgO-doped stoichiometric | MgO:SLN | MgO:LiNbO₃ | 3m | negative uniaxial, e-ray only |
| Lithium tantalate, 1% MgO-doped stoichiometric | MgO:SLT | MgO:LiTaO₃ | 3m | negative uniaxial |

### Birefringent optical crystals

Centrosymmetric, so no second-order nonlinearity; dispersion, walk-off and
thermo-optics for windows, polarizers and compensators.

| Material | Abbreviation | Formula | Point group | Optical class |
|---|---|---|---|---|
| α-Barium borate | α-BBO | α-BaB₂O₄ | 3̄m | negative uniaxial |
| Calcite | — | CaCO₃ | 3̄m | negative uniaxial |
| Sapphire | — | α-Al₂O₃ | 3̄m | negative uniaxial |

### Optically isotropic media

One refractive index, no angle or polarization argument: methods take `(wl_um, T_degC)`.

| Material | Formula |
|---|---|
| Fused silica | SiO₂ (amorphous) |
| Calcium fluoride | CaF₂ (cubic, m3̄m) |
| Yttrium aluminium garnet (YAG) | Y₃Al₅O₁₂ (cubic, m3̄m) |
| N-BK7 (SCHOTT) | borosilicate crown glass |

For biaxial crystals one class per principal dielectric plane is provided; the
angle argument is the one that varies in that plane (φ in xy, θ in yz and zx),
and an instance's `theta_rad` / `phi_rad` attributes say which. Media whose
Sellmeier equation carries no temperature term (α-BBO, calcite, sapphire,
quartz, 5% MgO:LiNbO₃, BiBO, the mid-infrared crystals, YAG, N-BK7) still take the temperature argument for a uniform
signature and ignore it; their `dndT` returns 0.

## Parallel processing

Medium objects are picklable, including after use, so they can be passed directly to `multiprocessing`, `concurrent.futures` or `joblib` workers, and stored with `pickle`/`shelve`. Lambdified dispersion functions are cached per class, so constructing many instances of the same crystal is cheap.

## Documentation

- **[Browser apps](https://akihiko-shimura.github.io/ndispers/)** — a refractive-index explorer and a phase-matching calculator, running client-side with nothing to install.
- [Tutorial notebook](examples/basic_usage.ipynb) — the basic workflow, from inspecting a crystal's Sellmeier equation to phase-matching curves, viewable directly on GitHub.
- [ndispers.readthedocs.io](https://ndispers.readthedocs.io/en/latest/) — conventions (units, angles, sign conventions) and the [media catalog](https://ndispers.readthedocs.io/en/latest/api/crystals/) of every crystal and glass with its formula, validity range and references.

## Development

```
git clone https://github.com/akihiko-shimura/ndispers.git
cd ndispers
uv sync
uv run pytest
```

To run the tutorial notebook, add the `notebook` group (matplotlib, IPython, JupyterLab, marimo):

```
uv sync --group notebook
uv run jupyter lab examples/basic_usage.ipynb
```

Releases are published to PyPI by pushing a version tag; see [docs/RELEASING.md](docs/RELEASING.md).

## License

MIT — see [LICENSE](LICENSE).
