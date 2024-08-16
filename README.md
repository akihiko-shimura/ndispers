# ndispers
*ndispers* is a Python package for calculating refractive index dispersion of various crystals and glasses used in the field of nonlinear/ultrafast optics. It is based on Sellmeier equartions and thermo-optic coefficients (dn/dT) reported in literature.

You can easily compute
- Refractive index
- Group delay
- Group velocity
- Group index
- Group velocity dispersion
- Third-order dispersion
- Walk-off angles
- dn/dT
- d^2n/dT^2

as a function of
1. Wavelength of light
2. Polar (theta) or azimuthal (phi) angles of wavevector with respect to dielectric principal axes of anisotropic crystals
3. Temperature of crystal
4. Polarization of light (ordinary- or extraordinary-ray)

The crystals have nonlinear-optics methods:
- Phase-mismacth, dk
- Phase-matching angles
- Phase-mathcing factor, sinc^2(dk*L/2)
- Effective nonlinear coefficient, deff


## Installation

In terminal,
```zsh
pip install ndispers
```

## Simple example

At first, make an object of β-BBO crystal.

```python
>>> import ndispers as nd
>>> bbo = nd.media.crystals.BetaBBO_Eimerl1987()
```

To look into the material information, 

```
>>> bbo.help
β-BBO (β-Ba B_2 O_4) crystal

    - Point group : 3m  (C_{3v})
    - Crystal system : Trigonal
    - Dielectic principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Tranparency range : 0.19 µm to 2.6 µm

    Sellmeier equation
    ------------------
    n(wl) = sqrt(A_i + B_i/(wl**2 - C_i) - D_i * wl**2) + dn/dT * (T - 20)  for i = o, e
    
    Validity range
    ---------------
    0.22 to 1.06 µm

    Ref
    ---
    - Eimerl, David, et al. "Optical, mechanical, and thermal properties of barium borate." Journal of applied physics 62.5 (1987): 1968-1983.
    - Nikogosyan, D. N. "Beta barium borate (BBO)." Applied Physics A 52.6 (1991): 359-368.

    Example
    -------
    >>> bbo = ndispers.media.crystals.BetaBBO_Eimerl1987()
    >>> bbo.n(0.6, 0.5*pi, 25, pol='e') # args: (wl_um, theta_rad, T_degC, pol)
```

To compute a refractive index,

```python
>>> bbo.n(0.532, 0, 25, pol='o')
array(1.67488405)
>>> bbo.n(0.532, 3.1416/2, 25, pol='e')
array(1.55546588)
```

where the four arguments are, respectively,
1. wavelength (in micrometer), 
2. theta angle (in radian),
3. temperature (in degree Celsius), 
4. polarization (`pol='o' or 'e'`, ordinary or extraordinary ray). 

Default is `pol='o'`. Note that `pol='e'` corresponds to `pol='o'` in index surface when theta angle is 0 radians. 
Output values are generically of `numpy.ndarray` type. You can input an array to each argument, getting an output array of the same shape, 

```python
>>> import numpy as np
>>> wl_ar = np.arange(0.2, 1.5, 0.2)
>>> wl_ar
array([0.2, 0.4, 0.6, 0.8, 1. , 1.2, 1.4])
>>> bbo.n(wl_ar, 0, 25, pol='o')
array([1.89625189, 1.692713, 1.66892613, 1.66039556, 1.65560236, 1.65199986, 1.64874414])
```

See [documentation](https://ndispers.readthedocs.io/en/latest/) for more features and examples.
