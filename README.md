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

as a function of
1. Wavelength of light
2. Polar (theta) or azimuthal (phi) angles of wave vector with respect to dielectric principal axes of anisotropic crystals
3. Temperature
4. Polarization of light (ordinary- or extraordinary-ray)


## Installation

In terminal,
```zsh
pip install ndispers
```

## A simple example

Firstly, make an object of beta-BBO crystal.

```python
>>> import ndispers as nd
>>> bbo = nd.media.crystals.BetaBBO_Eimerl1987()
```

For material information, 

```
>>> bbo.help
beta-BBO (beta-Ba B_2 O_4) crystal

    - Point group : 3m
    - Crystal system : Trigonal
    - Dielectic principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Tranparency range : 1.9 to 2.6 um

    Dispersion formula for refractive index
    ---------------------------------------
    n(wl) = sqrt(A_i + B_i/(wl**2 - C_i) - D_i * wl**2)  for i = o, e
    
    Validity range
    ---------------
    0.22 to 1.06 um

    Ref
    ----
    Eimerl, David, et al. "Optical, mechanical, and thermal properties of barium borate." Journal of applied physics 62.5 (1987): 1968-1983.
    Nikogosyan, D. N. "Beta barium borate (BBO)." Applied Physics A 52.6 (1991): 359-368.
```

To compute refractive indices, use a method of the `bbo` instance,

```python
>>> bbo.n(0.532, 0, 25, pol='o')
array(1.67488405)
>>> bbo.n(0.532, 3.1416/2, 25, pol='e')
array(1.55546588)
```

where the four arguments are
1. wavelength (in micrometer), 
2. theta angle (in radians),
3. temperature (in deg.C), 
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

Complete [documentation](https://ndispers.readthedocs.io/en/latest/) is under construction.
