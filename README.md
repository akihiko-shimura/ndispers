# ndispers
*ndispers* is a Python3 package for calculating refractive index dispersion of various crystals and glasses commonly used in the field of nonlinear/ultrafast optics. This is based on Sellmeier equartions and thermo-optics coefficients (dn/dT) reported in literatures.

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
3. Crystal temperature
4. Polarization of light (ordinary- or extraordinary-ray)

This package was created as a general-purpose tool for computing phase-matching angles and bandwidths etc. in designing nonlinear frequency mixing.

## Installation

In terminal, just
```zsh
pip install ndispers
```

## Example

Firstly, make an object of beta-BBO crystal.

```python
>>> import ndispers as nd
>>> bbo = nd.media.crystals.BetaBBO1995()
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
    n(wl) = sqrt(A_i + B_i/(1 - C_i/wl**2) + D_i/(1 - E_i/wl**2)) + dn/dT * (T -20)

    Thermo-optic coefficient
    -------------------------
    dn/dT = (G_i * R_i + H_i * R_i**2) / 2*n_i for i = o, e
    (R depends on wl)
    
    Validity range
    ---------------

    Ref
    ----
    Ghosh, Gorachand. "Temperature dispersion of refractive indices in β‐BaB2O4 and LiB3O5 crystals for nonlinear optical devices." Journal of applied physics 78.11 (1995): 6752-6760.
```

To compute refractive indices, use a method of the `bbo` instance,

```python
>>> bbo.n(0.532, 25, 0, pol='o')
array(1.6748653)
```

where the four arguments are
1. wavelength (in micro-meter), 
2. temperature (in deg.C), 
3. theta angle (in radians),
4. polarization (`pol='o' or 'e'`) as a keyword argument ('o' as default). 
Note that `pol='e'` corresponds to `pol='o'` when theta angle is 0 radians. 

Output value is generically of `numpy.ndarray` type. You can input an array to each argument, getting an output array of the same shape.

```python
>>> import numpy as np
>>> wl_ar = np.arange(0.2, 1.5, 0.1)
>>> wl_ar
array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1. , 1.1, 1.2, 1.3, 1.4])
>>> bbo.n(wl_ar, 25, 0, pol='o')
array([1.89369555, 1.73105003, 1.69347308, 1.67798949, 1.66981073,
       1.66476577, 1.6612739 , 1.65862123, 1.65644462, 1.65454248,
       1.65279501, 1.65112696, 1.64948873])
```


## Call for Contributions

With Sellmeier-equation and thermo-optic coefficient provided, this package can be extended to include new crystals or glasses. 
