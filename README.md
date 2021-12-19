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

This package is a general-purpose tool intended for computing phase-matching angles and bandwidths in designing nonlinear frequency mixing.

## Installation

In terminal, just
```zsh
pip install ndispers
```

## Example




## Call for Contributions

With Sellmeier-equation and thermo-optic coefficient provided, this package can be extended to include new crystals or glasses. 
