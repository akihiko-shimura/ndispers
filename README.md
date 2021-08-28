# ndispers
Ndispers is a Python3 package for calculating refractive index dispersion of various media such as crystals and glasses commonly used in the field of nonlinear/ultrafast  optics.
It allows you to calculate refractive index and its derivative quantities, such as group velocity, group velocity dispersion, group index, and third-order dispersion, as a function of wavelength of light, polar and azimuth angles of wave vector with respect to dielectric principal axes of anisotropic crystals.
Since this package provides those functions with the methods of an medium object, these calculations can be easily implemented to other numerical simulations in nonlinear/ultrafast optics.

Ndisperse uses SymPy (https://github.com/sympy/sympy) to simply implement Sellmeier's dispersion formula and its derivatives. Although the first call of a higher-oder derivatives (e.g., TOD) is rather slow, this drawback is avoided by caching the lambdified function from the second call onwards.
