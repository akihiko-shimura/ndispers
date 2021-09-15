# ndispers
*ndispers* is a Python3 package for calculating refractive index dispersion of various media such as crystals and glasses commonly used in the field of nonlinear/ultrafast optics.
It allows you to calculate refractive index and its derivative quantities, such as group velocity, group velocity dispersion, group index, third-order dispersion, walkoff angles, as a function of 1) wavelength of light, 2) polar (theta) or azimuthal (phi) angles of wave vector with respect to dielectric principal axes of anisotropic crystals, 3) crystal temperature, and 4) polarization of light (ordinary- or extraordinary-ray).
Since this package provides those functions with the methods of an medium object, you can easily implement it to other numerical simulations in nonlinear/ultrafast optics.

Ndisperse uses [SymPy](https://github.com/sympy/sympy) to simply implement Sellmeier's equations and its derivatives. Although the first call of a higher-oder derivatives (e.g., TOD) is rather slow, this drawback is avoided by caching the lambdified function from the second call onwards.
