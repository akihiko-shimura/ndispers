from ndispers._sym import sympy
from ndispers._baseclass import Medium, wl, T
from ndispers.helper import vars2

class LiF(Medium):
    """
    LiF (lithium fluoride) crystal

    - Point group : m3̄m  (Oh); space group Fm3̄m
    - Crystal system : cubic
    - Transparency range : 0.12 to 6.6 µm (Tropf et al., Table 20); the shortest-wavelength window material in common use

    Sellmeier equation
    ------------------
        n(wl)**2 = 1 + A1 * wl**2 / (wl**2 - B1**2) + A2 * wl**2 / (wl**2 - B2**2)

    Validity range
    ---------------
    0.1 to 10 µm

    Note
    ----
    The Sellmeier equation has no temperature term: T_degC is accepted for
    signature uniformity and ignored, and dndT returns 0. Tropf et al. (Table 20, from Feldman et al. 1978) list dn/dT = -16.0e-6 /K at 0.458 µm, -16.9e-6 at 1.15 µm and -14.5e-6 at 3.39 µm.

    Ref
    ---
    Sellmeier equation (as tabulated in):
      Tropf, W. J., Thomas, M. E., & Harris, T. J. (1995). Properties of crystals and glasses. In M. Bass (Ed.), Handbook of Optics, Vol. II, Chapter 33, Table 22. McGraw-Hill.
    Underlying measurements:
      Li, H. H. (1976). Refractive index of alkali halides and its wavelength and temperature derivatives. Journal of Physical and Chemical Reference Data, 5(2), 329-528. https://doi.org/10.1063/1.555536
    """
    _wl_range = (0.1, 10)  # um, Sellmeier validity (see docstring)

    __slots__ = ["_A1", "_B1", "_A2", "_B2"]

    def __init__(self):
        super().__init__()

        """ Constants of dispersion formula """
        self._A1 = 0.92549
        self._B1 = 0.07376
        self._A2 = 6.96747
        self._B2 = 32.79

    @property
    def symbols(self):
        return [wl, T]

    def n_expr(self, pol='o'):
        """ Sympy expression, dispersion formula """
        return sympy.sqrt(1 + self._A1 * wl**2 / (wl**2 - self._B1**2) + self._A2 * wl**2 / (wl**2 - self._B2**2))
