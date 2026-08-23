import sympy
from ndispers._baseclass import Medium, wl, T
from ndispers.helper import vars2

class ZnS(Medium):
    """
    ZnS (zinc sulfide) crystal, cubic (CVD / Cleartran grade)

    - Point group : 4̄3m  (Td); cubic (zincblende) - polycrystalline CVD material is optically isotropic
    - Crystal system : cubic
    - Transparency range : 0.4 to 13 µm (the range of the dispersion fit)

    Sellmeier equation
    ------------------
        n(wl)**2 = A + B / (wl**2 - C**2) + D / (wl**2 - E**2)

    Validity range
    ---------------
    0.405 to 13 µm, at 20 degC

    Note
    ----
    The Sellmeier equation has no temperature term: T_degC is accepted for
    signature uniformity and ignored, and dndT returns 0. Coefficients as tabulated by refractiveindex.info (main/ZnS/nk/Debenham; accessed 2026-08-23) from the source; n = 2.288 at 1.064 µm.

    Ref
    ---
    Sellmeier equation:
      Debenham, M. (1984). Refractive indices of zinc sulfide in the 0.405-13-µm wavelength range. Applied Optics, 23(14), 2238-2239. https://doi.org/10.1364/ao.23.002238
    """
    __slots__ = ["_A", "_B", "_C", "_D", "_E"]

    def __init__(self):
        super().__init__()

        """ Constants of dispersion formula """
        self._A = 8.393
        self._B = 0.14383
        self._C = 0.2421
        self._D = 4430.99
        self._E = 36.71

    @property
    def symbols(self):
        return [wl, T]

    def n_expr(self, pol='o'):
        """ Sympy expression, dispersion formula """
        return sympy.sqrt(self._A + self._B / (wl**2 - self._C**2) + self._D / (wl**2 - self._E**2))
