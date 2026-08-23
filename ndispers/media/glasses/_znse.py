from ndispers._sym import sympy
from ndispers._baseclass import Medium, wl, T
from ndispers.helper import vars2

class ZnSe(Medium):
    """
    ZnSe (zinc selenide) crystal, CVD polycrystalline grade

    - Point group : 4̄3m  (Td); cubic (zincblende) - polycrystalline CVD material is optically isotropic
    - Crystal system : cubic
    - Transparency range : 0.55 to 18 µm (the range of the dispersion fit); the standard mid-infrared window and lens material

    Sellmeier equation
    ------------------
        n(wl)**2 = 1 + A1 * wl**2 / (wl**2 - B1**2) + A2 * wl**2 / (wl**2 - B2**2) + A3 * wl**2 / (wl**2 - B3**2)

    Validity range
    ---------------
    0.55 to 18 µm

    Note
    ----
    The Sellmeier equation has no temperature term: T_degC is accepted for
    signature uniformity and ignored, and dndT returns 0. The fit agrees with the Connolly et al. 1979 (Raytran) set on refractiveindex.info to 4e-4 from 0.63 to 10.6 µm. Zincblende ZnSe is non-centrosymmetric (d14 ~ 30 pm/V) but, being cubic, cannot phase-match birefringently; it is listed here as a window material.

    Ref
    ---
    Sellmeier equation (as tabulated in):
      Tropf, W. J., Thomas, M. E., & Harris, T. J. (1995). Properties of crystals and glasses. In M. Bass (Ed.), Handbook of Optics, Vol. II, Chapter 33, Table 22. McGraw-Hill.
    Underlying measurements:
      Feldman, A., Horowitz, D., Waxler, R. M., & Dodge, M. J. (1978). Optical materials characterization: final technical report, February 1, 1978 - September 30, 1978. NBS Technical Note 993, February 1979 (not registered with Crossref).
    """
    __slots__ = ["_A1", "_B1", "_A2", "_B2", "_A3", "_B3"]

    def __init__(self):
        super().__init__()

        """ Constants of dispersion formula """
        self._A1 = 4.2980149
        self._B1 = 0.1920630
        self._A2 = 0.62776557
        self._B2 = 0.37878260
        self._A3 = 2.8955633
        self._B3 = 46.994595

    @property
    def symbols(self):
        return [wl, T]

    def n_expr(self, pol='o'):
        """ Sympy expression, dispersion formula """
        return sympy.sqrt(1 + self._A1 * wl**2 / (wl**2 - self._B1**2) + self._A2 * wl**2 / (wl**2 - self._B2**2) + self._A3 * wl**2 / (wl**2 - self._B3**2))
