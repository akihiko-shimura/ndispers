import sympy
from ndispers._baseclass import Medium, wl, T
from ndispers.helper import vars2

class NBK7(Medium):
    """
    N-BK7 (SCHOTT borosilicate crown glass)

    - Amorphous, optically isotropic
    - Transparency range : about 0.33 to 2.5 µm (internal transmittance > 0.9 per 10 mm)

    Sellmeier equation
    ------------------
        n(wl)**2 = 1 + B1 * wl**2 / (wl**2 - C1) + B2 * wl**2 / (wl**2 - C2) + B3 * wl**2 / (wl**2 - C3)
        (C in µm**2; the three-term form of the SCHOTT catalog, n_d = 1.5168)

    Validity range
    ---------------
    0.3 to 2.5 µm, at 20 degC (catalog data, relative to air)

    Note
    ----
    The Sellmeier equation has no temperature term: T_degC is accepted for
    signature uniformity and ignored, and dndT returns 0. SCHOTT's dn/dT
    dispersion formula (TIE-19) is not implemented.

    Ref
    ---
    Sellmeier equation:
      SCHOTT N-BK7 (vendor data, Zemax catalog 2017-01-20b, accessed 2026-08-23 via refractiveindex.info specs/schott/optical/N-BK7). https://www.schott.com/en-us/products/optical-glass-p1000267/downloads/
    """
    __slots__ = ["_B1", "_C1", "_B2", "_C2", "_B3", "_C3"]

    def __init__(self):
        super().__init__()

        """ Constants of dispersion formula """
        self._B1 = 1.03961212
        self._C1 = 0.00600069867
        self._B2 = 0.231792344
        self._C2 = 0.0200179144
        self._B3 = 1.01046945
        self._C3 = 103.560653

    @property
    def symbols(self):
        return [wl, T]

    def n_expr(self, pol='o'):
        """ Sympy expression, dispersion formula """
        return sympy.sqrt(1 + self._B1 * wl**2 / (wl**2 - self._C1) + self._B2 * wl**2 / (wl**2 - self._C2) + self._B3 * wl**2 / (wl**2 - self._C3))
