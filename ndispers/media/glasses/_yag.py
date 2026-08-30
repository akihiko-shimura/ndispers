from ndispers._sym import sympy
from ndispers._baseclass import Medium, wl, T
from ndispers.helper import vars2

class YAG(Medium):
    """
    YAG (Y₃Al₅O₁₂, yttrium aluminium garnet) crystal, undoped

    - Point group : m3̄m  (Oh); space group Ia3̄d
    - Crystal system : cubic
    - Transparency range : about 0.21 to 5.5 µm

    Sellmeier equation
    ------------------
        n(wl)**2 = 1 + A * wl**2 / (wl**2 - B) + C * wl**2 / (wl**2 - D)
        (B, D in µm**2)

    Validity range
    ---------------
    0.4 to 5.0 µm

    Note
    ----
    The Sellmeier equation has no temperature term: T_degC is accepted for
    signature uniformity and ignored, and dndT returns 0. Coefficients as
    tabulated by refractiveindex.info (main/Y3Al5O12/nk/Zelmon; accessed
    2026-08-23) from the source; n = 1.8328 at 587.6 nm and 1.8147 at 1064 nm,
    within 3e-4 of the Hrabovský et al. 2021 fit on the same site.

    Ref
    ---
    Sellmeier equation:
      Zelmon, D. E., Small, D. L., & Page, R. (1998). Refractive-index measurements of undoped yttrium aluminum garnet from 0.4 to 5.0 µm. Applied Optics, 37(21), 4933-4935. https://doi.org/10.1364/ao.37.004933
    """
    _wl_range = (0.4, 5.0)  # um, Sellmeier validity (see docstring)

    __slots__ = ["_A", "_B", "_C", "_D"]

    def __init__(self):
        super().__init__()

        """ Constants of dispersion formula """
        self._A = 2.28200
        self._B = 0.01185
        self._C = 3.27644
        self._D = 282.734

    @property
    def symbols(self):
        return [wl, T]

    def n_expr(self, pol='o'):
        """ Sympy expression, dispersion formula """
        return sympy.sqrt(1 + self._A * wl**2 / (wl**2 - self._B) + self._C * wl**2 / (wl**2 - self._D))
