import sympy
from ndispers._baseclass import Medium, wl, T
from ndispers.helper import vars2

class Ge(Medium):
    """
    Ge (crystalline germanium)

    - Point group : m3̄m  (Oh); diamond structure, centrosymmetric
    - Crystal system : cubic
    - Transparency range : 1.8 to 15 µm (Tropf et al., Table 20; vendors quote to 23 µm)

    Sellmeier equation
    ------------------
        n(wl)**2 = A + B * wl**2 / (wl**2 - C) + D * wl**2 / (wl**2 - E)
        (C, E in µm**2)

    Validity range
    ---------------
    2 to 12 µm

    Note
    ----
    The Sellmeier equation has no temperature term: T_degC is accepted for
    signature uniformity and ignored, and dndT returns 0. n = 4.025 at 4 µm. The source also gives the temperature dependence, which is not implemented here (Ge's dn/dT is about +4e-4 /K).

    Ref
    ---
    Sellmeier equation (as tabulated in):
      Tropf, W. J., Thomas, M. E., & Harris, T. J. (1995). Properties of crystals and glasses. In M. Bass (Ed.), Handbook of Optics, Vol. II, Chapter 33, Table 22. McGraw-Hill.
    Underlying measurements:
      Barnes, N. P., & Piltch, M. S. (1979). Temperature-dependent Sellmeier coefficients and nonlinear optics average power limit for germanium. JOSA, 69(1), 178-180. https://doi.org/10.1364/josa.69.000178
    """
    __slots__ = ["_A", "_B", "_C", "_D", "_E"]

    def __init__(self):
        super().__init__()

        """ Constants of dispersion formula """
        self._A = 9.28156
        self._B = 6.7288
        self._C = 0.44105
        self._D = 0.21307
        self._E = 3870.1

    @property
    def symbols(self):
        return [wl, T]

    def n_expr(self, pol='o'):
        """ Sympy expression, dispersion formula """
        return sympy.sqrt(self._A + self._B * wl**2 / (wl**2 - self._C) + self._D * wl**2 / (wl**2 - self._E))
