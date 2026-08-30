from ndispers._sym import sympy
from ndispers._baseclass import Medium, wl, T
from ndispers.helper import vars2

class Si(Medium):
    """
    Si (crystalline silicon)

    - Point group : m3̄m  (Oh); diamond structure, centrosymmetric
    - Crystal system : cubic
    - Transparency range : 1.1 to 6.5 µm (Tropf et al., Table 20; lattice bands beyond, with a further window past 8 µm in vendor data)

    Sellmeier equation
    ------------------
        n(wl)**2 = 1 + A1 * wl**2 / (wl**2 - B1**2) + A2 * wl**2 / (wl**2 - B2**2) + A3 * wl**2 / (wl**2 - B3**2)

    Validity range
    ---------------
    1.36 to 11 µm

    Note
    ----
    The Sellmeier equation has no temperature term: T_degC is accepted for
    signature uniformity and ignored, and dndT returns 0. n = 3.4527 at 2 µm. Silicon's dn/dT is large (about +1.6e-4 /K) and is not included.

    Ref
    ---
    Sellmeier equation (as tabulated in):
      Tropf, W. J., Thomas, M. E., & Harris, T. J. (1995). Properties of crystals and glasses. In M. Bass (Ed.), Handbook of Optics, Vol. II, Chapter 33, Table 22. McGraw-Hill.
    Underlying fit:
      Tatian, B. (1984). Fitting refractive-index data with the Sellmeier dispersion formula. Applied Optics, 23(24), 4477-4485. https://doi.org/10.1364/ao.23.004477
    """
    _wl_range = (1.36, 11)  # um, Sellmeier validity (see docstring)

    __slots__ = ["_A1", "_B1", "_A2", "_B2", "_A3", "_B3"]

    def __init__(self):
        super().__init__()

        """ Constants of dispersion formula """
        self._A1 = 10.6684293
        self._B1 = 0.301516485
        self._A2 = 0.003043475
        self._B2 = 1.13475115
        self._A3 = 1.54133408
        self._B3 = 1104.0

    @property
    def symbols(self):
        return [wl, T]

    def n_expr(self, pol='o'):
        """ Sympy expression, dispersion formula """
        return sympy.sqrt(1 + self._A1 * wl**2 / (wl**2 - self._B1**2) + self._A2 * wl**2 / (wl**2 - self._B2**2) + self._A3 * wl**2 / (wl**2 - self._B3**2))
