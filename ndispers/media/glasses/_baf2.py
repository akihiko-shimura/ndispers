from ndispers._sym import sympy
from ndispers._baseclass import Medium, wl, T
from ndispers.helper import vars2

class BaF2(Medium):
    """
    BaF₂ (barium fluoride) crystal

    - Point group : m3̄m  (Oh); space group Fm3̄m
    - Crystal system : cubic
    - Transparency range : about 0.15 to 12 µm (the mid-infrared window among the fluorides; the source's summary-table page was not legible, so this is the usual handbook figure)

    Sellmeier equation
    ------------------
        n(wl)**2 = 1 + A1 * wl**2 / (wl**2 - B1**2) + A2 * wl**2 / (wl**2 - B2**2) + A3 * wl**2 / (wl**2 - B3**2)

    Validity range
    ---------------
    0.27 to 10.3 µm

    Note
    ----
    The Sellmeier equation has no temperature term: T_degC is accepted for
    signature uniformity and ignored, and dndT returns 0. 

    Ref
    ---
    Sellmeier equation (as tabulated in):
      Tropf, W. J., Thomas, M. E., & Harris, T. J. (1995). Properties of crystals and glasses. In M. Bass (Ed.), Handbook of Optics, Vol. II, Chapter 33, Table 22. McGraw-Hill.
    Underlying measurements:
      Malitson, I. H. (1964). Refractive properties of barium fluoride. JOSA, 54(5), 628-632. https://doi.org/10.1364/josa.54.000628
    """
    _wl_range = (0.27, 10.3)  # um, Sellmeier validity (see docstring)

    __slots__ = ["_A1", "_B1", "_A2", "_B2", "_A3", "_B3"]

    def __init__(self):
        super().__init__()

        """ Constants of dispersion formula """
        self._A1 = 0.643356
        self._B1 = 0.057789
        self._A2 = 0.506762
        self._B2 = 0.10968
        self._A3 = 3.8261
        self._B3 = 46.3864

    @property
    def symbols(self):
        return [wl, T]

    def n_expr(self, pol='o'):
        """ Sympy expression, dispersion formula """
        return sympy.sqrt(1 + self._A1 * wl**2 / (wl**2 - self._B1**2) + self._A2 * wl**2 / (wl**2 - self._B2**2) + self._A3 * wl**2 / (wl**2 - self._B3**2))
