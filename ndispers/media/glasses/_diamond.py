import sympy
from ndispers._baseclass import Medium, wl, T
from ndispers.helper import vars2

class Diamond(Medium):
    """
    Diamond (crystalline carbon)

    - Point group : m3̄m  (Oh); diamond structure, centrosymmetric
    - Crystal system : cubic
    - Transparency range : 0.225 µm into the far infrared (two-phonon absorption between about 2.5 and 6.5 µm)

    Sellmeier equation
    ------------------
        n(wl)**2 = 1 + A1 * wl**2 / (wl**2 - B1**2) + A2 * wl**2 / (wl**2 - B2**2)

    Validity range
    ---------------
    0.226 to 0.644 µm measured (Peter 1923; Tropf tabulates the fit as 0.225 µm upward); the two-UV-pole fit has no infrared term and is commonly extrapolated through the infrared, where diamond's dispersion is small

    Note
    ----
    The Sellmeier equation has no temperature term: T_degC is accepted for
    signature uniformity and ignored, and dndT returns 0. n = 2.4175 at 587.6 nm.

    Ref
    ---
    Sellmeier equation (as tabulated in):
      Tropf, W. J., Thomas, M. E., & Harris, T. J. (1995). Properties of crystals and glasses. In M. Bass (Ed.), Handbook of Optics, Vol. II, Chapter 33, Table 22. McGraw-Hill.
    Underlying measurements:
      Peter, F. (1923). Über Brechungsindizes und Absorptionskonstanten des Diamanten zwischen 644 und 226 mµ. Zeitschrift für Physik, 15(1), 358-368. https://doi.org/10.1007/bf01330487
    """
    __slots__ = ["_A1", "_B1", "_A2", "_B2"]

    def __init__(self):
        super().__init__()

        """ Constants of dispersion formula """
        self._A1 = 4.3356
        self._B1 = 0.106
        self._A2 = 0.3306
        self._B2 = 0.175

    @property
    def symbols(self):
        return [wl, T]

    def n_expr(self, pol='o'):
        """ Sympy expression, dispersion formula """
        return sympy.sqrt(1 + self._A1 * wl**2 / (wl**2 - self._B1**2) + self._A2 * wl**2 / (wl**2 - self._B2**2))
