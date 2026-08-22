import sympy
from ndispers._baseclass import Medium, wl, T
from ndispers.helper import vars2

class CaF2(Medium):
    """
    CaF₂ (calcium fluoride) crystal

    - Point group : m3̄m  (Oh); space group Fm3̄m
    - Crystal system : cubic
    - Transparency range : 0.18 to 8 µm (depends on material grade)
    - Transmission Range : 0.13 to 10 µm (depends on material grade)
    
    Sellmeier equation
    ---------------------------------------
        n(wl) = sqrt(1 + A1 * wl**2 / (wl**2 - B1**2) + A2 * wl**2 / (wl**2 - B2**2) + A3 * wl**2 / (wl**2 - B3**2))

    Thermo-optic coefficient
    ------------------------
        dn/dT = -10.6e-6 /K around T=24 degC
    
    Validity range
    ---------------
    0.23 to 9.7 µm

    Ref
    ---
    Sellmeier equation:
      Malitson, I. H. (1963). A redetermination of some optical properties of calcium fluoride. Applied Optics, 2(11), 1103-1107. https://doi.org/10.1364/ao.2.001103

    Thermo-optic coefficient:
      Rocha, A. C. P., Silva, J. R., Lima, S. M., et al. (2016). Measurements of refractive indices and thermo-optical coefficients using a white-light Michelson interferometer. Applied Optics, 55(24), 6639-6643. https://doi.org/10.1364/ao.55.006639
    """
    __slots__ = ["_A1", "_B1", "_A2", "_B2", "_A3", "_B3", "_dndT"]

    def __init__(self):
        super().__init__()
        
        """ Constants of dispersion formula """
        # For ordinary ray
        self._A1 = 0.5675888
        self._B1 = 0.050263605
        self._A2 = 0.4710914
        self._B2 = 0.1003909
        self._A3 = 3.8484723
        self._B3 = 34.649040
        self._dndT = -10.6e-6 #1/K
    
    @property
    def symbols(self):
        return [wl, T]
     
    
    def n_expr(self, pol='o'):
        """ Sympy expression, dispersion formula """
        return sympy.sqrt(1 + self._A1 * wl**2 / (wl**2 - self._B1**2) + self._A2 * wl**2 / (wl**2 - self._B2**2) + self._A3 * wl**2 / (wl**2 - self._B3**2)) + self._dndT * (T - 24)
    

    

    
    
    
    
    