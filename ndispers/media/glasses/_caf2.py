import sympy
from ndispers._baseclass import Medium, wl, T
from ndispers.helper import vars2

class CaF2(Medium):
    """
    Ca F_2 (Calcium fluolide) crystal

    - Point group : Fm3m
    - Crystal system : cubic
    - Tranparency range : 0.18 to 8 µm (depends on material grade)
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
    ----
    Malitson, Irving H. "A redetermination of some optical properties of calcium fluoride." Applied Optics 2.11 (1963): 1103-1107.

    Example
    -------
    >>> caf2 = nd.media.glasses.CaF2()
    >>> caf2.n(0.532, 25) # refractive index at 0.532µm and at 25 degC
    array(1.43535685)
    
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
     
    @property
    def constants(self):
        print(vars2(self))
    
    def n_expr(self, pol='o'):
        """ Sympy expression, dispersion formula """
        return sympy.sqrt(1 + self._A1 * wl**2 / (wl**2 - self._B1**2) + self._A2 * wl**2 / (wl**2 - self._B2**2) + self._A3 * wl**2 / (wl**2 - self._B3**2)) + self._dndT * (T - 24)
    
    def n(self, wl_um, T_deg):
        """
        Refractive index as a function of wavelength

        input
        ------
        wl_um   :  float or array_like, wavelength in µm
        T_degC  :  float or array_like, temperature of crystal in degree C.
        
        return
        -------
        Refractive index, float
        """
        return super().n(wl_um, T_deg, pol='o')

    def dn_wl(self, wl_um, T_deg):
        return super().dn_wl(wl_um, T_deg, pol='o')
    
    def d2n_wl(self, wl_um, T_deg):
        return super().d2n_wl(wl_um, T_deg, pol='o')

    def d3n_wl(self, wl_um, T_deg):
        return super().d3n_wl(wl_um, T_deg, pol='o')
    
    def GD(self, wl_um, T_deg):
        """Group Delay [fs/mm]"""
        return super().GD(wl_um, T_deg, pol='o')
    
    def GV(self, wl_um, T_deg):
        """Group Velocity [µm/fs]"""
        return super().GV(wl_um, T_deg, pol='o')
    
    def ng(self, wl_um, T_deg):
        """Group index, c/Group velocity"""
        return super().ng(wl_um, T_deg, pol='o')
    
    def GVD(self, wl_um, T_deg):
        """Group Delay Dispersion [fs^2/mm]"""
        return super().GVD(wl_um, T_deg, pol='o')
    
    def TOD(self, wl_um, T_deg):
        """Third Order Dispersion [fs^3/mm]"""
        return super().TOD(wl_um, T_deg, pol='o')