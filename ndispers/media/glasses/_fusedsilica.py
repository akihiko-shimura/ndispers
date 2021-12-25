import sympy
from sympy.utilities import lambdify
from ndispers._baseclass import Medium, wl, T
from ndispers.helper import vars2

class FusedSilica(Medium):
    """
    Fused Silica glass
    
    Dispersion formula for refractive index
    ---------------------------------------
    n(wl_um) = sqrt(1 + B1 * wl_um**2/(wl_um**2 - C1) + B2 * wl_um**2/(wl_um**2 - C2) + B3 * wl_um**2/(wl_um**2 - C3))

    Thermo-optic coefficient
    ------------------------
    dn/dT = 11.3e-6 /K
    
    Validity range
    ---------------
    0.21 to 3.71 um

    Ref
    ----
    - W. S. Rodney and R. J. Spindler, "Index of Refraction of Fused-quartz Glass for Ultraviolet, Visible, and Infrared Wavelengths" J. Res. Nat. Bur. Stand. 53:185–189 (1954)
    - I. H. Malitson, "Interspecimen Comparison of the Refractive Index of Fused Silica" J. Opt. Soc. Am. 55 :1205-1209 (1965)
    - Rocha, A. C. P., et al. "Measurements of refractive indices and thermo-optical coefficients using a white-light Michelson interferometer." Applied optics 55.24 (2016): 6639-6643.

    Examples
    ---------
    
    
    """
    __slots__ = ["_B1", "_C1", "_B2", "_C2", "_B3", "_C3", "_dndT"]

    def __init__(self):
        super().__init__()

        """ Constants of dispersion formula """
        # For ordinary ray
        self._B1 = 0.6961663
        self._C1 = 0.0684043**2
        self._B2 = 0.4079426
        self._C2 = 0.1162414**2
        self._B3 = 0.8974794
        self._C3 = 9.896161**2
        self._dndT = 11.3e-6 #/K
    
    @property
    def symbols(self):
        return [wl, T]
     
    @property
    def constants(self):
        print(vars2(self))
    
    def n_expr(self, pol='o'):
        """ Sympy expression, dispersion formula """
        return sympy.sqrt(1 + self._B1 * wl**2 / (wl**2 - self._C1) + self._B2 * wl**2 / (wl**2 - self._C2) + self._B3 * wl**2 / (wl**2 - self._C3)) + self._dndT * (T - 24)
    
    def n(self, wl_um, T_degC):
        """
        Refractive index as a function of wavelength

        input
        ------
        wl_um   :  float or array_like, wavelength in um
        T_degC  :  float or array_like, temperature of crystal in degree C.

        return
        -------
        Refractive index, float
        """
        return super().n(wl_um, T_degC, pol='o')

    def dn_wl(self, wl_um, T_degC):
        return super().dn_wl(wl_um, T_degC, pol='o')
    
    def d2n_wl(self, wl_um, T_degC):
        return super().d2n_wl(wl_um, T_degC, pol='o')

    def d3n_wl(self, wl_um, T_degC):
        return super().d3n_wl(wl_um, T_degC, pol='o')
    
    def GD(self, wl_um, T_degC):
        """Group Delay [fs/mm]"""
        return super().GD(wl_um, T_degC, pol='o')
    
    def GV(self, wl_um, T_degC):
        """Group Velocity [um/fs]"""
        return super().GV(wl_um, T_degC, pol='o')
    
    def ng(self, wl_um, T_degC):
        """Group index, c/Group velocity"""
        return super().ng(wl_um, T_degC, pol='o')
    
    def GVD(self, wl_um, T_degC):
        """Group Delay Dispersion [fs^2/mm]"""
        return super().GVD(wl_um, T_degC, pol='o')
    
    def TOD(self, wl_um, T_degC):
        """Third Order Dispersion [fs^3/mm]"""
        return super().TOD(wl_um, T_degC, pol='o')