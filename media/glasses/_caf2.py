import sympy
from sympy.utilities import lambdify
from ndispers._baseclass import Medium, wl, phi, theta, pi
from helper import vars2

class CaF2(Medium):
    """
    Fused Silica glass
    
    Dispersion formula of refractive index
    ---------------------------------------
    n(wl_um) = 
    
    Validity range
    ---------------
    

    Ref
    ----
    日本結晶光学株式会社カタログ，光学結晶CaF2

    Usage
    ------
    
    @author: Akihiko Shimura
    """
    __slots__ = ["_A1", "_B1", "_A2", "_B2", "_A3", "_B3"]

    def __init__(self):
        super().__init__()

        """ Constants of dispersion formula """
        # For ordinary ray
        self._A1 = 6.254288046e-1
        self._B1 = 2.813183822e-3
        self._A2 = 4.132684951e-1
        self._B2 = 1.066206606e-2
        self._A3 = 3.409193892
        self._B3 = 1.065596428e-3
    
    @property
    def symbols(self):
        return [wl, theta, phi]
     
    @property
    def constants(self):
        print(vars2(self))
    
    def n_expr(self, pol='o'):
        """ Sympy expression, dispersion formula """
        return sympy.sqrt(1 + self._A1 * wl**2 / (wl**2 - self._B1) + self._A2 * wl**2 / (wl**2 - self._B2) + self._A3 * wl**2 / (wl**2 - self._B3))
    
    def n(self, wl_um, pol='o'):
        """
        Refractive index as a function of wavelength, theta and phi angles for each eigen polarization of light.

        input
        ------
        wl_um     :  float, wavelength in um

        return
        -------
        Refractive index, float
        """
        return super().n(wl_um, 0, 0, pol=pol)

    def dn_wl(self, wl_um, pol='o'):
        return super().dn_wl(wl_um, 0, 0, pol=pol)
    
    def d2n_wl(self, wl_um, pol='o'):
        return super().d2n_wl(wl_um, 0, 0, pol=pol)

    def d3n_wl(self, wl_um, pol='o'):
        return super().d3n_wl(wl_um, 0, 0, pol=pol)
    
    def GD(self, wl_um, pol='o'):
        """Group Delay [fs/mm]"""
        return super().GD(wl_um, 0, 0, pol=pol)
    
    def GV(self, wl_um, pol='o'):
        """Group Velocity [um/fs]"""
        return super().GV(wl_um, 0, 0, pol=pol)
    
    def ng(self, wl_um, pol='o'):
        """Group index, c/Group velocity"""
        return super().ng(wl_um, 0, 0, pol=pol)
    
    def GVD(self, wl_um, pol='o'):
        """Group Delay Dispersion [fs^2/mm]"""
        return super().GVD(wl_um, 0, 0, pol=pol)
    
    def TOD(self, wl_um, pol='o'):
        """Third Order Dispersion [fs^3/mm]"""
        return super().TOD(wl_um, 0, 0, pol=pol)