import sympy
from sympy.utilities import lambdify
from ndispers._baseclass import Medium, wl, phi, theta, pi

class FusedSilica(Medium):
    """
    Fused Silica glass
    
    Dispersion formula of refractive index
    ---------------------------------------
    n(wl_um) = sqrt(1 + B1 * wl_um**2/(wl_um**2 - C1) + B2 * wl_um**2/(wl_um**2 - C2) + B3 * wl_um**2/(wl_um**2 - C3))
    
    Validity range
    ---------------
    0.21 – 3.71 um
    T = 20 degC

    Ref
    ----
    - W. S. Rodney and R. J. Spindler, "Index of Refraction of Fused-quartz Glass for Ultraviolet, Visible, and Infrared Wavelengths" J. Res. Nat. Bur. Stand. 53:185–189 (1954)
    - I. H. Malitson, "Interspecimen Comparison of the Refractive Index of Fused Silica" J. Opt. Soc. Am. 55 :1205-1209 (1965)

    Usage
    ------
    
    @author: Akihiko Shimura
    """
    __slots__ = ["_B1", "_C1", "_B2", "_C2", "_B3", "_C3"]

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
     
    @property
    def property(self):
        msg  = ["B1 = %g" % self._B1]
        msg += ["C1 = %g" % self._C1]
        msg += ["B2 = %g" % self._B2]
        msg += ["C2 = %g" % self._C2]
        msg += ["B3 = %g" % self._B3]
        msg += ["C3 = %g" % self._C3]
        print("\n".join(msg))
    
    def n_expr(self, pol='o'):
        """ Sympy expression, dispersion formula """
        return sympy.sqrt(1 + self._B1 * wl**2 / (wl**2 - self._C1) + self._B2 * wl**2 / (wl**2 - self._C2) + self._B3 * wl**2 / (wl**2 - self._C3))
    
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