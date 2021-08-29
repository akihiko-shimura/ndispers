import sympy
from ndispers._baseclass import Medium, wl, phi, theta, pi
import numpy

class BetaBBO(Medium):
    """
    beta-BBO (beta-Ba B_2 O_4) crystal

    - Point group : 3m
    - Crystal ststem : Trigonal
    - Dielectic principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Tranparency range : 1.9 - 2.6 um

    Dispersion formula of refractive index
    ---------------------------------------
    n(wl_um) = sqrt(1 + B1_i*wl**2/(wl**2 - C1_i) + B2_i*wl**2/(wl**2 - C2_i) + B3_i*wl**2/(wl**2 - C3_i))  for i = o, e
    
    Validity range
    ---------------
    0.188 - 5.2 um

    Ref
    ----
    Tamošauskas, Gintaras, et al. "Transmittance and phase matching of BBO crystal in the 3−5 μm range and its application for the characterization of mid-infrared laser pulses." Optical Materials Express 8.6 (2018): 1410-1418.

    Usage
    ------
    # Create an instance of BetaBBO object.
    bbo = ndispers.media.crystals.BetaBBO()
    # Get a refractive index for e-ray as a function of wavelength (um), thetaand phi angles (rad).
    bbo.n(0.6, 0, 0, pol='o') # here, n does not depend on theta and phi for o-ray
    # For e-ray, n depends on theta,
    bbo.n(0.6, 0.5*pi, 0, pol='e') # along z-axis
    bbo.n(0.6, 0.23*pi, 0, pol='e') # for theta = 0.23*pi rad
    bbo.n(0.6, 0*pi, 0, pol='e') # for theta = 0 rad, this corresponds to the case of o-ray
    # Derivative dispersion quantities are also easily obtained.
    bbo.GVD(0.6, 0.23*pi, 0, pol='e')

    @author: Akihiko Shimura
    """
    __slots__ = ["_B1_o", "_C1_o", "_B2_o", "_C2_o", "_B3_o", "_C3_o",
                 "_B1_e", "_C1_e", "_B2_e", "_C2_e", "_B3_e", "_C3_e"]

    def __init__(self):
        super().__init__()

        """ Constants of dispersion formula """
        # For ordinary ray
        self._B1_o = 0.90291
        self._C1_o = 0.003926
        self._B2_o = 0.83155
        self._C2_o = 0.018786
        self._B3_o = 0.76536
        self._C3_o = 60.01
        # For extraordinary ray
        self._B1_e = 1.151075
        self._C1_e = 0.007142
        self._B2_e = 0.21803
        self._C2_e = 0.02259
        self._B3_e = 0.656
        self._C3_e = 263
    
    @property
    def property(self):
        msg = ["B1_o = %g" % self._B1_o]
        msg += ["C1_o = %g" % self._C1_o]
        msg += ["B2_o = %g" % self._B2_o]
        msg += ["C2_o = %g" % self._C2_o]
        msg += ["B3_o = %g" % self._B3_o]
        msg += ["C3_o = %g" % self._C3_o]
        msg += ["B1_e = %g" % self._B1_e]
        msg += ["C1_e = %g" % self._C1_e]
        msg += ["B2_e = %g" % self._B2_e]
        msg += ["C2_e = %g" % self._C2_e]
        msg += ["B3_e = %g" % self._B3_e]
        msg += ["C3_e = %g" % self._C3_e]
        print("\n".join(msg))
    
    def n_o_expr(self):
        """ Sympy expression, dispersion formula for o-ray """
        return sympy.sqrt(1.0 + self._B1_o * wl**2/ (wl**2 - self._C1_o) + self._B2_o * wl**2/ (wl**2 - self._C2_o) + self._B3_o * wl**2/ (wl**2 - self._C3_o))
    
    def n_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-ray """
        return sympy.sqrt(1.0 + self._B1_e * wl**2/ (wl**2 - self._C1_e) + self._B2_e * wl**2/ (wl**2 - self._C2_e) + self._B3_e * wl**2/ (wl**2 - self._C3_e))

    def n_expr(self, pol):
        """"
        Sympy expression, 
        dispersion formula of a general ray with an angle theta to optic axis. If theta = 0, this expression reduces to 'no_expre'.

        n(theta) = n_e / sqrt( sin(theta)**2 + (n_e/n_o)**2 * cos(theta)**2 )
        """
        if pol == 'o':
            return self.n_o_expr()
        elif pol == 'e':
            return self.n_e_expr() / sympy.sqrt( sympy.sin(theta)**2 + (self.n_e_expr()/self.n_o_expr())**2 * sympy.cos(theta)**2 )
        else:
            raise ValueError("pol = '%s' must be 'o' or 'e'" % pol)

    
    def n(self, wl_um, theta_rad, phi_rad, pol='o'):
        """
        Refractive index as a function of wavelength, theta and phi angles for each eigen polarization of light.

        input
        ------
        wl_um     :  float, wavelength in um
        pol       :  str, 'o' or 'e', polarization of light
        theta_rad :  float, 0 to pi radians
        phi_rad   :  float, option, this is arbitray and the result does not depend on this value.
s
        return
        -------
        Refractive index, float
        """
        return super().n(wl_um, theta_rad, phi_rad, pol=pol)

    def dn_wl(self, wl_um, theta_rad, phi_rad, pol='o'):
        return super().dn_wl(wl_um, theta_rad, phi_rad, pol=pol)
    
    def d2n_wl(self, wl_um, theta_rad, phi_rad, pol='o'):
        return super().d2n_wl(wl_um, theta_rad, phi_rad, pol=pol)

    def d3n_wl(self, wl_um, theta_rad, phi_rad, pol='o'):
        return super().d3n_wl(wl_um, theta_rad, phi_rad, pol=pol)
    
    def GD(self, wl_um, theta_rad, phi_rad, pol='o'):
        """Group Delay [fs/mm]"""
        return super().GD(wl_um, theta_rad, phi_rad, pol=pol)
    
    def GV(self, wl_um, theta_rad, phi_rad, pol='o'):
        """Group Velocity [um/fs]"""
        return super().GV(wl_um, theta_rad, phi_rad, pol=pol)
    
    def ng(self, wl_um, theta_rad, phi_rad, pol='o'):
        """Group index, c/Group velocity"""
        return super().ng(wl_um, theta_rad, phi_rad, pol=pol)
    
    def GVD(self, wl_um, theta_rad, phi_rad, pol='o'):
        """Group Delay Dispersion [fs^2/mm]"""
        return super().GVD(wl_um, theta_rad, phi_rad, pol=pol)
    
    def TOD(self, wl_um, theta_rad, phi_rad, pol='o'):
        """Third Order Dispersion [fs^3/mm]"""
        return super().TOD(wl_um, theta_rad, phi_rad, pol=pol)