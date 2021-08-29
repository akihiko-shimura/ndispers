import sympy
from ndispers._baseclass import Medium, wl, phi, theta, pi
import numpy

class LB4(Medium):
    """
    LB4 or LTB (Li_2 B_4 O_7, lithium tetraborate) crystal

    - Point group : 4mm
    - Crystal ststem : Tetragonal
    - Dielectic principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - - Tranparency range : 0.16 - 3.6 um

    Dispersion formula of refractive index
    ---------------------------------------
    n(wl_um) = sqrt(A_i + B_i/(wl**2 - C_i) - D_i * wl**2)  for i = o, e
    
    Validity range
    ---------------
    0.18 - 2.3 um
    T = 25 degC
    P = 101325 Pa

    Ref
    ----
    Sugawara, Tamotsu, Ryuichi Komatsu, and Satoshi Uda. "Linear and nonlinear optical properties of lithium tetraborate." Solid state communications 107.5 (1998): 233-237.

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
    __slots__ = ["_A_o", "_B_o", "_C_o", "_D_o", 
                 "_A_e", "_B_e", "_C_e", "_D_e"]

    def __init__(self):
        super().__init__()

        """ Constants of dispersion formula """
        # For ordinary ray
        self._A_o = 2.564310
        self._B_o = 0.012337
        self._C_o = 0.114467**2
        self._D_o = 0.019075
        # For extraordinary ray
        self._A_e = 2.386510
        self._B_e = 0.010664
        self._C_e = 0.113483**2
        self._D_e = 0.012813
    
    @property
    def property(self):
        msg  = ["A_o = %g" % self._A_o]
        msg += ["B_o = %g" % self._B_o]
        msg += ["C_o = %g" % self._C_o]
        msg += ["D_o = %g" % self._D_o]
        msg += ["A_e = %g" % self._A_e]
        msg += ["B_e = %g" % self._B_e]
        msg += ["C_e = %g" % self._C_e]
        msg += ["D_e = %g" % self._D_e]
        print("\n".join(msg))
    
    def n_o_expr(self):
        """ Sympy expression, dispersion formula for o-ray """
        return sympy.sqrt(self._A_o + self._B_o / (wl**2 - self._C_o) - self._D_o * wl**2)
    
    def n_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-ray """
        return sympy.sqrt(self._A_e + self._B_e / (wl**2 - self._C_e) - self._D_e * wl**2)

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
