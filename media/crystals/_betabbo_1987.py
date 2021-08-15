import sympy
from sympy.utilities import lambdify
from ndispers._baseclass import Medium, wl, phi, theta, pi

class BetaBBO(Medium):
    """
    beta-BBO (beta-Ba B_2 O_4) crystal

    - Point group : 3m
    - Crystal class : Trigonal
    - Dielectic principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis

    Dispersion formula of refractive index
    ---------------------------------------
    n(wl_um) = sqrt(A_i + B_i/(wl_um**2 - C_i) - D_i * wl_um**2)  for i = o, e
    
    Validity range
    ---------------
    0.22 - 1.06 um

    Ref
    ----
    Eimerl, David, et al. "Optical, mechanical, and thermal properties of barium borate." Journal of applied physics 62.5 (1987): 1968-1983.

    Usage
    ------
    # Create an instance of BetaBBO object.
    bbo = ndispers.media.crystals.BetaBBO()
    # Get a refractive index for e-ray as a function of wavelength (um) and theta angle.
    bbo.n(0.6, pol='e', theta_rad=0.2*pi, phi_rad='arb')
    # Derivative dispersion quantities are also easily obtained.
    bbo.GVD(0.6, pol='e', theta_rad=0.2*pi, phi_rad='arb')

    @author: Akihiko Shimura
    """
    # __slots__ = ["_A_o", "_B_o", "_C_o", "_D_o", 
                #  "_A_e", "_B_e", "_C_e", "_D_e"]
    def __init__(self):
        super().__init__(self)

        """ Constants of dispersion formula """
        # For ordinary ray
        self._A_o = 2.7405
        self._B_o = 0.0184
        self._C_o = 0.0179
        self._D_o = 0.0155
        # For extraordinary ray
        self._A_e = 2.3730
        self._B_e = 0.0128
        self._C_e = 0.0156
        self._D_e = 0.0044
    
    @property
    def property(self):
        msg = ["A_o = %.4f" % self._A_o]
        msg += ["B_o = %.4f" % self._B_o]
        msg += ["C_o = %.4f" % self._C_o]
        msg += ["D_o = %.4f" % self._D_o]
        msg += ["A_e = %.4f" % self._A_e]
        msg += ["B_e = %.4f" % self._B_e]
        msg += ["C_e = %.4f" % self._C_e]
        msg += ["D_e = %.4f" % self._D_e]
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
    
    def n(self, wl_um, pol, theta_rad):
        """
        Refractive index as a function of wavelength, theta and phi angles for each eigen polarization of light.

        input
        ------
        wl_um     :  float, wavelength in um, valid from 0.22 to 1.06 um
        pol       :  str, 'o' or 'e', polarization of light
        theta_rad :  float, 0 to pi radians

        return
        -------
        Refractive index, float
        """
        return lambdify([wl, theta], self.n_expr(pol=pol), 'numpy')(wl_um, theta_rad)

    def dn_wl(self, wl_um, pol, theta_rad, phi_rad=0.0):
        return super().dn_wl(wl_um, pol, theta_rad, phi_rad)
    
    def d2n_wl(self, wl_um, pol, theta_rad, phi_rad=0.0):
        return super().d2n_wl(wl_um, pol, theta_rad, phi_rad)

    def d3n_wl(self, wl_um, pol, theta_rad, phi_rad=0.0):
        return super().d3n_wl(wl_um, pol, theta_rad, phi_rad)

    def GD(self, wl_um, pol, theta_rad, phi_rad=0.0):
        """Group Delay [fs/mm]"""
        return super().GD(wl_um, pol, theta_rad, phi_rad)
    
    def GV(self, wl_um, pol, theta_rad, phi_rad=0.0):
        """Group Velocity [um/fs]"""
        return super().GV(wl_um, pol, theta_rad, phi_rad)
    
    def ng(self, wl_um, pol, theta_rad, phi_rad=0.0):
        """Group index, c/Group velocity"""
        return super().ng(wl_um, pol, theta_rad, phi_rad)

    def GVD(self, wl_um, pol, theta_rad, phi_rad=0.0):
        """Group Delay Dispersion [fs^2/mm]"""
        return super().GVD(wl_um, pol, theta_rad, phi_rad)

    def TOD(self, wl_um, pol, theta_rad, phi_rad=0.0):
        """Third Order Dispersion [fs^3/mm]"""
        return super().TOD(wl_um, pol, theta_rad, phi_rad)