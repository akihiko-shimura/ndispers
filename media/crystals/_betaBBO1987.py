import sympy
from sympy.utilities import lambdify
from ndispers._baseclass import Medium, wl, phi, theta, T, pi

class BetaBBO(Medium):
    """
    beta-BBO (beta-Ba B_2 O_4) crystal

    - Point group : 3m
    - Crystal system : Trigonal
    - Dielectic principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Tranparency range : 1.9 - 2.6 um

    Dispersion formula of refractive index
    ---------------------------------------
    n(wl_um) = sqrt(A_i + B_i/(wl_um**2 - C_i) - D_i * wl_um**2)  for i = o, e
    
    Validity range
    ---------------
    0.22 - 1.06 um


    Ref
    ----
    Eimerl, David, et al. "Optical, mechanical, and thermal properties of barium borate." Journal of applied physics 62.5 (1987): 1968-1983.
    Nikogosyan, D. N. "Beta barium borate (BBO)." Applied Physics A 52.6 (1991): 359-368.

    Usage
    ------
    >>> bbo = ndispers.media.crystals.BetaBBO1987()
    >>> bbo.n(0.6, 0, pol='o') # for o-ray, n does not depend on theta.
    >>> bbo.n(0.6, 0.5*pi, pol='e') # along z-axis, it is pure e-ray.
    >>> bbo.n(0.6, 0.23*pi, pol='e')
    >>> bbo.n(0.6, 0*pi, pol='e') # for theta = 0 rad, it corresponds to o-ray.
    >>> bbo.GVD(0.6, 0.23*pi, pol='e')

    @author: Akihiko Shimura
    """
    __slots__ = ["_BetaBBO__plane", "_BetaBBO__theta_rad", "_BetaBBO__phi_rad",
                 "_A_o", "_B_o", "_C_o", "_D_o",
                 "_A_e", "_B_e", "_C_e", "_D_e",
                 "_dndT_o", "_dndT_e"]
                 
    def __init__(self):
        super().__init__()
        self._BetaBBO__plane = 'arb'
        self._BetaBBO__theta_rad = 'var'
        self._BetaBBO__phi_rad = 'arb'

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
        # dn/dT
        self._dndT_o = -16.6e-6 #/degC
        self._dndT_e = -9.3e-6 #/degC
    
    @property
    def plane(self):
        return self._BetaBBO__plane
    @property
    def theta_rad(self):
        return self._BetaBBO__theta_rad
    @property
    def phi_rad(self):
        return self._BetaBBO__phi_rad
    @property
    def angles(self):
        msg =  ["plane = %s" % self._BetaBBO__plane]
        msg += ["theta_rad = %s" % self._BetaBBO__theta_rad]
        msg += ["phi_rad = %s" % self._BetaBBO__phi_rad]
        print("\n".join(msg))
    @property
    def symbols(self):
        return [wl, theta, phi, T]
    
    @property
    def constants(self):
        msg = ["A_o = %g" % self._A_o]
        msg += ["B_o = %g" % self._B_o]
        msg += ["C_o = %g" % self._C_o]
        msg += ["D_o = %g" % self._D_o]
        msg += ["A_e = %g" % self._A_e]
        msg += ["B_e = %g" % self._B_e]
        msg += ["C_e = %g" % self._C_e]
        msg += ["D_e = %g" % self._D_e]
        msg += ["dn_o/dT = %g" % self._dndT_o]
        msg += ["dn_e/dT = %g" % self._dndT_e]
        print("\n".join(msg))
    
    def n_o_expr(self):
        """ Sympy expression, dispersion formula for o-ray """
        return sympy.sqrt(self._A_o + self._B_o / (wl**2 - self._C_o) - self._D_o * wl**2) + self._dndT_o * (T - 20)
    
    def n_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-ray """
        return sympy.sqrt(self._A_e + self._B_e / (wl**2 - self._C_e) - self._D_e * wl**2) + self._dndT_e * (T - 20)

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
    
    def n(self, wl_um, theta_rad, T_degC, pol='o'):
        """
        Refractive index as a function of wavelength, theta and phi angles for each eigen polarization of light.

        input
        ------
        wl_um     :  float or array_like, wavelength in um
        theta_rad :  float or array_like, 0 to pi radians
        T_degC    :  float or array_like, temperature of crystal in degree C.
        pol       :  {'o', 'e'}, optional, polarization of light

        return
        -------
        Refractive index, float
        """
        return super().n(wl_um, theta_rad, 0.0, T_degC, pol=pol)

    def dn_wl(self, wl_um, theta_rad, T_degC, pol='o'):
        return super().dn_wl(wl_um, theta_rad, 0.0, T_degC, pol=pol)
    
    def d2n_wl(self, wl_um, theta_rad, T_degC, pol='o'):
        return super().d2n_wl(wl_um, theta_rad, 0.0, T_degC, pol=pol)

    def d3n_wl(self, wl_um, theta_rad, T_degC, pol='o'):
        return super().d3n_wl(wl_um, theta_rad, 0.0, T_degC, pol=pol)
    
    def GD(self, wl_um, theta_rad, T_degC, pol='o'):
        """Group Delay [fs/mm]"""
        return super().GD(wl_um, theta_rad, 0.0, T_degC, pol=pol)
    
    def GV(self, wl_um, theta_rad, T_degC, pol='o'):
        """Group Velocity [um/fs]"""
        return super().GV(wl_um, theta_rad, 0.0, T_degC, pol=pol)
    
    def ng(self, wl_um, theta_rad, T_degC, pol='o'):
        """Group index, c/Group velocity"""
        return super().ng(wl_um, theta_rad, 0.0, T_degC, pol=pol)
    
    def GVD(self, wl_um, theta_rad, T_degC, pol='o'):
        """Group Delay Dispersion [fs^2/mm]"""
        return super().GVD(wl_um, theta_rad, 0.0, T_degC, pol=pol)
    
    def TOD(self, wl_um, theta_rad, T_degC, pol='o'):
        """Third Order Dispersion [fs^3/mm]"""
        return super().TOD(wl_um, theta_rad, 0.0, T_degC, pol=pol)