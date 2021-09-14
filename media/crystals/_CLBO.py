import sympy
from sympy.utilities import lambdify
from ndispers._baseclass import Medium, wl, phi, theta, T, pi

class CLBO(Medium):
    """
    CLBO (Cs Li B_6 O_10, Cesium Lithium Borate) crystal

    - Point group : 4m2
    - Crystal system : Tetragonal
    - Dielectic principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Tranparency range : 0.18 - 2.75 um

    Dispersion formula of refractive index
    ---------------------------------------
    n(wl_um) = sqrt(A_i + B_i/(wl_um**2 - C_i) - D_i * wl_um**2)  for i = o, e
    dn_o/dT = (At_o + Bt_o/wl)*1e-6
    dn_e/dT = (At_e + Bt_e/wl + Ct_e/wl**2 + Dt_e/wl**3)*1e-6
    
    Validity range
    ---------------
    0.1914 - 2.09 um
    dn/dT : 0.2128 - 1.3382 um

    Ref
    ----
    Umemura, N., Yoshida, K., Kamimura, T., Mori, Y., Sasaki, T., & Kato, K. "New data on the phase-matching properties of CsLiB6O10." Advanced Solid State Lasers. Optical Society of America, 1999.

    Usage
    ------
    >>> clbo = ndispers.media.crystals.CLBO()
    >>> clbo.n(0.6, 0, pol='o') # for o-ray, n does not depend on theta.
    >>> clbo.n(0.6, 0.5*pi, pol='e') # along z-axis, it is pure e-ray.
    >>> clbo.n(0.6, 0.23*pi, pol='e')
    >>> clbo.n(0.6, 0*pi, pol='e') # for theta = 0 rad, it corresponds to o-ray.
    >>> clbo.GVD(0.6, 0.23*pi, pol='e')

    @author: Akihiko Shimura
    """
    __slots__ = ["_CLBO__plane", "_CLBO__theta_rad", "_CLBO__phi_rad",
                 "_A_o", "_B_o", "_C_o", "_D_o", 
                 "_A_e", "_B_e", "_C_e", "_D_e",
                 "_At_o", "_Bt_o",
                 "_At_e", "_Bt_e", "_Ct_e", "_Dt_e"]

    def __init__(self):
        super().__init__()
        self._CLBO__plane = 'arb'
        self._CLBO__theta_rad = 'var'
        self._CLBO__phi_rad = 'arb'

        """ Constants of dispersion formula """
        # For ordinary ray
        self._A_o = 2.2104
        self._B_o = 0.01018
        self._C_o = 0.01424
        self._D_o = 0.01258
        # For extraordinary ray
        self._A_e = 2.0588
        self._B_e = 0.00838
        self._C_e = 0.01363
        self._D_e = 0.00607
        # dn/dT
        self._At_o = -12.48 #1/K
        self._Bt_o = -0.328 #um/K
        self._At_o = -8.36 #1/K
        self._Bt_e = 0.047 #um/K
        self._Ct_e = 0.039 #um^2/K
        self._Dt_e = 0.014 #um^3/K
    
    @property
    def plane(self):
        return self._CLBO__plane
    @property
    def theta_rad(self):
        return self._CLBO__theta_rad
    @property
    def phi_rad(self):
        return self._CLBO__phi_rad
    @property
    def angles(self):
        msg =  ["plane = %s" % self._CLBO__plane]
        msg += ["theta_rad = %s" % self._CLBO__theta_rad]
        msg += ["phi_rad = %s" % self._CLBO__phi_rad]
        print("\n".join(msg))
    @property
    def symbols(self):
        return [wl, theta, phi, T]

    @property
    def constants(self):
        msg  = ["A_o = %g" % self._A_o]
        msg += ["B_o = %g" % self._B_o]
        msg += ["C_o = %g" % self._C_o]
        msg += ["D_o = %g" % self._D_o]
        msg += ["A_e = %g" % self._A_e]
        msg += ["B_e = %g" % self._B_e]
        msg += ["C_e = %g" % self._C_e]
        msg += ["D_e = %g" % self._D_e]
        msg += ["At_o = %g" % self._At_o]
        msg += ["Bt_o = %g" % self._Bt_o]
        msg += ["At_e = %g" % self._At_e]
        msg += ["Bt_e = %g" % self._Bt_e]
        msg += ["Ct_e = %g" % self._Ct_e]
        msg += ["Dt_e = %g" % self._Dt_e]
        print("\n".join(msg))

    def dndT_o(self):
        return  (self._At_o + self.Bt_o/wl)*1e-6
    
    def dndT_e(self):
        return  (self._At_e + self._Bt_e/wl + self._Ct_e/wl**2 + self._Dt_e/wl**3)*1e-6
    
    def n_o_expr(self):
        """ Sympy expression, dispersion formula for o-ray """
        return sympy.sqrt(self._A_o + self._B_o / (wl**2 - self._C_o) - self._D_o * wl**2) + self.dndT_o() * (T - 20.0)
    
    def n_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-ray """
        return sympy.sqrt(self._A_e + self._B_e / (wl**2 - self._C_e) - self._D_e * wl**2) + self.dndT_e() * (T - 20.0)

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