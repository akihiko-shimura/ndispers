import sympy
from sympy.utilities import lambdify
from ndispers._baseclass import Medium, wl, phi, theta, T
from ndispers.helper import vars2

class BetaBBO(Medium):
    """
    beta-BBO (beta-Ba B_2 O_4) crystal

    - Point group : 3m
    - Crystal system : Trigonal
    - Dielectic principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Tranparency range : 1.9 to 2.6 um

    Dispersion formula for refractive index
    ---------------------------------------
    n(wl) = sqrt(A_i + B_i/(1 - C_i/wl**2) + D_i/(1 - E_i/wl**2)) + dn/dT * (T -20)

    Thermo-optic coefficient
    ------------------------
    dn/dT = (G_i * R_i + H_i * R_i**2) / 2*n_i for i=o, e
    R_i = wl**2/(wl**2 - wl0_i**2), wl0_i = 0.0652 for i=o, wl0=0.0730 for i=e
    
    Validity range
    --------------

    Ref
    ----
    Ghosh, Gorachand. "Temperature dispersion of refractive indices in β‐BaB2O4 and LiB3O5 crystals for nonlinear optical devices." Journal of applied physics 78.11 (1995): 6752-6760.

    Example
    -------
    >>> bbo = ndispers.media.crystals.BetaBBO_Eimerl1987()
    >>> bbo.n(0.6, 0, 40, pol='o') # args: (wl_um, theta_rad, T_degC, pol)
    >>> bbo.n(0.6, 0.5*pi, 40, pol='e') # along z-axis, it is pure e-ray.
    >>> bbo.n(0.6, 0*pi, 40, pol='e') # for theta = 0 rad, it corresponds to o-ray.
    >>> bbo.GVD(0.6, 0.3*pi, 40, pol='e')
    >>> bbo.pmAngles_sfg(1.064, 1.064, 40, deg=True)
    {'ooe': {'theta': [22.895], 'phi': None},
     'eeo': {'theta': [], 'phi': None},
     'oee': {'theta': [32.575], 'phi': None},
     'eoe': {'theta': [32.575], 'phi': None},
     'eoo': {'theta': [], 'phi': None},
     'oeo': {'theta': [], 'phi': None}}

    
    """
    __slots__ = ["_BetaBBO__plane", "_BetaBBO__theta_rad", "_BetaBBO__phi_rad",
                 "_A_o", "_B_o", "_C_o", "_D_o", "_E_o",
                 "_A_e", "_B_e", "_C_e", "_D_e", "_E_e",
                 "_G_o", "_H_o", "_R_o",
                 "_G_e", "_H_e", "_R_e"]
                 
    def __init__(self):
        super().__init__()
        self._BetaBBO__plane = 'arb'
        self._BetaBBO__theta_rad = 'var'
        self._BetaBBO__phi_rad = 'arb'

        """ Constants of dispersion formula """
        # For ordinary ray
        self._A_o = 1.7018379
        self._B_o = 1.0357554
        self._C_o = 0.018003440
        self._D_o = 1.2479989
        self._E_o = 91
        # For extraordinary ray
        self._A_e = 1.5920433
        self._B_e = 0.7816893
        self._C_e = 0.016067891
        self._D_e = 0.8403893
        self._E_e = 91
        # dn/dT
        self._G_o = -19.3007e-6
        self._G_e = -141.421e-6
        self._H_o = -34.9683e-6
        self._H_e = 110.8630e-6
        self._R_o = wl**2/(wl**2 - 0.0652**2)
        self._R_e = wl**2/(wl**2 - 0.0730**2)
    
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
    def symbols(self):
        return [wl, theta, phi, T]
    
    @property
    def constants(self):
        print(vars2(self))
    
    def _n_T20_o_expr(self):
        """ Sympy expression, dispersion formula for o-ray at 20degC """
        return sympy.sqrt(self._A_o + self._B_o / (1 - self._C_o/wl**2) + self._D_o/(1 - self._E_o/wl**2))
    
    def _n_T20_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-ray at 20degC """
        return sympy.sqrt(self._A_e + self._B_e / (1 - self._C_e/wl**2) + self._D_e/(1 - self._E_e/wl**2))
    
    def dndT_o_expr(self):
        return (self._G_o * self._R_o + self._H_o * self._R_o**2) / (2*self._n_T20_o_expr())
    
    def dndT_e_expr(self):
        return (self._G_e * self._R_e + self._H_e * self._R_e**2) / (2*self._n_T20_e_expr())
    
    def n_o_expr(self):
        return self._n_T20_o_expr() + self.dndT_o_expr() * (T - 20)
    
    def n_e_expr(self):
        return self._n_T20_e_expr() + self.dndT_e_expr() * (T - 20)

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
        Refractive index, float or array_like
        """
        return super().n(wl_um, theta_rad, 0, T_degC, pol=pol)

    def dn_wl(self, wl_um, theta_rad, T_degC, pol='o'):
        return super().dn_wl(wl_um, theta_rad, 0, T_degC, pol=pol)
    
    def d2n_wl(self, wl_um, theta_rad, T_degC, pol='o'):
        return super().d2n_wl(wl_um, theta_rad, 0, T_degC, pol=pol)

    def d3n_wl(self, wl_um, theta_rad, T_degC, pol='o'):
        return super().d3n_wl(wl_um, theta_rad, 0, T_degC, pol=pol)
    
    def GD(self, wl_um, theta_rad, T_degC, pol='o'):
        """Group Delay [fs/mm]"""
        return super().GD(wl_um, theta_rad, 0, T_degC, pol=pol)
    
    def GV(self, wl_um, theta_rad, T_degC, pol='o'):
        """Group Velocity [um/fs]"""
        return super().GV(wl_um, theta_rad, 0, T_degC, pol=pol)
    
    def ng(self, wl_um, theta_rad, T_degC, pol='o'):
        """Group index, c/Group velocity"""
        return super().ng(wl_um, theta_rad, 0, T_degC, pol=pol)
    
    def GVD(self, wl_um, theta_rad, T_degC, pol='o'):
        """Group Delay Dispersion [fs^2/mm]"""
        return super().GVD(wl_um, theta_rad, 0, T_degC, pol=pol)
    
    def TOD(self, wl_um, theta_rad, T_degC, pol='o'):
        """Third Order Dispersion [fs^3/mm]"""
        return super().TOD(wl_um, theta_rad, 0, T_degC, pol=pol)
    
    def woa_theta(self, wl_um, theta_rad, T_degC, pol='e'):
        """ Polar walk-off angle [rad] """
        return super().woa_theta(wl_um, theta_rad, 0, T_degC, pol=pol)
    
    def woa_phi(self, wl_um, theta_rad, T_degC, pol='e'):
        """ Azimuthal walk-off angle [rad] """
        return super().woa_phi(wl_um, theta_rad, 0, T_degC, pol=pol)
    
    def dndT(self, wl_um, theta_rad, T_degC, pol='o'):
        return super().dndT(wl_um, theta_rad, 0, T_degC, pol=pol)