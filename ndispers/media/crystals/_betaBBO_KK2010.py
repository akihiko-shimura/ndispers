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
    n(wl) = sqrt(A_i + B_i/(wl**2 - C_i) + D_i/(wl**2 - E_i)) + dn/dT * (T -20)

    Thermo-optic coefficient
    ------------------------
    dn/dT = F_i/wl**3 + G_i/wl**2 + H_i/wl + I_i  for i = o, e
    
    Validity range
    --------------
    o: 0.2048 to 3.22 um
    e: 0.1916 to 0.2048 um
    dn/dT: 0.195 to 1.618 um

    Ref
    ---
    Kato, K., N. Umemura, and T. Mikami. "Sellmeier and thermo-optic dispersion formulas for Beta-Ba B_2 O_4 (revisited)." Nonlinear Frequency Generation and Conversion: Materials, Devices, and Applications IX. Vol. 7582. International Society for Optics and Photonics, 2010.

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
                 "_F_o", "_G_o", "_H_o", "_I_o",
                 "_F_e", "_G_e", "_H_e", "_I_e"]
                 
    def __init__(self):
        super().__init__()
        self._BetaBBO__plane = 'arb'
        self._BetaBBO__theta_rad = 'var'
        self._BetaBBO__phi_rad = 'arb'

        """ Constants of dispersion formula """
        # For ordinary ray
        self._A_o = 3.63357
        self._B_o = 0.01878
        self._C_o = 0.01822
        self._D_o = 60.9129
        self._E_o = 67.8505
        # For extraordinary ray (best fit from 0.2048 to 3.22 um)
        self._A_e = 3.33469
        self._B_e = 0.01237
        self._C_e = 0.01647
        self._D_e = 79.0672
        self._E_e = 82.2919
        # For extraordinary ray (best fit from 0.1916 to 0.2048 um)
        # self._A_e = 3.38630
        # self._B_e = 0.00921
        # self._C_e = 0.02073
        # self._D_e = 79.0672
        # self._E_e = 82.2919
        # dn/dT (best fit from 0.195 to 1.618 um)
        self._F_o = -0.0137e-5
        self._F_e = 0.0413e-5
        self._G_o = 0.0607e-5
        self._G_e = -0.2119e-5
        self._H_o = -0.1334e-5
        self._H_e = 0.4408e-5
        self._I_o = -1.5287e-5
        self._I_e = -1.2749e-5
    
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
        return sympy.sqrt(self._A_o + self._B_o/(wl**2 - self._C_o) + self._D_o/(wl**2 - self._E_o))
    
    def _n_T20_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-ray at 20degC """
        return sympy.sqrt(self._A_e + self._B_e/(wl**2 - self._C_e) + self._D_e/(wl**2 - self._E_e))
    
    def dndT_o_expr(self):
        return self._F_o/wl**3 + self._G_o/wl**2 + self._I_o/wl + self._I_o

    def dndT_e_expr(self):
        return self._F_e/wl**3 + self._G_e/wl**2 + self._I_e/wl + self._I_e
    
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