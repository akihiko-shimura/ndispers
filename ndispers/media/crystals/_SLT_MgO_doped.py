import sympy

from ndispers._baseclass import T, phi, theta, wl
from ndispers.groups import Uniax_neg_3m
from ndispers.helper import vars2


class SLT(Uniax_neg_3m):
    """
    0.5% MgO-doped stoichiometric Lithium tantalate (Li Ta O_3) crystal

    - Point group : 3m  (C_{3v})
    - Crystal system : Trigonal
    - Dielectic principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Tranparency range : 0.28 µm to 5.6 µm

    Sellmeier equation
    ------------------
    n(wl) = sqrt(a1_i + b1_i * f + (a2_i + b2_i * f)/(wl**2 - (a3_i + b3_i * f)**2) + (a4_i + b4_i * f)/(wl**2 - (a5_i + b5_i * f)**2) - a6_i * wl**2) for i=o,e
    f = (T - T0) * (T + T0 + 2 * 273.16) with T0 = 24.5 degC

    Ref
    ---
    - Sellmeier equation, ratio of d22 vs d33:
    Dolev, I., et al. "Linear and nonlinear optical properties of MgO: LiTaO 3." Applied Physics B 96 (2009): 423-432.
    - Second order nonlinear-optical coefficients, d33, d31:
    Shoji, I., et al. "Absolute scale of second-order nonlinear-optical coefficients," J. Opt. Soc. Am. B 14 (1997): 2268-2294
    
    Note
    ----
    d coefficients are of *congrunet* LT crystal, not of stoichiometric.

    Example
    -------
    >>> bbo = ndispers.media.crystals.BetaBBO_Eimerl1987()
    >>> bbo.n(0.6, 0.5*pi, 25, pol='e') # args: (wl_um, theta_rad, T_degC, pol)
    
    """
    __slots__ = ["_LT__plane", "_LT__theta_rad", "_LT__phi_rad",
                 "_a1_o", "_a2_o", "_a3_o", "_a4_o",  "_a5_o", "_a6_o",
                 "_a1_e", "_a2_e", "_a3_e", "_a4_e",  "_a5_e", "_a6_e",
                 "_b1_o", "_b2_o", "_b3_o", "_b4_o", "_b5_o", 
                 "_b1_e", "_b2_e", "_b3_e", "_b4_e", "_b5_e"]

    def __init__(self):
        super().__init__()
        self._LT__plane = 'arb'
        self._LT__theta_rad = 'var'
        self._LT__phi_rad = 'arb'

        """ Constants of dispersion formula """
        # 0.5% MgO-doped SLT
        ## ordinary ray
        self._a1_o = 4.5082
        self._a2_o = 0.084888
        self._a3_o = 0.19552
        self._a4_o = 1.1570
        self._a5_o = 8.2517
        self._a6_o = 0.0237
        self._b1_o = 2.0704e-8
        self._b2_o = 1.4449e-8
        self._b3_o = 1.5978e-8
        self._b4_o = 4.7686e-6
        self._b5_o = 1.1127e-5
        ## extraordinary ray
        self._a1_e = 4.5615
        self._a2_e = 0.08488
        self._a3_e = 0.1927
        self._a4_e = 5.5832
        self._a5_e = 8.3067
        self._a6_e = 0.021696
        self._b1_e = 4.782e-7
        self._b2_e = 3.0913e-8
        self._b3_e = 2.7326e-8
        self._b4_e = 1.4837e-5
        self._b5_e = 1.3647e-7
        # Second-order nonlinear optical coefficients
        self._d33_1064shg = 13.8 #pm/V
        self._d31_1064shg = 0.85 #pm/V
        self._d22_1064shg = 0.12 * 13.8 #pm/V
    
    @property
    def plane(self):
        return self._LT__plane

    @property
    def theta_rad(self):
        return self._LT__theta_rad

    @property
    def phi_rad(self):
        return self._LT__phi_rad

    @property
    def symbols(self):
        return [wl, theta, phi, T]
    
    @property
    def constants(self):
        print(vars2(self))
    
    def n_o_expr(self):
        """ Sympy expression, dispersion formula for o-wave """
        return sympy.sqrt( self._a1_o + self._b1_o * self.f_expr() + \
            (self._a2_o + self._b2_o * self.f_expr()) / (wl**2 - (self._a3_o + self._b3_o * self.f_expr())**2) + \
                (self._a4_o + self._b4_o * self.f_expr()) / (wl**2 - (self._a5_o + self._b5_o * self.f_expr())**2) - self._a6_o * wl**2 )
                
    def n_e_expr(self):
        """ Sympy expression, dispersion formula for o-wave """
        return sympy.sqrt( self._a1_e + self._b1_e * self.f_expr() + \
            (self._a2_e + self._b2_e * self.f_expr()) / (wl**2 - (self._a3_e + self._b3_e * self.f_expr())**2) + \
                (self._a4_e + self._b4_e * self.f_expr()) / (wl**2 - (self._a5_e + self._b5_e * self.f_expr())**2) - self._a6_e * wl**2 )

    def f_expr(self):
        return (T - 24.5) * (T + 24.5 + 2 * 273.16)

    def n_expr(self, pol):
        """"
        Sympy expression, 
        dispersion formula of a general ray with an angle theta to optic axis. If theta = 0, this expression reduces to 'n_o_expr'.

        n(theta) = n_e / sqrt( sin(theta)**2 + (n_e/n_o)**2 * cos(theta)**2 )
        """
        if pol == 'o':
            return self.n_o_expr()
        elif pol == 'e':
            return self.n_e_expr() / sympy.sqrt( sympy.sin(theta)**2 + (self.n_e_expr()/self.n_o_expr())**2 * sympy.cos(theta)**2 )
        else:
            raise ValueError("pol = '%s' must be 'o' or 'e'" % pol)
    
    def n(self, wl_um, theta_rad, T_degC, pol='e'):
        """
        Refractive index as a function of wavelength, theta or phi angles for each eigen polarization of light.

        input
        -----
        wl_um     :  float or array_like, wavelength in µm
        theta_rad :  float or array_like, 0 to pi radians
        T_degC    :  float or array_like, temperature of crystal in degree C.
        pol       :  {'e'}, optional, polarization of light

        return
        -------
        numpy.array
        
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
        """Group Velocity [µm/fs]"""
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
    
    #------------------------------------------------------------------------------------------
    # Wavelength dependence of second-order nonlinear coefficients estimated from Miller's rule
    #------------------------------------------------------------------------------------------
    def d22_sfg(self, wl1o, wl2o, T_degC):
        return super().d22_sfg(wl1o, wl2o, T_degC, delta22=self.delta22(self._d22_1064shg, 1.064, 1.064, T_degC))

    def d31_sfg(self, wl1o, wl2o, T_degC):
        return super().d31_sfg(wl1o, wl2o, T_degC, delta31=self.delta31(self._d31_1064shg, 1.064, 1.064, T_degC))