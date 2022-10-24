import sympy
from ndispers._baseclass import wl, phi, theta, T
from ndispers.groups import Uniax_neg_3m
from ndispers.helper import vars2

class BetaBBO(Uniax_neg_3m):
    """
    β-BBO (β-Ba B_2 O_4) crystal

    - Point group : 3m  (C_{3v})
    - Crystal system : Trigonal
    - Dielectic principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Tranparency range : 0.19 µm to 2.6 µm
    
    Sellmeier equation
    ------------------
    n(wl_um) = sqrt(1 + B1_i*wl**2/(wl**2 - C1_i) + B2_i*wl**2/(wl**2 - C2_i) + B3_i*wl**2/(wl**2 - C3_i)) + dn/dT * (T - 20)  for i = o, e
    
    Validity range
    --------------
    0.188 - 5.2 µm

    Ref
    ---
    Sellmeier equation:
    - Tamošauskas, Gintaras, et al. "Transmittance and phase matching of BBO crystal in the 3-5 μm range and its application for the characterization of mid-infrared laser pulses." Optical Materials Express 8.6 (2018): 1410-1418.
    dn/dT from Nikogosyan, D. N. "Beta barium borate (BBO)." Applied Physics A 52.6 (1991): 359-368.

    Nonlinear optical coefficients:
    - Shoji, Ichiro, et al. "Absolute measurement of second-order nonlinear-optical coefficients of β-BaB2O4 for visible to ultraviolet second-harmonic wavelengths." JOSA B 16.4 (1999): 620-624.

    Example
    -------
    >>> bbo = ndispers.media.crystals.BetaBBO_Eimerl1987()
    >>> bbo.n(0.6, 0.5*pi, 25, pol='e') # args: (wl_um, theta_rad, T_degC, pol)
    
    """
    __slots__ = ["_BetaBBO__plane", "_BetaBBO__theta_rad", "_BetaBBO__phi_rad",
                 "_B1_o", "_C1_o", "_B2_o", "_C2_o", "_B3_o", "_C3_o",
                 "_B1_e", "_C1_e", "_B2_e", "_C2_e", "_B3_e", "_C3_e",
                 "_dndT_o", "_dndT_e",
                 "_d31_1064shg", "_d22_1064shg"]

    def __init__(self):
        super().__init__()
        self._BetaBBO__plane = 'arb'
        self._BetaBBO__theta_rad = 'var'
        self._BetaBBO__phi_rad = 'arb'

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
        # dn/dT
        self._dndT_o = -16.6e-6 #/degC
        self._dndT_e = -9.3e-6 #/degC
         # Second-order nonlinear optical coefficients
        self._d31_1064shg = 0.04 #pm/V
        self._d22_1064shg = 2.2 #pm/V
    
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
    def constants(self):
        print(vars2(self))

    @property
    def symbols(self):
        return [wl, theta, phi, T]
    
    @property
    def constants(self):
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
        msg += ["dn_o/dT = %g" % self._dndT_o]
        msg += ["dn_e/dT = %g" % self._dndT_e]
        print("\n".join(msg))
    
    def n_o_expr(self):
        """ Sympy expression, dispersion formula for o-wave """
        return sympy.sqrt(1.0 + self._B1_o * wl**2/ (wl**2 - self._C1_o) + self._B2_o * wl**2/ (wl**2 - self._C2_o) + self._B3_o * wl**2/ (wl**2 - self._C3_o)) + self._dndT_o * (T - 20)
    
    def n_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-wave """
        return sympy.sqrt(1.0 + self._B1_e * wl**2/ (wl**2 - self._C1_e) + self._B2_e * wl**2/ (wl**2 - self._C2_e) + self._B3_e * wl**2/ (wl**2 - self._C3_e)) + self._dndT_e * (T - 20)


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
        wl_um     :  float or array_like, wavelength in µm
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