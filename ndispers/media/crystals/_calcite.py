import sympy
from ndispers._baseclass import Medium, wl, phi, theta
from ndispers.helper import vars2

class Calcite(Medium):
    """
    calcite (Ca C O_3) crystal

    - Point group : -3m  (D_{3d})
    - Crystal system : Trigonal
    - Dielectic principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Tranparency range : 

    Sellmeier equation
    ------------------
    n(wl) 
    = sqrt(1 + A1_o * wl**2 / (wl**2 - B1_o**2) + A2_o * wl**2 / (wl**2 - B2_o**2) + A3_o * wl**2 / (wl**2 - B3_o**2) + A4_o * wl**2 / (wl**2 - B4_o**2)) for o-wave
    = sqrt(1 + A1_e * wl**2 / (wl**2 - B1_e**2) + A2_e * wl**2 / (wl**2 - B2_e**2) + A3_e * wl**2 / (wl**2 - B3_e**2)) for e-wave
    
    Validity range
    ---------------
    0.2 to 2.2 µm for o-wave
    0.2 to 3.3 µm for e-wave

    Ref
    ----
    Handbook of Optics: Devices, Measurements, and Properties, Volume II, by Michael Bass (ed),
    Chapter 33: PROPERTIES OF CRYSTALS AND GLASSES, William J. Tropf, Michael E. Thomas, and Terry J. Harris
    
    Example
    -------
    >>> bbo = ndispers.media.crystals.BetaBBO_Eimerl1987()
    >>> bbo.n(0.6, 0.5*pi, 25, pol='e') # args: (wl_um, theta_rad, T_degC, pol)
    
    """
    __slots__ = ["_Calcite__plane", "_Calcite__theta_rad", "_Calcite__phi_rad",
                 "_A1_o", "_B1_o", "_A2_o", "_B2_o", "_A3_o", "_B3_o", "_A4_o", "_B4_o",
                 "_A1_e", "_B1_e", "_A2_e", "_B2_e", "_A3_e", "_B3_e"]

    def __init__(self):
        super().__init__()
        self._Calcite__plane = 'arb'
        self._Calcite__theta_rad = 'var'
        self._Calcite__phi_rad = 'arb'

        """ Constants of dispersion formula """
        # For ordinary ray
        self._A1_o = 0.8559
        self._B1_o = 0.0588
        self._A2_o = 0.8391
        self._B2_o = 0.141
        self._A3_o = 0.0009
        self._B3_o = 0.197
        self._A4_o = 0.6845
        self._B4_o = 7.005
        # For extraordinary ray
        self._A1_e = 1.0856
        self._B1_e = 0.07897
        self._A2_e = 0.0988
        self._B2_e = 0.142
        self._A3_e = 0.317
        self._B3_e = 11.468
    
    @property
    def plane(self):
        return self._Calcite__plane

    @property
    def theta_rad(self):
        return self._Calcite__theta_rad

    @property
    def phi_rad(self):
        return self._Calcite__phi_rad

    @property
    def symbols(self):
        return [wl, theta, phi]
    
    @property
    def constants(self):
        print(vars2(self))
    
    def n_o_expr(self):
        """ Sympy expression, dispersion formula for o-wave """
        return sympy.sqrt(1 + self._A1_o * wl**2 / (wl**2 - self._B1_o**2) + self._A2_o * wl**2 / (wl**2 - self._B2_o**2) + self._A3_o * wl**2 / (wl**2 - self._B3_o**2) + self._A4_o * wl**2 / (wl**2 - self._B4_o**2))
    
    def n_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-wave """
        return sympy.sqrt(1 + self._A1_e * wl**2 / (wl**2 - self._B1_e**2) + self._A2_e * wl**2 / (wl**2 - self._B2_e**2) + self._A3_e * wl**2 / (wl**2 - self._B3_e**2))

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
    
    def n(self, wl_um, theta_rad, pol='o'):
        """
        Refractive index as a function of wavelength, theta and phi angles for each eigen polarization of light.

        input
        ------
        wl_um     :  float, wavelength in µm
        theta_rad :  float, 0 to pi radians
        pol       :  str, 'o' or 'e', polarization of light

        return
        -------
        Refractive index, float or array_like
        
        """
        return super().n(wl_um, theta_rad, 0, pol=pol)

    def dn_wl(self, wl_um, theta_rad, pol='o'):
        return super().dn_wl(wl_um, theta_rad, 0, pol=pol)
    
    def d2n_wl(self, wl_um, theta_rad, pol='o'):
        return super().d2n_wl(wl_um, theta_rad, 0, pol=pol)

    def d3n_wl(self, wl_um, theta_rad, pol='o'):
        return super().d3n_wl(wl_um, theta_rad, 0, pol=pol)
    
    def GD(self, wl_um, theta_rad, pol='o'):
        """Group Delay [fs/mm]"""
        return super().GD(wl_um, theta_rad, 0, pol=pol)
    
    def GV(self, wl_um, theta_rad, pol='o'):
        """Group Velocity [µm/fs]"""
        return super().GV(wl_um, theta_rad, 0, pol=pol)
    
    def ng(self, wl_um, theta_rad, pol='o'):
        """Group index, c/Group velocity"""
        return super().ng(wl_um, theta_rad, 0, pol=pol)
    
    def GVD(self, wl_um, theta_rad, pol='o'):
        """Group Delay Dispersion [fs^2/mm]"""
        return super().GVD(wl_um, theta_rad, 0, pol=pol)
    
    def TOD(self, wl_um, theta_rad, pol='o'):
        """Third Order Dispersion [fs^3/mm]"""
        return super().TOD(wl_um, theta_rad, 0, pol=pol)
    
    def woa_theta(self, wl_um, theta_rad, T_degC, pol='e'):
        return super().woa_theta(wl_um, theta_rad, 0, T_degC, pol=pol)
    
    def woa_phi(self, wl_um, theta_rad, T_degC, pol='e'):
        return super().woa_phi(wl_um, theta_rad, 0, T_degC, pol=pol)