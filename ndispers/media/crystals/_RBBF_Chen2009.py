import sympy

from ndispers._baseclass import Medium, T, phi, theta, wl
from ndispers.helper import vars2


class RBBF(Medium):
    """
    RBBF (RbBe2BO3F2, Rubidium Beryllium Borate Fluoride) crystal

    - Point group : 32  (D_3)
    - Crystal system : Trigonal
    - Dielectic principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Tranparency range : 0.16 to 3.55 µm
    - Hardness: 2.9 on the Mohs scale
    - Highly stable in air and even in hot water at 100°C or in acids such as HNO3 and HCl

    Sellmeier equation
    ------------------
    n_o^2 = 1 + 1.18675λ²/(λ² - 0.00750) - 0.00910λ²  (λ is in µm)
    n_e^2 = 1 + 0.97530λ²/(λ² - 0.00665) - 0.00145λ²  (λ is in µm)
    
    Validity range
    ---------------
    Deep UV to near infrared

    Ref
    ----
    Chen, C., Wu, Y., Li, Y., Wang, J., Wu, B., Jiang, M., Zhang, G., & Ye, N. (2009). Growth, properties, and application to nonlinear optics of a nonlinear optical crystal: RbBe2BO3F2. Journal of the Optical Society of America B, 26(8), 1519-1525. https://opg.optica.org/josab/abstract.cfm?uri=josab-26-8-1519
    
    Example
    -------
    >>> rbbf = ndispers.media.crystals.RBBF()
    >>> rbbf.n(0.6, 0.5*pi, 25, pol='e') # args: (wl_um, theta_rad, T_degC, pol)
    
    """
    __slots__ = ["_RBBF__plane", "_RBBF__theta_rad", "_RBBF__phi_rad",
                 "_B_o", "_C_o", "_D_o", 
                 "_B_e", "_C_e", "_D_e",
                 "_dndT_o", "_dndT_e"]

    def __init__(self):
        super().__init__()
        self._RBBF__plane = 'arb'
        self._RBBF__theta_rad = 'var'
        self._RBBF__phi_rad = 'arb'

        """ Constants of dispersion formula """
        # For ordinary ray (from the paper's equation: n_o^2 = 1 + 1.18675λ²/(λ² - 0.00750) - 0.00910λ²)
        self._B_o = 1.18675  # Numerator of first term
        self._C_o = 0.00750  # Denominator constant of first term
        self._D_o = 0.00910  # Coefficient of λ² term
        
        # For extraordinary ray (from the paper's equation: n_e^2 = 1 + 0.97530λ²/(λ² - 0.00665) - 0.00145λ²)
        self._B_e = 0.97530  # Numerator of first term
        self._C_e = 0.00665  # Denominator constant of first term
        self._D_e = 0.00145  # Coefficient of λ² term
        
        # dn/dT
        # not accessible from the paper
        self._dndT_o = 0.0 #/degC
        self._dndT_e = 0.0 #/degC
    
    @property
    def plane(self):
        return self._RBBF__plane

    @property
    def theta_rad(self):
        return self._RBBF__theta_rad

    @property
    def phi_rad(self):
        return self._RBBF__phi_rad

    @property
    def constants(self):
        print(vars2(self))
    
    @property
    def symbols(self):
        return [wl, theta, phi, T]

    def n_o_expr(self):
        """ Sympy expression, dispersion formula for o-wave """
        return sympy.sqrt(1 + self._B_o * wl**2 / (wl**2 - self._C_o) - self._D_o * wl**2) + self._dndT_o * (T - 20)
    
    def n_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-wave """
        return sympy.sqrt(1 + self._B_e * wl**2 / (wl**2 - self._C_e) - self._D_e * wl**2) + self._dndT_e * (T - 20)

    def n_expr(self, pol):
        """
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
        return super().woa_theta(wl_um, theta_rad, 0, T_degC, pol=pol)
    
    def woa_phi(self, wl_um, theta_rad, T_degC, pol='e'):
        return super().woa_phi(wl_um, theta_rad, 0, T_degC, pol=pol)

    def dndT(self, wl_um, theta_rad, T_degC, pol='o'):
        return super().dndT(wl_um, theta_rad, 0, T_degC, pol=pol)
