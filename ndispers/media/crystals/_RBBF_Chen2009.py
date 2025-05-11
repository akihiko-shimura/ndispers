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
    
    Thermo-optic dispersion equation
    ------------------------------
    dn/dT = (A/λ³ + B/λ² + C/λ + D) × 10⁻⁶  (λ is in µm)
    
    For ordinary ray (n_o):
    A = 0.099911, B = -0.553474, C = 1.454609, D = -13.260115
    
    For extraordinary ray (n_e):
    A = 0.285633, B = -2.482927, C = 6.916728, D = -16.153736
    
    Validity range
    ---------------
    Sellmeier equation: Deep UV to near infrared
    Thermo-optic equation: 0.194µm ≤ λ ≤ 1.014µm
    
    Implementation notes
    ------------------
    RBBF is a uniaxial crystal where the refractive index doesn't depend on the azimuthal angle phi.
    In all method calls to the parent class, a value of 0 is passed for the phi_rad parameter.

    Ref
    ----
    - Chen, C., Wu, Y., Li, Y., Wang, J., Wu, B., Jiang, M., Zhang, G., & Ye, N. (2009). Growth, properties, and application to nonlinear optics of a nonlinear optical crystal: RbBe2BO3F2. Journal of the Optical Society of America B, 26(8), 1519-1525. https://opg.optica.org/josab/abstract.cfm?uri=josab-26-8-1519
    - Zhai, N., Wang, L., Liu, L., Wang, X., Zhu, Y., & Chen, C. (2013). Measurement of thermal refractive index coefficients of nonlinear optical crystal RbBe2BO3F2. Optical Materials, 36(2), 333-336.

    Example
    -------
    >>> rbbf = ndispers.media.crystals.RBBF()
    >>> rbbf.n(0.6, 0.5*pi, 25, pol='e') # args: (wl_um, theta_rad, T_degC, pol)
    
    """
    __slots__ = ["_RBBF__plane", "_RBBF__theta_rad", "_RBBF__phi_rad",
                 "_B_o", "_C_o", "_D_o", 
                 "_B_e", "_C_e", "_D_e",
                 "_dndT_o_A", "_dndT_o_B", "_dndT_o_C", "_dndT_o_D",
                 "_dndT_e_A", "_dndT_e_B", "_dndT_e_C", "_dndT_e_D"]

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
        
        # dn/dT coefficients from Zhai et al. (2013)
        # Formula: dn/dT = (A/λ³ + B/λ² + C/λ + D) × 10⁻⁶
        # Valid for 0.194μm ≤ λ ≤ 1.014μm
        # Ordinary ray
        self._dndT_o_A = 0.099911  # coefficient of 1/λ³
        self._dndT_o_B = -0.553474  # coefficient of 1/λ²
        self._dndT_o_C = 1.454609  # coefficient of 1/λ
        self._dndT_o_D = -13.260115  # coefficient of constant term
        
        # Extraordinary ray
        self._dndT_e_A = 0.285633  # coefficient of 1/λ³
        self._dndT_e_B = -2.482927  # coefficient of 1/λ²
        self._dndT_e_C = 6.916728  # coefficient of 1/λ
        self._dndT_e_D = -16.153736  # coefficient of constant term
    
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

    def dndT_o_expr(self):
        """ Sympy expression for thermo-optic coefficient of o-wave (dn/dT) """
        # Formula: dn/dT = (A/λ³ + B/λ² + C/λ + D) × 10⁻⁶
        return (self._dndT_o_A / wl**3 + self._dndT_o_B / wl**2 + self._dndT_o_C / wl + self._dndT_o_D) * 1e-6
    
    def n_o_expr(self):
        """ Sympy expression, dispersion formula for o-wave """
        n_o = sympy.sqrt(1 + self._B_o * wl**2 / (wl**2 - self._C_o) - self._D_o * wl**2)
        # Add temperature dependence
        return n_o + self.dndT_o_expr() * (T - 24)  # T_ref = 24°C in the paper
    
    def dndT_e_expr(self):
        """ Sympy expression for thermo-optic coefficient of e-wave (dn/dT) """
        # Formula: dn/dT = (A/λ³ + B/λ² + C/λ + D) × 10⁻⁶
        return (self._dndT_e_A / wl**3 + self._dndT_e_B / wl**2 + self._dndT_e_C / wl + self._dndT_e_D) * 1e-6
        
    def n_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-wave """
        n_e = sympy.sqrt(1 + self._B_e * wl**2 / (wl**2 - self._C_e) - self._D_e * wl**2)
        # Add temperature dependence
        return n_e + self.dndT_e_expr() * (T - 24)  # T_ref = 24°C in the paper

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
        
        Implementation notes
        -----
        When calling the parent class method, we pass 0 for phi_rad since RBBF is a uniaxial crystal
        where the refractive index doesn't depend on the azimuthal angle phi.
        The argument order in the parent class is (wl_um, theta_rad, phi_rad, T_degC, pol).
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
        """Thermo-optic coefficient (dn/dT) [1/°C]
        
        Formula: dn/dT = (A/λ³ + B/λ² + C/λ + D) × 10⁻⁶
        Valid for 0.194μm ≤ λ ≤ 1.014μm
        
        input
        ------
        wl_um     :  float or array_like, wavelength in µm
        theta_rad :  float or array_like, 0 to pi radians
        T_degC    :  float or array_like, temperature of crystal in degree C.
        pol       :  {'o', 'e'}, optional, polarization of light
        
        return
        -------
        dn/dT in 1/°C, float or array_like
        """
        return super().dndT(wl_um, theta_rad, 0, T_degC, pol=pol)
