import sympy
from ndispers._baseclass import wl, phi, theta, T
from ndispers.groups import Uniax_neg_32
from ndispers.helper import vars2

class KBBF(Uniax_neg_32):
    """
    KBBF (K Be_2 B O_3 F_2) crystal

    - Point group : 32  (D_2)
    - Crystal system : Trigonal
    - Dielectric principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Transparency range : 0.155 µm to 3.66 µm

    Sellmeier equation
    ------------------
    n(wl) = sqrt(A_i + B1_i * wl**2 / (wl**2 - C1_i**2) + B2_i * wl**2 / (wl**2 - C2_i**2) - D_i * wl**2)  + dndT_i * (T - 22) for i = o, e

    Validity range
    --------------
    Fitted to refractive indices measured from 0.2537 to 1.530 µm;
    validated against observed phase-matching angles for harmonic
    generation down to about 0.15 µm in the deep UV.

    Note
    ----
    No thermo-optic coefficients (dn/dT) have been reported in the literature
    for KBBF, so the temperature dependence cannot be computed: dndT_i above
    is set to zero, n is independent of T, and the dndT method returns 0.

    Ref
    ---
    Sellmeier equation:
    Li, R., Wang, L., Wang, X., Wang, G., & Chen, C. (2016). Dispersion relations of refractive indices suitable for KBe 2 BO 3 F 2 crystal deep-ultraviolet applications. Applied Optics, 55(36), 10423-10426.
    
    Nonlinear optical coefficients:
    Chen, C. T., Wang, G. L., Wang, X. Y., & Xu, Z. Y. (2009). Deep-UV nonlinear optical crystal KBe2BO3F2—discovery, growth, optical properties and applications. Applied Physics B, 97(1), 9-25.
    """
    __slots__ = ["_A_o", "_B1_o", "_B2_o", "_C1_o", "_C2_o",  "_D_o",
                 "_A_e", "_B1_e", "_B2_e", "_C1_e", "_C2_e",  "_D_e",
                 "_d11_1064shg", "_d14_1064shg",
                 "_dndT_o", "_dndT_e"]
                 
    def __init__(self):
        super().__init__()
        self._plane = 'arb'
        self._theta_rad = 'var'
        self._phi_rad = 'arb'

        """ Constants of dispersion formula """
        # For ordinary ray
        self._A_o =  1.024248
        self._B1_o = 0.9502782
        self._B2_o = 0.1960298
        self._C1_o = 0.0738546
        self._C2_o = 0.1298386
        self._D_o =  0.0113908
        # For extraordinary ray
        self._A_e =  0.9411543
        self._B1_e = 0.8684699
        self._B2_e = 0.1256642
        self._C1_e = 0.0646955
        self._C2_e = 0.1196215
        self._D_e =  0.0044736
        # dn/dT
        # (not known yet)
        self._dndT_o = 0.0 #/degC
        self._dndT_e = 0.0 #/degC
        # Second-order nonlinear optical coefficients
        self._d11_1064shg = 0.47 #pm/V
        self._d14_1064shg = 0.0 #pm/V
    

    
    
    def n_o_expr(self):
        """ Sympy expression, dispersion formula for o-wave """
        return sympy.sqrt(self._A_o + self._B1_o * wl**2 / (wl**2 - self._C1_o**2) + self._B2_o * wl**2 / (wl**2 - self._C2_o**2) - self._D_o * wl**2) + self._dndT_o * (T - 22)
    
    def n_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-wave """
        return sympy.sqrt(self._A_e + self._B1_e * wl**2 / (wl**2 - self._C1_e**2) + self._B2_e * wl**2 / (wl**2 - self._C2_e**2) - self._D_e * wl**2) + self._dndT_e * (T - 22)

    def n_expr(self, pol):
        """
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
    

    

    
    
    
    
    
    
    
    
    
    #------------------------------------------------------------------------------------------
    # Wavelength dependence of second-order nonlinear coefficients estimated from Miller's rule
    #------------------------------------------------------------------------------------------
    def d11_sfg(self, wl1o, wl2o, T_degC):
        return super().d11_sfg(wl1o, wl2o, T_degC, delta11=self.delta11(self._d11_1064shg, 1.064, 1.064, T_degC))

    def d14_sfg(self, wl1o, wl2o, T_degC):
        return super().d14_sfg(wl1o, wl2o, T_degC, delta14=self.delta14(self._d14_1064shg, 1.064, 1.064, T_degC))
    