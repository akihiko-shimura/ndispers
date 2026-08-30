from ndispers._sym import sympy

from ndispers._baseclass import T, phi, theta, wl
from ndispers.groups import Uniax_32
from ndispers.helper import vars2

class RBBF(Uniax_32):
    """
    RBBF (RbBe₂BO₃F₂, rubidium beryllium fluoroborate) crystal

    - Point group : 32  (D3)
    - Crystal system : Trigonal
    - Dielectric principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Transparency range : 0.16 to 3.55 µm
    - Hardness: 2.9 on the Mohs scale
    - Highly stable in air and even in hot water at 100°C or in acids such as HNO3 and HCl

    Sellmeier equation
    ------------------
        n_o^2 = 1 + 1.18675λ²/(λ² - 0.00750) - 0.00910λ²  (λ is in µm)
        n_e^2 = 1 + 0.97530λ²/(λ² - 0.00665) - 0.00145λ²  (λ is in µm)

    Thermo-optic coefficient
    ------------------------
        dn/dT = (A/λ³ + B/λ² + C/λ + D) × 10⁻⁶  (λ is in µm)

    For ordinary ray (n_o):
    A = 0.099911, B = -0.553474, C = 1.454609, D = -13.260115

    For extraordinary ray (n_e):
    A = 0.285633, B = -2.482927, C = 6.916728, D = -16.153736

    The temperature correction dn/dT * (T - 24) vanishes at 24 degC.

    Validity range
    ---------------
    Sellmeier equation: Deep UV to near infrared
    Thermo-optic equation: 0.194µm ≤ λ ≤ 1.014µm

    Ref
    ---
    Sellmeier equation:
      Chen, C., Wu, Y., Li, Y., et al. (2009). Deep UV nonlinear optical crystal: RbBe2(BO3)F2. JOSA B, 26(8), 1519-1525. https://doi.org/10.1364/josab.26.001519  (Sellmeier equation and, Section 5, the d11 coefficient)

    Thermo-optic coefficients:
      Zhai, N., Wang, L., Liu, L., et al. (2013). Measurement of thermal refractive index coefficients of nonlinear optical crystal RbBe2BO3F2. Optical Materials, 36(2), 333-336. https://doi.org/10.1016/j.optmat.2013.09.017
    """
    _wl_range = (0.16, 3.55)  # um, Sellmeier validity (see docstring)

    __slots__ = ["_B_o", "_C_o", "_D_o",
                 "_B_e", "_C_e", "_D_e",
                 "_dndT_o_A", "_dndT_o_B", "_dndT_o_C", "_dndT_o_D",
                 "_dndT_e_A", "_dndT_e_B", "_dndT_e_C", "_dndT_e_D"]

    _d_ref = {"d11": (0.45, 1.064, 1.064)}
    _d_note = ("Chen et al. 2009 (JOSA B), Maker fringes at 1.064 um relative to "
               "d36(KDP) = 0.39 pm/V. d14 = 0 under Kleinman symmetry. Miller scaling "
               "untested for RBBF.")

    def __init__(self):
        super().__init__()
        self._plane = 'arb'
        self._theta_rad = 'var'
        self._phi_rad = 'arb'

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
