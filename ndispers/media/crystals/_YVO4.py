from ndispers._sym import sympy
from ndispers._baseclass import Medium, wl, phi, theta, T
from ndispers.helper import vars2

class YVO4(Medium):
    """
    YVO₄ (yttrium orthovanadate) crystal, undoped

    - Point group : 4/m  (C4h); centrosymmetric, so no second-order nonlinearity
    - Crystal system : Tetragonal (zircon structure)
    - Dielectric principal axis, z // c-axis (x, y-axes are arbitrary)
    - Positive uniaxial, with optic axis parallel to z-axis
    - Transparency range : about 0.4 to 5 µm

    Sellmeier equation
    ------------------
        n(wl)**2 = A_i + B_i / (wl**2 - C_i) - D_i * wl**2   for i = o, e
        (C in µm**2)

    Validity range
    ---------------
    0.48 to 1.34 µm, at 20 degC

    Note
    ----
    The Sellmeier equation has no temperature term: T_degC is accepted for
    signature uniformity and ignored, and dndT returns 0 (the source also
    reports thermo-optic coefficients, not implemented). The large positive
    birefringence (n_e - n_o = 0.21 at 1.064 µm) is what makes YVO₄ the usual
    material for walk-off polarizers and displacers. Coefficients as tabulated
    by refractiveindex.info (main/YVO4/nk/Shi-o-20C, Shi-e-20C; accessed
    2026-08-23) from the source; the earlier Birnbaum & DeShazer set on the same
    site agrees to 3e-3 from 0.6 to 1.2 µm (6e-3 at 0.5 µm).

    Ref
    ---
    Sellmeier equation:
      Shi, H. S., Zhang, G., & Shen, H. Y. (2001). Measurement of principal refractive indices and the thermal refractive index coefficients of yttrium vanadate. Journal of Synthetic Crystals, 30(1), 85-88 (not registered with Crossref).
    """
    _wl_range = (0.48, 1.34)  # um, Sellmeier validity (see docstring)

    __slots__ = ["_A_o", "_B_o", "_C_o", "_D_o",
                 "_A_e", "_B_e", "_C_e", "_D_e"]

    def __init__(self):
        super().__init__()
        self._plane = 'arb'
        self._theta_rad = 'var'
        self._phi_rad = 'arb'

        """ Constants of dispersion formula """
        # For ordinary ray
        self._A_o = 3.778790
        self._B_o = 0.070479
        self._C_o = 0.045731
        self._D_o = 0.009701
        # For extraordinary ray
        self._A_e = 4.607200
        self._B_e = 0.108087
        self._C_e = 0.052495
        self._D_e = 0.014305

    def n_o_expr(self):
        """ Sympy expression, dispersion formula for o-wave """
        return sympy.sqrt(self._A_o + self._B_o / (wl**2 - self._C_o) - self._D_o * wl**2)

    def n_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-wave """
        return sympy.sqrt(self._A_e + self._B_e / (wl**2 - self._C_e) - self._D_e * wl**2)

    def n_expr(self, pol):
        """ Sympy expression, dispersion formula of a general ray with an angle theta to optic axis. If theta = 0, this expression reduces to 'no_expre'. """
        if pol == 'o':
            return self.n_o_expr()
        elif pol == 'e':
            return self.n_e_expr() / sympy.sqrt( sympy.sin(theta)**2 + (self.n_e_expr()/self.n_o_expr())**2 * sympy.cos(theta)**2 )
        else:
            raise ValueError("pol = '%s' must be 'o' or 'e'" % pol)
