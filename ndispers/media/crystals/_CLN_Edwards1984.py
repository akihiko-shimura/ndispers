from ndispers._sym import sympy
from ndispers._baseclass import wl, phi, theta, T
from ndispers.groups import Uniax_3m
from ndispers.helper import vars2

class CLN(Uniax_3m):
    """
    Congruent lithium niobate (LiNbO₃), undoped, both rays, temperature-dependent

    - Point group : 3m  (C3v)
    - Crystal system : Trigonal
    - Dielectric principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Transparency range : 0.4 to 4.5 µm

    Sellmeier equation
    ------------------
        n(wl, T)**2 = A1_i + (A2_i + B1_i * F) / (wl**2 - (A3_i + B2_i * F)**2) + B3_i * F - A4_i * wl**2   for i = o, e
        F = (T - 24.5) * (T + 24.5 + 546)
        (wl in µm, T in degC)

    The temperature dependence is part of the equation: dndT is its
    derivative, and F vanishes at 24.5 degC.

    Validity range
    ---------------
    0.4 to 3.4 µm, at 0 to 500 degC

    Note
    ----
    Fitted to the room-temperature indices of Nelson & Mikulyak 1974 (0.4 to
    3.1 µm, 24.5 degC) and to the temperature dependence measured by Smith et
    al. 1976 at 0.633 and 3.39 µm. The source expects it to serve over 0.4 to
    4.5 µm, but Jundt 1997 found it inadequate for phase matching beyond
    3.3 µm, more so at high temperature: the infrared absorption is a single
    term in wl**2. The range above is that of the fitted data. For the
    extraordinary ray in the mid-infrared use CLN_Jundt1997; for both rays at
    room temperature to 5 µm, CLN_Zelmon1997.

    This is the only undoped congruent set with a temperature-dependent
    ordinary index, which birefringent phase matching needs. The source's
    Table II (noncritical type-I difference-frequency mixing of 488 nm with a
    dye laser, 2.2 to 3.2 µm, 176 to 365 degC) is reproduced within 1 degC.

    Ref
    ---
    Sellmeier equation (Eq. 1, Table I):
      Edwards, G. J., & Lawrence, M. (1984). A temperature-dependent dispersion equation for congruently grown lithium niobate. Optical and Quantum Electronics, 16(4), 373-375. https://doi.org/10.1007/bf00620081
    Nonlinear optical coefficients:
      Shoji, I., Kondo, T., Kitamoto, A., Shirane, M., & Ito, R. (1997). Absolute scale of second-order nonlinear-optical coefficients. JOSA B, 14(9), 2268-2294. https://doi.org/10.1364/josab.14.002268
      Roberts, D. A. (1992). Simplified characterization of uniaxial and biaxial nonlinear optical crystals: a plea for standardization of nomenclature and conventions. IEEE Journal of Quantum Electronics, 28(10), 2057-2074. https://doi.org/10.1109/3.159516
    """
    _wl_range = (0.4, 3.4)  # um, Sellmeier validity (see docstring)

    __slots__ = ["_A1_o", "_A2_o", "_A3_o", "_A4_o", "_B1_o", "_B2_o", "_B3_o",
                 "_A1_e", "_A2_e", "_A3_e", "_A4_e", "_B1_e", "_B2_e", "_B3_e"]

    _d_ref = {"d33": (25.2, 1.064, 1.064),
              "d31": (4.6, 1.064, 1.064),
              "d22": (-2.1, 1.064, 1.064)}
    _d_note = ("d33 and d31: Shoji et al. 1997 (congruent LiNbO3, 1.064 um SHG, absolute; "
               "Table 6 also lists 19.5 and 3.2 at 1.313 um, 25.7 and 4.8 at 0.852 um). "
               "d22: Roberts 1992 (congruent undoped LiNbO3); its sign is opposite to "
               "d31's (Alford & Smith 2001 use d_yyy/d_zxx = -0.49). Alford & Smith 2001 "
               "find Miller scaling good for LiNbO3.")

    def __init__(self):
        super().__init__()
        self._plane = 'arb'
        self._theta_rad = 'var'
        self._phi_rad = 'arb'

        """ Constants of dispersion formula, Edwards & Lawrence 1984 Table I """
        # For ordinary ray
        self._A1_o = 4.9048
        self._A2_o = 0.11775
        self._A3_o = 0.21802
        self._A4_o = 0.027153
        self._B1_o = 2.2314e-8
        self._B2_o = -2.9671e-8
        self._B3_o = 2.1429e-8
        # For extraordinary ray
        self._A1_e = 4.5820
        self._A2_e = 0.09921
        self._A3_e = 0.21090
        self._A4_e = 0.021940
        self._B1_e = 5.2716e-8
        self._B2_e = -4.9143e-8
        self._B3_e = 2.2971e-7

    def F_expr(self):
        return (T - 24.5) * (T + 24.5 + 546)

    def n_o_expr(self):
        """ Sympy expression, dispersion formula for o-wave """
        return sympy.sqrt(self._A1_o + (self._A2_o + self._B1_o * self.F_expr()) / (wl**2 - (self._A3_o + self._B2_o * self.F_expr())**2) + self._B3_o * self.F_expr() - self._A4_o * wl**2)

    def n_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-wave """
        return sympy.sqrt(self._A1_e + (self._A2_e + self._B1_e * self.F_expr()) / (wl**2 - (self._A3_e + self._B2_e * self.F_expr())**2) + self._B3_e * self.F_expr() - self._A4_e * wl**2)

    def n_expr(self, pol):
        """ Sympy expression, dispersion formula of a general ray with an angle theta to optic axis. If theta = 0, this expression reduces to 'no_expre'. """
        if pol == 'o':
            return self.n_o_expr()
        elif pol == 'e':
            return self.n_e_expr() / sympy.sqrt( sympy.sin(theta)**2 + (self.n_e_expr()/self.n_o_expr())**2 * sympy.cos(theta)**2 )
        else:
            raise ValueError("pol = '%s' must be 'o' or 'e'" % pol)
