import sympy
from ndispers._baseclass import wl, phi, theta, T
from ndispers.groups import Uniax_3m
from ndispers.helper import vars2

class MgOLN(Uniax_3m):
    """
    5 mol% MgO-doped congruent lithium niobate (LiNbO₃) crystal, both rays

    - Point group : 3m  (C3v)
    - Crystal system : Trigonal
    - Dielectric principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Transparency range : 0.4 to 5.0 µm (measured range of the source)

    Sellmeier equation
    ------------------
        n(wl)**2 = 1 + A_i * wl**2 / (wl**2 - B_i) + C_i * wl**2 / (wl**2 - D_i) + E_i * wl**2 / (wl**2 - F_i)   for i = o, e
        (B, D, F in µm**2)

    Validity range
    ---------------
    0.4 to 5.0 µm, at 21 degC

    Note
    ----
    The column labels of Table 2 in the source are interchanged in print:
    read literally they give n_e > n_o, while LiNbO₃ is negative uniaxial and
    the paper's Fig. 2 shows n_o above n_e. The coefficients are assigned
    here so that n_o > n_e (n_o = 2.229, n_e = 2.147 at 1.064 µm).

    This is the class to use for birefringent phase matching in LiNbO₃: it
    has both the ordinary and the extraordinary ray, which SLN (e-ray only)
    lacks, so ooe/oee d_eff can be evaluated. The Sellmeier equation has no
    temperature term - it was fitted at 21 degC - so T_degC is accepted and
    ignored, and dndT returns 0; for the temperature dependence of the
    extraordinary index of MgO-doped LiNbO₃ use Gayer et al. 2008 (the SLN
    class implements their 1% MgO stoichiometric set).

    Ref
    ---
    Sellmeier equation:
      Zelmon, D. E., Small, D. L., & Jundt, D. (1997). Infrared corrected Sellmeier coefficients for congruently grown lithium niobate and 5 mol. % magnesium oxide-doped lithium niobate. JOSA B, 14(12), 3319-3322. https://doi.org/10.1364/josab.14.003319
    Nonlinear optical coefficients:
      Shoji, I., Kondo, T., Kitamoto, A., Shirane, M., & Ito, R. (1997). Absolute scale of second-order nonlinear-optical coefficients. JOSA B, 14(9), 2268-2294. https://doi.org/10.1364/josab.14.002268
      Roberts, D. A. (1992). Simplified characterization of uniaxial and biaxial nonlinear optical crystals: a plea for standardization of nomenclature and conventions. IEEE Journal of Quantum Electronics, 28(10), 2057-2074. https://doi.org/10.1109/3.159516
    """
    __slots__ = ["_A_o", "_B_o", "_C_o", "_D_o", "_E_o", "_F_o",
                 "_A_e", "_B_e", "_C_e", "_D_e", "_E_e", "_F_e"]

    _d_ref = {"d33": (25.0, 1.064, 1.064),
              "d31": (4.4, 1.064, 1.064),
              "d22": (-2.1, 1.064, 1.064)}
    _d_note = ("d33 and d31: Shoji et al. 1997 (5%MgO:LiNbO3, 1.064 um SHG, absolute). "
               "d22: Roberts 1992 (congruent undoped LiNbO3); its sign is opposite to "
               "d31's (Alford & Smith 2001 use d_yyy/d_zxx = -0.49). Alford & Smith 2001 "
               "find Miller scaling good for LiNbO3.")

    def __init__(self):
        super().__init__()
        self._plane = 'arb'
        self._theta_rad = 'var'
        self._phi_rad = 'arb'

        """ Constants of dispersion formula, Zelmon et al. 1997 Table 2 """
        # The printed Table 2 heads these two columns "n_e" and "n_o" in that
        # order, which would make n_e > n_o; LiNbO3 is negative uniaxial, the
        # paper's own Fig. 2 shows n_o above n_e, and the undoped crystal's
        # Table 1 has the larger A in the n_o column. The labels are taken as
        # interchanged in print and assigned here by the physics.
        # For ordinary ray (printed under "n_e")
        self._A_o = 2.4272
        self._B_o = 0.01478
        self._C_o = 1.4617
        self._D_o = 0.05612
        self._E_o = 9.6536
        self._F_o = 371.216
        # For extraordinary ray (printed under "n_o")
        self._A_e = 2.2454
        self._B_e = 0.01242
        self._C_e = 1.3005
        self._D_e = 0.05313
        self._E_e = 6.8972
        self._F_e = 331.33

    def n_o_expr(self):
        """ Sympy expression, dispersion formula for o-wave """
        return sympy.sqrt(1 + self._A_o * wl**2 / (wl**2 - self._B_o) + self._C_o * wl**2 / (wl**2 - self._D_o) + self._E_o * wl**2 / (wl**2 - self._F_o))

    def n_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-wave """
        return sympy.sqrt(1 + self._A_e * wl**2 / (wl**2 - self._B_e) + self._C_e * wl**2 / (wl**2 - self._D_e) + self._E_e * wl**2 / (wl**2 - self._F_e))

    def n_expr(self, pol):
        """ Sympy expression, dispersion formula of a general ray with an angle theta to optic axis. If theta = 0, this expression reduces to 'no_expre'. """
        if pol == 'o':
            return self.n_o_expr()
        elif pol == 'e':
            return self.n_e_expr() / sympy.sqrt( sympy.sin(theta)**2 + (self.n_e_expr()/self.n_o_expr())**2 * sympy.cos(theta)**2 )
        else:
            raise ValueError("pol = '%s' must be 'o' or 'e'" % pol)
