from ndispers._sym import sympy
from ndispers._baseclass import wl, phi, theta, T
from ndispers.groups import Uniax_3m
from ndispers.helper import vars2

class CLN(Uniax_3m):
    """
    Congruent lithium niobate (LiNbO₃), undoped, both rays, at 21 degC

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
    Minimum-deviation measurement on a crystal grown from a congruent melt
    (48.38 mol% Li₂O); the fit reproduces the measured indices within 2e-4.
    The third oscillator is in the infrared, which is what carries the fit to
    5 µm. Unlike Table 2 of the same paper (the 5% MgO-doped crystal, see
    MgOLN_Zelmon1997), the column labels of Table 1 are right as printed:
    n_o = 2.232, n_e = 2.156 at 1.064 µm.

    The Sellmeier equation has no temperature term - it was fitted at
    21 degC - so T_degC is accepted and ignored, and dndT returns 0. With
    temperature: CLN_Edwards1984 (both rays, to 3.4 µm) and CLN_Jundt1997
    (extraordinary ray, to 5 µm).

    Ref
    ---
    Sellmeier equation (Table 1):
      Zelmon, D. E., Small, D. L., & Jundt, D. (1997). Infrared corrected Sellmeier coefficients for congruently grown lithium niobate and 5 mol. % magnesium oxide-doped lithium niobate. JOSA B, 14(12), 3319-3322. https://doi.org/10.1364/josab.14.003319
    Nonlinear optical coefficients:
      Shoji, I., Kondo, T., Kitamoto, A., Shirane, M., & Ito, R. (1997). Absolute scale of second-order nonlinear-optical coefficients. JOSA B, 14(9), 2268-2294. https://doi.org/10.1364/josab.14.002268
      Roberts, D. A. (1992). Simplified characterization of uniaxial and biaxial nonlinear optical crystals: a plea for standardization of nomenclature and conventions. IEEE Journal of Quantum Electronics, 28(10), 2057-2074. https://doi.org/10.1109/3.159516
    """
    _wl_range = (0.4, 5.0)  # um, Sellmeier validity (see docstring)

    __slots__ = ["_A_o", "_B_o", "_C_o", "_D_o", "_E_o", "_F_o",
                 "_A_e", "_B_e", "_C_e", "_D_e", "_E_e", "_F_e"]

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

        """ Constants of dispersion formula, Zelmon et al. 1997 Table 1 """
        # For ordinary ray
        self._A_o = 2.6734
        self._B_o = 0.01764
        self._C_o = 1.2290
        self._D_o = 0.05914
        self._E_o = 12.614
        self._F_o = 474.6
        # For extraordinary ray
        self._A_e = 2.9804
        self._B_e = 0.02047
        self._C_e = 0.5981
        self._D_e = 0.0666
        self._E_e = 8.9543
        self._F_e = 416.08

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
