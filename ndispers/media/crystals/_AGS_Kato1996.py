from ndispers._sym import sympy
from ndispers._baseclass import wl, phi, theta, T
from ndispers.groups import Uniax_42m
from ndispers.helper import vars2

class AGS(Uniax_42m):
    """
    AGS (AgGaS₂, silver thiogallate) crystal, Kato & Shirahata 1996 parameterisation

    - Point group : 4̄2m  (D2d)
    - Crystal system : Tetragonal
    - Dielectric principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Transparency range : about 0.47 to 13 µm

    Sellmeier equation
    ------------------
        n(wl)**2 = A_i + B_i / (wl**2 - C_i) + D_i * wl**2 + E_i * wl**4 + F_i * wl**6   for i = o, e
        (C in µm**2)

    Validity range
    ---------------
    0.54 to 12.9 µm

    Note
    ----
    The Sellmeier equation has no temperature term: T_degC is accepted for
    signature uniformity and ignored, and dndT returns 0.
    Coefficients as tabulated by refractiveindex.info (main/AgGaS2/nk/Kato-o, -e;
    accessed 2026-08-23) from the source; cross-checked against the Takaoka &
    Kato 1999 ordinary-ray set (AGS_Takaoka1999) to 1e-3 over 1-10 µm.

    Ref
    ---
    Sellmeier equation:
      Kato, K., & Shirahata, H. (1996). Nonlinear IR generation in AgGaS2. Japanese Journal of Applied Physics, 35(Part 1, No. 9A), 4645-4648. https://doi.org/10.1143/JJAP.35.4645
    Nonlinear optical coefficient:
      Roberts, D. A. (1992). Simplified characterization of uniaxial and biaxial nonlinear optical crystals: a plea for standardization of nomenclature and conventions. IEEE Journal of Quantum Electronics, 28(10), 2057-2074. https://doi.org/10.1109/3.159516
    """
    _wl_range = (0.54, 12.9)  # um, Sellmeier validity (see docstring)

    __slots__ = ["_A_o", "_B_o", "_C_o", "_D_o", "_E_o", "_F_o",
                 "_A_e", "_B_e", "_C_e", "_D_e", "_E_e", "_F_e"]

    _d_ref = {"d36": (17.5, 1.064, 1.064)}
    _d_note = ("Roberts 1992, Table V: 17.5 pm/V for 1.064 um SHG (relative to quartz), and "
               "11.2 pm/V for 10.6 um SHG, on his scale d36(KDP) = 0.39 pm/V. The 1.064 um "
               "value is held; Miller scaling from it gives 11.1 pm/V for 10.6 um SHG, against "
               "his 11.2 pm/V entry.")

    def __init__(self):
        super().__init__()
        self._plane = 'arb'
        self._theta_rad = 'var'
        self._phi_rad = 'arb'

        """ Constants of dispersion formula """
        # For ordinary ray
        self._A_o = 5.79419
        self._B_o = 0.23114
        self._C_o = 0.06882
        self._D_o = -2.4534e-3
        self._E_o = 3.1814e-7
        self._F_o = -9.7051e-9
        # For extraordinary ray
        self._A_e = 5.54120
        self._B_e = 0.22041
        self._C_e = 0.09824
        self._D_e = -2.5240e-3
        self._E_e = 3.6214e-7
        self._F_e = -8.3605e-9

    def n_o_expr(self):
        """ Sympy expression, dispersion formula for o-wave """
        return sympy.sqrt(self._A_o + self._B_o / (wl**2 - self._C_o) + self._D_o * wl**2 + self._E_o * wl**4 + self._F_o * wl**6)

    def n_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-wave """
        return sympy.sqrt(self._A_e + self._B_e / (wl**2 - self._C_e) + self._D_e * wl**2 + self._E_e * wl**4 + self._F_e * wl**6)

    def n_expr(self, pol):
        """ Sympy expression, dispersion formula of a general ray with an angle theta to optic axis. If theta = 0, this expression reduces to 'no_expre'. """
        if pol == 'o':
            return self.n_o_expr()
        elif pol == 'e':
            return self.n_e_expr() / sympy.sqrt( sympy.sin(theta)**2 + (self.n_e_expr()/self.n_o_expr())**2 * sympy.cos(theta)**2 )
        else:
            raise ValueError("pol = '%s' must be 'o' or 'e'" % pol)
