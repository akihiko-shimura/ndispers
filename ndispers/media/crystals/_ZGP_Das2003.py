from ndispers._sym import sympy
from ndispers._baseclass import wl, phi, theta, T
from ndispers.groups import Uniax_42m
from ndispers.helper import vars2

class ZGP(Uniax_42m):
    """
    ZGP (ZnGeP₂, zinc germanium phosphide) crystal, Das2003 parameterisation

    - Point group : 4̄2m  (D2d)
    - Crystal system : Tetragonal
    - Dielectric principal axis, z // c-axis (x, y-axes are arbitrary)
    - Positive uniaxial, with optic axis parallel to z-axis
    - Transparency range : about 0.74 to 12 µm (strong absorption below 2 µm)

    Sellmeier equation
    ------------------
        n(wl)**2 = A_i + B_i * wl**2 / (wl**2 - C_i) + D_i * wl**2 / (wl**2 - E_i)   for i = o, e
        (C, E in µm**2)

    Validity range
    ---------------
    1.32 to 11 µm

    Note
    ----
    Its type-I limit, 10.74 um, is close to the 10.78 um of Kato 1997, and it does phase-match SHG of the 10.6 um CO2 line (81.4 deg) where ZGP_Zelmon2001 falls 5e-4 in index short. Published ZGP Sellmeier sets differ by 10-20 degrees near the range boundaries.

    The Sellmeier equation has no temperature term: T_degC is accepted for
    signature uniformity and ignored, and dndT returns 0.
    Coefficients as tabulated by refractiveindex.info (main/ZnGeP2/nk/Das-o, -e; accessed 2026-08-23) from the source; cross-checked against the Zelmon 2001 set (ZGP_Zelmon2001) to 5e-3 over 2.5-8 µm.

    Ref
    ---
    Sellmeier equation:
      Das, S., Bhar, G. C., Gangopadhyay, S., & Ghosh, C. (2003). Linear and nonlinear optical properties of ZnGeP2 crystal for infrared laser device applications: revisited. Applied Optics, 42(21), 4335-4340. https://doi.org/10.1364/ao.42.004335
    Nonlinear optical coefficient:
      Roberts, D. A. (1992). Simplified characterization of uniaxial and biaxial nonlinear optical crystals: a plea for standardization of nomenclature and conventions. IEEE Journal of Quantum Electronics, 28(10), 2057-2074. https://doi.org/10.1109/3.159516
    """
    _wl_range = (1.32, 11)  # um, Sellmeier validity (see docstring)

    __slots__ = ["_A_o", "_B_o", "_C_o", "_D_o", "_E_o",
                 "_A_e", "_B_e", "_C_e", "_D_e", "_E_e"]

    _d_ref = {"d36": (69.0, 10.6, 10.6)}
    _d_note = ("Roberts 1992, Table VI: 69 pm/V for 10.6 um SHG, on his scale d36(KDP) = "
               "0.39 pm/V. Miller scaling untested for ZnGeP2.")

    def __init__(self):
        super().__init__()
        self._plane = 'arb'
        self._theta_rad = 'var'
        self._phi_rad = 'arb'

        """ Constants of dispersion formula """
        # For ordinary ray
        self._A_o = 1 + 4.67491   # RII 'formula 2' lists n**2 - 1 = c1 + ...
        self._B_o = 4.077926
        self._C_o = 0.159328
        self._D_o = 1.896005
        self._E_o = 900.0
        # For extraordinary ray
        self._A_e = 1 + 2.65014
        self._B_e = 6.310153
        self._C_e = 0.125099
        self._D_e = 1.731381
        self._E_e = 900.0

    def n_o_expr(self):
        """ Sympy expression, dispersion formula for o-wave """
        return sympy.sqrt(self._A_o + self._B_o * wl**2 / (wl**2 - self._C_o) + self._D_o * wl**2 / (wl**2 - self._E_o))

    def n_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-wave """
        return sympy.sqrt(self._A_e + self._B_e * wl**2 / (wl**2 - self._C_e) + self._D_e * wl**2 / (wl**2 - self._E_e))

    def n_expr(self, pol):
        """ Sympy expression, dispersion formula of a general ray with an angle theta to optic axis. If theta = 0, this expression reduces to 'no_expre'. """
        if pol == 'o':
            return self.n_o_expr()
        elif pol == 'e':
            return self.n_e_expr() / sympy.sqrt( sympy.sin(theta)**2 + (self.n_e_expr()/self.n_o_expr())**2 * sympy.cos(theta)**2 )
        else:
            raise ValueError("pol = '%s' must be 'o' or 'e'" % pol)
