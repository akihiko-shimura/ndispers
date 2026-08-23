import sympy
from ndispers._baseclass import wl, phi, theta, T
from ndispers.groups import Uniax_42m
from ndispers.helper import vars2

class AGSe(Uniax_42m):
    """
    AGSe (AgGaSe₂, silver gallium selenide) crystal, Kato, Miyata & Petrov 2021 parameterisation

    - Point group : 4̄2m  (D2d)
    - Crystal system : Tetragonal
    - Dielectric principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Transparency range : about 0.71 to 18 µm

    Sellmeier equation
    ------------------
        n(wl)**2 = A_i + B_i / (wl**2 - C_i) + D_i / (wl**2 - E_i)   for i = o, e
        (C, E in µm**2)

    Validity range
    ---------------
    0.81 to 18 µm, at 20-22 degC

    Note
    ----
    The Sellmeier equation has no temperature term: T_degC is accepted for
    signature uniformity and ignored, and dndT returns 0.
    Coefficients as tabulated by refractiveindex.info (main/AgGaSe2/nk/Kato-o,
    -e; accessed 2026-08-23) from the source; they reproduce the earlier
    Harasaki & Kato set on refractiveindex.info to 2e-4 over 1-10 µm.

    Ref
    ---
    Sellmeier equation:
      Kato, K., Miyata, K., & Petrov, V. (2021). Refined Sellmeier equations for AgGaSe2 up to 18 µm. Applied Optics, 60(4), 805-808. https://doi.org/10.1364/ao.401828
    Nonlinear optical coefficient:
      Roberts, D. A. (1992). Simplified characterization of uniaxial and biaxial nonlinear optical crystals: a plea for standardization of nomenclature and conventions. IEEE Journal of Quantum Electronics, 28(10), 2057-2074. https://doi.org/10.1109/3.159516
    """
    __slots__ = ["_A_o", "_B_o", "_C_o", "_D_o", "_E_o",
                 "_A_e", "_B_e", "_C_e", "_D_e", "_E_e"]

    _d_ref = {"d36": (33.0, 10.6, 10.6)}
    _d_note = ("Roberts 1992, Table V: 33 pm/V for 10.6 um SHG (from a LiIO3-relative "
               "measurement), on his scale d36(KDP) = 0.39 pm/V. Miller scaling untested "
               "for AgGaSe2.")

    def __init__(self):
        super().__init__()
        self._plane = 'arb'
        self._theta_rad = 'var'
        self._phi_rad = 'arb'

        """ Constants of dispersion formula """
        # For ordinary ray
        self._A_o = 11.62264
        self._B_o = 0.43221
        self._C_o = 0.15597
        self._D_o = 18868.52
        self._E_o = 3953.34
        # For extraordinary ray
        self._A_e = 11.33168
        self._B_e = 0.45951
        self._C_e = 0.21284
        self._D_e = 17816.81
        self._E_e = 3828.78

    def n_o_expr(self):
        """ Sympy expression, dispersion formula for o-wave """
        return sympy.sqrt(self._A_o + self._B_o / (wl**2 - self._C_o) + self._D_o / (wl**2 - self._E_o))

    def n_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-wave """
        return sympy.sqrt(self._A_e + self._B_e / (wl**2 - self._C_e) + self._D_e / (wl**2 - self._E_e))

    def n_expr(self, pol):
        """ Sympy expression, dispersion formula of a general ray with an angle theta to optic axis. If theta = 0, this expression reduces to 'no_expre'. """
        if pol == 'o':
            return self.n_o_expr()
        elif pol == 'e':
            return self.n_e_expr() / sympy.sqrt( sympy.sin(theta)**2 + (self.n_e_expr()/self.n_o_expr())**2 * sympy.cos(theta)**2 )
        else:
            raise ValueError("pol = '%s' must be 'o' or 'e'" % pol)
