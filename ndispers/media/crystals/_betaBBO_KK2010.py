import sympy
from ndispers._baseclass import wl, phi, theta, T
from ndispers.groups import Uniax_3m
from ndispers.helper import vars2

class BetaBBO(Uniax_3m):
    """
    β-BBO (β-BaB₂O₄, barium borate) crystal

    - Point group : 3m  (C3v)
    - Crystal system : Trigonal
    - Dielectric principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Transparency range : 0.19 µm to 2.6 µm

    Sellmeier equation
    ------------------
        n(wl) = sqrt(A_i + B_i/(wl**2 - C_i) + D_i/(wl**2 - E_i)) + dn/dT * (T -20)  i = o, e

    Thermo-optic coefficient
    ------------------------
        dn/dT = F_i/wl**3 + G_i/wl**2 + H_i/wl + I_i  for i = o, e

    Validity range
    --------------
    o: 0.2048 to 3.22 µm
    e: 0.1916 to 0.2048 µm
        dn/dT: 0.195 to 1.618 µm

    Ref
    ---
    Sellmeier equation and thermo-optic coefficients:
      Kato, K., Umemura, N., & Mikami, T. (2010). Sellmeier and thermo-optic dispersion formulas for β-BaB2O4 (revisited). Nonlinear Frequency Generation and Conversion: Materials, Devices, and Applications IX, Proc. SPIE, 7582, 75820H.  (not registered with Crossref)

    Nonlinear optical coefficients:
      Shoji, I., Kondo, T., & Ito, R. (1999). Absolute measurement of second-order nonlinear-optical coefficients of β-BaB2O4 for visible to ultraviolet second-harmonic wavelengths. JOSA B, 16(4), 620-624. https://doi.org/10.1364/josab.16.000620
    """
    __slots__ = ["_A_o", "_B_o", "_C_o", "_D_o", "_E_o",
                 "_A_e", "_B_e", "_C_e", "_D_e", "_E_e",
                 "_F_o", "_G_o", "_H_o", "_I_o",
                 "_F_e", "_G_e", "_H_e", "_I_e"]

    # Second-order nonlinear coefficients, pm/V, at the (wl1, wl2) of their
    # measurement; scaled to other wavelengths by Miller's rule (see groups).
    _d_ref = {"d22": (2.2, 1.064, 1.064),
              "d31": (0.04, 1.064, 1.064),
              "d33": (0.04, 1.064, 1.064)}
    _d_note = ("Shoji et al. 1999: d22 absolute by the wedge technique, d31 and d33 "
               "relative to d22, all at 1.064 um SHG; d22 and d31 of the same sign "
               "(Alford & Smith 2001 take d_yyy/d_zxx = +55). Miller scaling of d22 "
               "is supported by Alford & Smith 2001 over 532-1319 nm SHG.")

    def __init__(self):
        super().__init__()
        self._plane = 'arb'
        self._theta_rad = 'var'
        self._phi_rad = 'arb'

        """ Constants of dispersion formula """
        # For ordinary ray
        self._A_o = 3.63357
        self._B_o = 0.01878
        self._C_o = 0.01822
        self._D_o = 60.9129
        self._E_o = 67.8505
        # For extraordinary ray (best fit from 0.2048 to 3.22 µm)
        self._A_e = 3.33469
        self._B_e = 0.01237
        self._C_e = 0.01647
        self._D_e = 79.0672
        self._E_e = 82.2919
        # For extraordinary ray (best fit from 0.1916 to 0.2048 µm)
        # self._A_e = 3.38630
        # self._B_e = 0.00921
        # self._C_e = 0.02073
        # self._D_e = 79.0672
        # self._E_e = 82.2919
        # dn/dT (best fit from 0.195 to 1.618 µm)
        self._F_o = -0.0137e-5
        self._F_e = 0.0413e-5
        self._G_o = 0.0607e-5
        self._G_e = -0.2119e-5
        self._H_o = -0.1334e-5
        self._H_e = 0.4408e-5
        self._I_o = -1.5287e-5
        self._I_e = -1.2749e-5


    def _n_T20_o_expr(self):
        """ Sympy expression, dispersion formula for o-wave at 20degC """
        return sympy.sqrt(self._A_o + self._B_o/(wl**2 - self._C_o) + self._D_o/(wl**2 - self._E_o))

    def _n_T20_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-wave at 20degC """
        return sympy.sqrt(self._A_e + self._B_e/(wl**2 - self._C_e) + self._D_e/(wl**2 - self._E_e))

    def dndT_o_expr(self):
        return self._F_o/wl**3 + self._G_o/wl**2 + self._H_o/wl + self._I_o

    def dndT_e_expr(self):
        return self._F_e/wl**3 + self._G_e/wl**2 + self._H_e/wl + self._I_e

    def n_o_expr(self):
        return self._n_T20_o_expr() + self.dndT_o_expr() * (T - 20)

    def n_e_expr(self):
        return self._n_T20_e_expr() + self.dndT_e_expr() * (T - 20)

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
