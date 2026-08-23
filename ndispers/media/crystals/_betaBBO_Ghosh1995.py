from ndispers._sym import sympy
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
        n(wl) = sqrt(A_i + B_i/(1 - C_i/wl**2) + D_i/(1 - E_i/wl**2)) + dn/dT * (T -20)

    Thermo-optic coefficient
    ------------------------
        dn/dT = (G_i * R_i + H_i * R_i**2) / 2*n_i  for i = o, e
        R_i = wl**2/(wl**2 - wl0_i**2), where wl0_i = 0.0652 for i = o or wl0 = 0.0730  for i = e

    Validity range
    --------------
    Fitted to refractive indices measured at 0.365 to 1.1 µm, supplemented
    by phase-matching-angle data in the UV and near infrared; the
    temperature dispersion is analyzed over about -40 to 320 degC.

    Ref
    ---
    Sellmeier equation and thermo-optic coefficients:
      Ghosh, G. (1995). Temperature dispersion of refractive indices in β-BaB2O4 and LiB3O5 crystals for nonlinear optical devices. Journal of Applied Physics, 78(11), 6752-6760. https://doi.org/10.1063/1.360499

    Nonlinear optical coefficients:
      Shoji, I., Kondo, T., & Ito, R. (1999). Absolute measurement of second-order nonlinear-optical coefficients of β-BaB2O4 for visible to ultraviolet second-harmonic wavelengths. JOSA B, 16(4), 620-624. https://doi.org/10.1364/josab.16.000620
    """
    __slots__ = ["_A_o", "_B_o", "_C_o", "_D_o", "_E_o",
                 "_A_e", "_B_e", "_C_e", "_D_e", "_E_e",
                 "_G_o", "_H_o",
                 "_G_e", "_H_e"]

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
        self._A_o = 1.7018379
        self._B_o = 1.0357554
        self._C_o = 0.018003440
        self._D_o = 1.2479989
        self._E_o = 91
        # For extraordinary ray
        self._A_e = 1.5920433
        self._B_e = 0.7816893
        self._C_e = 0.016067891
        self._D_e = 0.8403893
        self._E_e = 91
        # dn/dT
        self._G_o = -19.3007e-6
        self._G_e = -141.421e-6
        self._H_o = -34.9683e-6
        self._H_e = 110.8630e-6


    def _n_T20_o_expr(self):
        """ Sympy expression, dispersion formula for o-wave at 20degC """
        return sympy.sqrt(self._A_o + self._B_o / (1 - self._C_o/wl**2) + self._D_o/(1 - self._E_o/wl**2))

    def _n_T20_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-wave at 20degC """
        return sympy.sqrt(self._A_e + self._B_e / (1 - self._C_e/wl**2) + self._D_e/(1 - self._E_e/wl**2))

    def dndT_o_expr(self):
        return (self._G_o * self._R_o + self._H_o * self._R_o**2) / (2*self._n_T20_o_expr())

    def dndT_e_expr(self):
        return (self._G_e * self._R_e + self._H_e * self._R_e**2) / (2*self._n_T20_e_expr())

    @property
    def _R_o(self):
        # a sympy expression, built on demand so that constructing
        # the medium does not need sympy
        return wl**2/(wl**2 - 0.0652**2)

    @property
    def _R_e(self):
        # a sympy expression, built on demand so that constructing
        # the medium does not need sympy
        return wl**2/(wl**2 - 0.0730**2)

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
