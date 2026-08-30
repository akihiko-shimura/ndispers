from ndispers._sym import sympy

from ndispers._baseclass import T, phi, theta, wl
from ndispers.groups import Uniax_3m
from ndispers.helper import vars2

class SLT(Uniax_3m):
    """
    0.5% MgO-doped stoichiometric lithium tantalate (LiTaO₃) crystal

    - Point group : 3m  (C3v)
    - Crystal system : Trigonal
    - Dielectric principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Transparency range : 0.28 µm to 5.6 µm

    Sellmeier equation
    ------------------
        n(wl) = sqrt(a1_i + b1_i * f + (a2_i + b2_i * f)/(wl**2 - (a3_i + b3_i * f)**2) + (a4_i + b4_i * f)/(wl**2 - (a5_i + b5_i * f)**2) - a6_i * wl**2) for i=o,e
        f = (T - T0) * (T + T0 + 2 * 273.16) with T0 = 24.5 degC

    Validity range
    --------------
    Fitted to experimental data over 0.35 to 6 µm for the extraordinary
    wave and 0.375 to 3.75 µm for the ordinary wave, from room
    temperature up to 200 degC (Dolev 2009, abstract).

    Note
    ----
    d coefficients are of *congrunet* LT crystal, not of stoichiometric.

    Ref
    ---
    Sellmeier equation and ratio of d22 to d33:
      Dolev, I., Ganany-Padowicz, A., Gayer, O., et al. (2009). Linear and nonlinear optical properties of MgO:LiTaO3. Applied Physics B, 96(2-3), 423-432. https://doi.org/10.1007/s00340-009-3502-3

    Nonlinear optical coefficients d33, d31:
      Shoji, I., Kondo, T., Kitamoto, A., et al. (1997). Absolute scale of second-order nonlinear-optical coefficients. JOSA B, 14(9), 2268-2294. https://doi.org/10.1364/josab.14.002268
    """
    _wl_range = (0.35, 6)  # um, Sellmeier validity (see docstring)

    __slots__ = ["_a1_o", "_a2_o", "_a3_o", "_a4_o",  "_a5_o", "_a6_o",
                 "_a1_e", "_a2_e", "_a3_e", "_a4_e",  "_a5_e", "_a6_e",
                 "_b1_o", "_b2_o", "_b3_o", "_b4_o", "_b5_o",
                 "_b1_e", "_b2_e", "_b3_e", "_b4_e", "_b5_e"]

    _d_ref = {"d33": (13.8, 1.064, 1.064),
              "d31": (0.85, 1.064, 1.064),
              "d22": (0.12 * 13.8, 1.064, 1.064)}
    _d_note = ("d33 and d31: Shoji et al. 1997 (LiTaO3, 1.064 um SHG); d22 from the "
               "ratio d22/d33 = 0.12 of Dolev et al. 2009. Miller scaling untested "
               "for LiTaO3.")

    def __init__(self):
        super().__init__()
        self._plane = 'arb'
        self._theta_rad = 'var'
        self._phi_rad = 'arb'

        """ Constants of dispersion formula """
        # 0.5% MgO-doped SLT
        ## ordinary ray
        self._a1_o = 4.5082
        self._a2_o = 0.084888
        self._a3_o = 0.19552
        self._a4_o = 1.1570
        self._a5_o = 8.2517
        self._a6_o = 0.0237
        self._b1_o = 2.0704e-8
        self._b2_o = 1.4449e-8
        self._b3_o = 1.5978e-8
        self._b4_o = 4.7686e-6
        self._b5_o = 1.1127e-5
        ## extraordinary ray
        self._a1_e = 4.5615
        self._a2_e = 0.08488
        self._a3_e = 0.1927
        self._a4_e = 5.5832
        self._a5_e = 8.3067
        self._a6_e = 0.021696
        self._b1_e = 4.782e-7
        self._b2_e = 3.0913e-8
        self._b3_e = 2.7326e-8
        self._b4_e = 1.4837e-5
        self._b5_e = 1.3647e-7


    def n_o_expr(self):
        """ Sympy expression, dispersion formula for o-wave """
        return sympy.sqrt( self._a1_o + self._b1_o * self.f_expr() + \
            (self._a2_o + self._b2_o * self.f_expr()) / (wl**2 - (self._a3_o + self._b3_o * self.f_expr())**2) + \
                (self._a4_o + self._b4_o * self.f_expr()) / (wl**2 - (self._a5_o + self._b5_o * self.f_expr())**2) - self._a6_o * wl**2 )

    def n_e_expr(self):
        """ Sympy expression, dispersion formula for o-wave """
        return sympy.sqrt( self._a1_e + self._b1_e * self.f_expr() + \
            (self._a2_e + self._b2_e * self.f_expr()) / (wl**2 - (self._a3_e + self._b3_e * self.f_expr())**2) + \
                (self._a4_e + self._b4_e * self.f_expr()) / (wl**2 - (self._a5_e + self._b5_e * self.f_expr())**2) - self._a6_e * wl**2 )

    def f_expr(self):
        return (T - 24.5) * (T + 24.5 + 2 * 273.16)

    def n_expr(self, pol):
        """
        Sympy expression,
        dispersion formula of a general ray with an angle theta to optic axis. If theta = 0, this expression reduces to 'n_o_expr'.

        n(theta) = n_e / sqrt( sin(theta)**2 + (n_e/n_o)**2 * cos(theta)**2 )
        """
        if pol == 'o':
            return self.n_o_expr()
        elif pol == 'e':
            return self.n_e_expr() / sympy.sqrt( sympy.sin(theta)**2 + (self.n_e_expr()/self.n_o_expr())**2 * sympy.cos(theta)**2 )
        else:
            raise ValueError("pol = '%s' must be 'o' or 'e'" % pol)
