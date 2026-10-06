from ndispers._sym import sympy
from ndispers._baseclass import T, phi, theta, wl
from ndispers.groups import Uniax_3m
from ndispers.helper import vars2

class CLN(Uniax_3m):
    """
    Congruent lithium niobate (LiNbO₃), undoped, extraordinary ray, temperature-dependent

    - Point group : 3m  (C3v)
    - Crystal system : Trigonal
    - Dielectric principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Transparency range : 0.4 to 5 µm (absorption 0.08 /cm at 4 µm and 0.94 /cm at 5 µm)

    Sellmeier equation
    ------------------
        n_e(wl, T)**2 = a1 + b1 * f + (a2 + b2 * f) / (wl**2 - (a3 + b3 * f)**2) + (a4 + b4 * f) / (wl**2 - a5**2) - a6 * wl**2
        f = (T - 24.5) * (T + 570.82)
        (wl in µm, T in degC)

    The temperature dependence is part of the equation: dndT is its
    derivative, and f vanishes at 24.5 degC.

    Validity range
    ---------------
    0.4 to 5 µm, at room temperature to 250 degC

    Note
    ----
    Extraordinary ray only (pol='e' is the default of every method and 'o'
    raises): the equation of choice for quasi-phase-matching in undoped
    periodically poled LiNbO₃. It refits the data behind CLN_Edwards1984
    together with the tuning of a 1.064 µm-pumped PPLN optical parametric
    oscillator, adding an infrared oscillator at a5 so that the mid-infrared
    idler is predicted correctly. For the ordinary ray use CLN_Edwards1984
    (with temperature) or CLN_Zelmon1997 (21 degC).

    Ref
    ---
    Sellmeier equation (Eq. 4, Table 2):
      Jundt, D. H. (1997). Temperature-dependent Sellmeier equation for the index of refraction, n_e, in congruent lithium niobate. Optics Letters, 22(20), 1553-1555. https://doi.org/10.1364/ol.22.001553
    Nonlinear optical coefficient:
      Shoji, I., Kondo, T., Kitamoto, A., Shirane, M., & Ito, R. (1997). Absolute scale of second-order nonlinear-optical coefficients. JOSA B, 14(9), 2268-2294. https://doi.org/10.1364/josab.14.002268
    """
    _wl_range = (0.4, 5)  # um, Sellmeier validity (see docstring)

    __slots__ = ["_a1_e", "_a2_e", "_a3_e", "_a4_e", "_a5_e", "_a6_e",
                 "_b1_e", "_b2_e", "_b3_e", "_b4_e"]

    # Only d33 is usable: this class has no o-ray Sellmeier equation, so the
    # o-ray susceptibilities that Miller scaling of d22 and d31 (and any
    # interaction with an o-wave) would need are not available. eee
    # (quasi-phase-matching) works.
    _d_ref = {"d33": (25.2, 1.064, 1.064)}
    _d_note = ("Shoji et al. 1997 (congruent LiNbO3, 1.064 um SHG, absolute). Only d33 is "
               "held, so only eee (QPM) can be evaluated: this class has no o-ray Sellmeier "
               "equation. Alford & Smith 2001 find Miller scaling good for LiNbO3.")

    _default_pol = 'e'   # the source gives the extraordinary index only

    def __init__(self):
        super().__init__()
        self._plane = 'arb'
        self._theta_rad = 'var'
        self._phi_rad = 'arb'

        """ Constants of dispersion formula, Jundt 1997 Table 2 """
        self._a1_e = 5.35583
        self._a2_e = 0.100473
        self._a3_e = 0.20692
        self._a4_e = 100
        self._a5_e = 11.34927
        self._a6_e = 1.5334e-2
        self._b1_e = 4.629e-7
        self._b2_e = 3.862e-8
        self._b3_e = -0.89e-8
        self._b4_e = 2.657e-5

    def n_e_expr(self):
        """ Sympy expression, dispersion formula for e-wave """
        return sympy.sqrt( self._a1_e + self._b1_e * self.f_expr() + \
            (self._a2_e + self._b2_e * self.f_expr()) / (wl**2 - (self._a3_e + self._b3_e * self.f_expr())**2) + \
                (self._a4_e + self._b4_e * self.f_expr()) / (wl**2 - self._a5_e**2) - self._a6_e * wl**2 )

    def f_expr(self):
        return (T - 24.5) * (T + 570.82)

    def n_expr(self, pol):
        """
        Sympy expression,
        dispersion formula,
        only for e-wave

        """
        if pol == 'e':
            return self.n_e_expr()
        else:
            raise ValueError("pol = '%s' must be 'e'. Sellmeier equation for pol='o' is not implemented for this module." % pol)
