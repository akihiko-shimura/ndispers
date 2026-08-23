from ndispers._sym import sympy

from ndispers._baseclass import T, phi, theta, wl
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
        n(wl) = sqrt(A_i + B_i/(wl**2 - C_i) - D_i * wl**2) + dn/dT * (T - 20)  for i = o, e

    Validity range
    ---------------
    0.22 to 1.06 µm

    Ref
    ---
    Sellmeier equation:
      Eimerl, D., Davis, L., Velsko, S., et al. (1987). Optical, mechanical, and thermal properties of barium borate. Journal of Applied Physics, 62(5), 1968-1983. https://doi.org/10.1063/1.339536

    Thermo-optic coefficients:
      Nikogosyan, D. N. (1991). Beta barium borate (BBO): a review of its properties and applications. Applied Physics A, 52(6), 359-368. https://doi.org/10.1007/bf00323647

    Nonlinear optical coefficients:
      Shoji, I., Kondo, T., & Ito, R. (1999). Absolute measurement of second-order nonlinear-optical coefficients of β-BaB2O4 for visible to ultraviolet second-harmonic wavelengths. JOSA B, 16(4), 620-624. https://doi.org/10.1364/josab.16.000620
    """
    __slots__ = ["_A_o", "_B_o", "_C_o", "_D_o",
                 "_A_e", "_B_e", "_C_e", "_D_e",
                 "_dndT_o", "_dndT_e"]

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
        self._A_o = 2.7405
        self._B_o = 0.0184
        self._C_o = 0.0179
        self._D_o = 0.0155
        # For extraordinary ray
        self._A_e = 2.3730
        self._B_e = 0.0128
        self._C_e = 0.0156
        self._D_e = 0.0044
        # dn/dT
        self._dndT_o = -16.6e-6 #/degC
        self._dndT_e = -9.3e-6 #/degC


    def n_o_expr(self):
        """ Sympy expression, dispersion formula for o-wave """
        return sympy.sqrt(self._A_o + self._B_o / (wl**2 - self._C_o) - self._D_o * wl**2) + self._dndT_o * (T - 20)

    def n_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-wave """
        return sympy.sqrt(self._A_e + self._B_e / (wl**2 - self._C_e) - self._D_e * wl**2) + self._dndT_e * (T - 20)

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
