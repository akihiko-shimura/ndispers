from ndispers._sym import sympy
from ndispers._baseclass import wl, phi, theta, T
from ndispers.groups import Uniax_42m
from ndispers.helper import vars2

class CLBO(Uniax_42m):
    """
    CLBO (CsLiB₆O₁₀, caesium lithium borate) crystal

    - Point group : 4̄2m  (D2d)
    - Crystal system : Tetragonal
    - Dielectric principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Transparency range : 0.18 to 2.75 µm

    Sellmeier equation
    ------------------
        n(wl) = sqrt(A_i + B_i/(wl**2 - C_i) - D_i * wl**2)  for i = o, e

    Thermo-optic coefficient
    -------------------------
        dn_o/dT = (At_o + Bt_o/wl)*1e-6
        dn_e/dT = (At_e + Bt_e/wl + Ct_e/wl**2 + Dt_e/wl**3)*1e-6

    The temperature correction dn/dT * (T - 20) vanishes at 20 degC.

    Validity range
    ---------------
    0.1914 to 2.09 µm
        dn/dT : 0.2128 to 1.3382 µm

    Ref
    ---
    Umemura, N., Yoshida, K., Kamimura, T., et al. (1999). New data on the phase-matching properties of CsLiB6O10. Advanced Solid State Lasers, PD15. https://doi.org/10.1364/assl.1999.pd15
    Nonlinear optical coefficient:
      Mori, Y., Kuroda, I., Nakajima, S., Sasaki, T., & Nakai, S. (1995). Nonlinear optical properties of cesium lithium borate. Japanese Journal of Applied Physics, 34(Part 2, No. 3A), L296-L298. https://doi.org/10.1143/JJAP.34.L296
    """
    __slots__ = ["_A_o", "_B_o", "_C_o", "_D_o",
                 "_A_e", "_B_e", "_C_e", "_D_e",
                 "_At_o", "_Bt_o",
                 "_At_e", "_Bt_e", "_Ct_e", "_Dt_e"]

    _d_ref = {"d36": (0.95, 1.064, 1.064)}
    _d_note = ("Mori et al. 1995, 1.064 um SHG, relative to d14(CBO). Miller scaling "
               "untested for CLBO.")

    def __init__(self):
        super().__init__()
        self._plane = 'arb'
        self._theta_rad = 'var'
        self._phi_rad = 'arb'

        """ Constants of dispersion formula """
        # For ordinary ray
        self._A_o = 2.2104
        self._B_o = 0.01018
        self._C_o = 0.01424
        self._D_o = 0.01258
        # For extraordinary ray
        self._A_e = 2.0588
        self._B_e = 0.00838
        self._C_e = 0.01363
        self._D_e = 0.00607
        # dn/dT
        self._At_o = -12.48 #1/K
        self._Bt_o = -0.328 #µm/K
        self._At_e = -8.36 #1/K
        self._Bt_e = 0.047 #µm/K
        self._Ct_e = 0.039 #µm^2/K
        self._Dt_e = 0.014 #µm^3/K


    def dndT_o_expr(self):
        return  (self._At_o + self._Bt_o/wl)*1e-6

    def dndT_e_expr(self):
        return  (self._At_e + self._Bt_e/wl + self._Ct_e/wl**2 + self._Dt_e/wl**3)*1e-6

    def n_o_expr(self):
        """ Sympy expression, dispersion formula for o-wave """
        return sympy.sqrt(self._A_o + self._B_o / (wl**2 - self._C_o) - self._D_o * wl**2) + self.dndT_o_expr() * (T - 20)

    def n_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-wave """
        return sympy.sqrt(self._A_e + self._B_e / (wl**2 - self._C_e) - self._D_e * wl**2) + self.dndT_e_expr() * (T - 20)

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
