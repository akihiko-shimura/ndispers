from ndispers._sym import sympy
from ndispers._baseclass import wl, phi, theta, T
from ndispers.groups import Uniax_4mm
from ndispers.helper import vars2

class LB4(Uniax_4mm):
    """
    LB4 or LTB (Li₂B₄O₇, lithium tetraborate) crystal

    - Point group : 4mm  (C4v)
    - Crystal system : Tetragonal
    - Dielectric principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Transparency range : 0.16 to 3.6 µm

    Sellmeier equation
    ------------------
        n(wl) = sqrt(A_i + B_i/(wl**2 - C_i) - D_i * wl**2) for i = o,e

    Thermo-optic coefficient
    -------------------------
        dn/dT = At_i + Bt_i * wl + Ct_i * wl**2 + Dt_i * wl**3 + Et_i * (T - 25)  for i = o, e

    Validity range
    ---------------
    0.18 to 2.3 µm
    at 101325 Pa
        dn/dT : 0.436 to 0.644 µm, -40 to 60 degC

    Ref
    ---
    Sugawara, T., Komatsu, R., & Uda, S. (1998). Linear and nonlinear optical properties of lithium tetraborate. Solid State Communications, 107(5), 233-237. https://doi.org/10.1016/s0038-1098(98)00190-2
    Nonlinear optical coefficient:
      Komatsu, R., Sugawara, T., Sassa, K., et al. (1997). Growth and ultraviolet application of Li2B4O7 crystals: Generation of the fourth and fifth harmonics of Nd:Y3Al5O12 lasers. Applied Physics Letters, 70(26), 3492-3494. https://doi.org/10.1063/1.119210
      Petrov, V., Rotermund, F., Noack, F., Komatsu, R., Sugawara, T., & Uda, S. (1998). Vacuum ultraviolet application of Li2B4O7 crystals: Generation of 100 fs pulses down to 170 nm. Journal of Applied Physics, 84(11), 5887-5892. https://doi.org/10.1063/1.368904
    """
    __slots__ = ["_A_o", "_B_o", "_C_o", "_D_o",
                 "_A_e", "_B_e", "_C_e", "_D_e",
                 "_At_o", "_Bt_o", "_Ct_o", "_Dt_o", "_Et_o",
                 "_At_e", "_Bt_e", "_Ct_e", "_Dt_e", "_Et_e"]

    _d_ref = {"d31": (0.12, 1.064, 1.064)}
    _d_note = ("Komatsu et al. 1997, 1.064 um SHG, as tabulated by Petrov et al. 1998. "
               "d33 is not held (it enters only eee). Miller scaling untested for Li2B4O7.")

    def __init__(self):
        super().__init__()
        self._plane = 'arb'
        self._theta_rad = 'var'
        self._phi_rad = 'arb'

        """ Constants of dispersion formula """
        # For ordinary ray
        self._A_o = 2.564310
        self._B_o = 0.012337
        self._C_o = 0.114467**2
        self._D_o = 0.019075
        # For extraordinary ray
        self._A_e = 2.386510
        self._B_e = 0.010664
        self._C_e = 0.113483**2
        self._D_e = 0.012813
        # dn/dT coefficients
        # for o-wave
        self._At_o = 1.893e-5 #1/K
        self._Bt_o = -88.17e-6 #1/(µm*K)
        self._Ct_o = 1.497e-4 #1/(µm^2*K)
        self._Dt_o = -8.643e-5 #1/(µm^3*K)
        self._Et_o = -2.55e-8 #1/(K^2)
        # for e-wave
        self._At_e = 1.297e-5 #1/K
        self._Bt_e = -45.50e-6 #1/(µm*K)
        self._Ct_e = 0.714e-4 #1/(µm^2*K)
        self._Dt_e = -3.868e-5 #1/(µm^3*K)
        self._Et_e = -2.08e-8 #1/(K^2)


    def dndT_o_expr(self):
        return self._At_o + self._Bt_o * wl + self._Ct_o * wl**2 + self._Dt_o * wl**3 + self._Et_o * (T - 25)

    def dndT_e_expr(self):
        return self._At_e + self._Bt_e * wl + self._Ct_e * wl**2 + self._Dt_e * wl**3 + self._Et_e * (T - 25)

    def n_o_expr(self):
        """ Sympy expression, dispersion formula for o-wave """
        return sympy.sqrt(self._A_o + self._B_o / (wl**2 - self._C_o) - self._D_o * wl**2) + self.dndT_o_expr() * (T - 25)

    def n_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-wave """
        return sympy.sqrt(self._A_e + self._B_e / (wl**2 - self._C_e) - self._D_e * wl**2) + self.dndT_e_expr() * (T - 25)

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
