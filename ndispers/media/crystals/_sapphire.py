from ndispers._sym import sympy
from ndispers._baseclass import Medium, wl, phi, theta, T
from ndispers.helper import vars2

class Sapphire(Medium):
    """
    Sapphire (α-Al₂O₃, corundum) crystal

    - Point group : 3̄m  (D3d); centrosymmetric, so no second-order nonlinearity
    - Crystal system : Trigonal
    - Dielectric principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Transparency range : 0.2 to 5.5 µm (the range of the dispersion fit)

    Sellmeier equation
    ------------------
        n(wl)**2 = 1 + A1_i * wl**2 / (wl**2 - B1_i**2) + A2_i * wl**2 / (wl**2 - B2_i**2) + A3_i * wl**2 / (wl**2 - B3_i**2)   for i = o, e

    Validity range
    ---------------
    0.2 to 5.5 µm

    Note
    ----
    The Sellmeier equation has no temperature term: T_degC is accepted for
    signature uniformity and ignored, and dndT returns 0 (sapphire's dn/dT is
    positive, of order 1e-5 /K; use a handbook value if an estimate is
    needed).

    Ref
    ---
    Sellmeier equation (Malitson & Dodge's fit, as tabulated in):
      Tropf, W. J., Thomas, M. E., & Harris, T. J. (1995). Properties of crystals and glasses. In M. Bass (Ed.), Handbook of Optics, Vol. II, Chapter 33, Table 22. McGraw-Hill.
    Underlying measurements:
      Malitson, I. H. (1962). Refraction and dispersion of synthetic sapphire. JOSA, 52(12), 1377-1379. https://doi.org/10.1364/josa.52.001377
      Malitson, I. H., & Dodge, M. J. (1972). Refractive index and birefringence of synthetic sapphire. JOSA, 62, 1405A (abstract; not registered with Crossref).
    """
    _wl_range = (0.2, 5.5)  # um, Sellmeier validity (see docstring)

    __slots__ = ["_A1_o", "_B1_o", "_A2_o", "_B2_o", "_A3_o", "_B3_o",
                 "_A1_e", "_B1_e", "_A2_e", "_B2_e", "_A3_e", "_B3_e"]

    def __init__(self):
        super().__init__()
        self._plane = 'arb'
        self._theta_rad = 'var'
        self._phi_rad = 'arb'

        """ Constants of dispersion formula """
        # For ordinary ray
        self._A1_o = 1.4313493
        self._B1_o = 0.0726631
        self._A2_o = 0.65054713
        self._B2_o = 0.1193242
        self._A3_o = 5.3414021
        self._B3_o = 18.028251
        # For extraordinary ray
        self._A1_e = 1.5039759
        self._B1_e = 0.0740288
        self._A2_e = 0.55069141
        self._B2_e = 0.1216529
        self._A3_e = 6.5927379
        self._B3_e = 20.072248

    def n_o_expr(self):
        """ Sympy expression, dispersion formula for o-wave """
        return sympy.sqrt(1 + self._A1_o * wl**2 / (wl**2 - self._B1_o**2) + self._A2_o * wl**2 / (wl**2 - self._B2_o**2) + self._A3_o * wl**2 / (wl**2 - self._B3_o**2))

    def n_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-wave """
        return sympy.sqrt(1 + self._A1_e * wl**2 / (wl**2 - self._B1_e**2) + self._A2_e * wl**2 / (wl**2 - self._B2_e**2) + self._A3_e * wl**2 / (wl**2 - self._B3_e**2))

    def n_expr(self, pol):
        """ Sympy expression, dispersion formula of a general ray with an angle theta to optic axis. If theta = 0, this expression reduces to 'no_expre'. """
        if pol == 'o':
            return self.n_o_expr()
        elif pol == 'e':
            return self.n_e_expr() / sympy.sqrt( sympy.sin(theta)**2 + (self.n_e_expr()/self.n_o_expr())**2 * sympy.cos(theta)**2 )
        else:
            raise ValueError("pol = '%s' must be 'o' or 'e'" % pol)
