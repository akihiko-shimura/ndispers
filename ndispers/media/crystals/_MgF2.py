import sympy
from ndispers._baseclass import Medium, wl, phi, theta, T
from ndispers.helper import vars2

class MgF2(Medium):
    """
    MgF₂ (magnesium fluoride, sellaite) crystal

    - Point group : 4/mmm  (D4h); centrosymmetric, so no second-order nonlinearity
    - Crystal system : Tetragonal (rutile structure)
    - Dielectric principal axis, z // c-axis (x, y-axes are arbitrary)
    - Positive uniaxial, with optic axis parallel to z-axis
    - Transparency range : 0.13 to 7.7 µm (o-ray, Tropf et al. Table 20)

    Sellmeier equation
    ------------------
        n(wl)**2 = 1 + A1_i * wl**2 / (wl**2 - B1_i**2) + A2_i * wl**2 / (wl**2 - B2_i**2) + A3_i * wl**2 / (wl**2 - B3_i**2)   for i = o, e

    Validity range
    ---------------
    0.20 to 7.04 µm

    Note
    ----
    The Sellmeier equation has no temperature term: T_degC is accepted for
    signature uniformity and ignored, and dndT returns 0. Tropf et al.
    (Table 20, from Feldman et al. 1978) list dn_o/dT = +1.12e-6 /K at
    0.633 µm, +0.88e-6 at 1.15 µm and +1.19e-6 at 3.39 µm. The small positive
    birefringence (n_e - n_o = 0.012 at 589 nm) is what makes MgF₂ the usual
    material for UV/VUV wave plates and Brewster windows; woa_theta gives the
    walk-off of the e-ray.

    The e-ray's infrared pole is 23.771995 µm (Dodge 1984, as transcribed by
    refractiveindex.info, main/MgF2/nk/Dodge-e). The Handbook table prints
    12.771995 there (a typo in the Handbook); that value gives n_e = 1.3882 at
    589 nm and a birefringence of 0.0105, against the measured 1.3895 and
    0.0118 that 23.771995 reproduces, and the o-ray's matching pole is
    23.793604 - so 23.771995 is taken as the source value.

    Ref
    ---
    Sellmeier equation (as tabulated in):
      Tropf, W. J., Thomas, M. E., & Harris, T. J. (1995). Properties of crystals and glasses. In M. Bass (Ed.), Handbook of Optics, Vol. II, Chapter 33, Table 22. McGraw-Hill.
    Underlying measurements:
      Dodge, M. J. (1984). Refractive properties of magnesium fluoride. Applied Optics, 23(12), 1980-1985. https://doi.org/10.1364/ao.23.001980
    """
    __slots__ = ["_A1_o", "_B1_o", "_A2_o", "_B2_o", "_A3_o", "_B3_o",
                 "_A1_e", "_B1_e", "_A2_e", "_B2_e", "_A3_e", "_B3_e"]

    def __init__(self):
        super().__init__()
        self._plane = 'arb'
        self._theta_rad = 'var'
        self._phi_rad = 'arb'

        """ Constants of dispersion formula """
        # For ordinary ray
        self._A1_o = 0.48755108
        self._B1_o = 0.04338408
        self._A2_o = 0.39875031
        self._B2_o = 0.09461442
        self._A3_o = 2.3120353
        self._B3_o = 23.793604
        # For extraordinary ray
        self._A1_e = 0.41344023
        self._B1_e = 0.03684262
        self._A2_e = 0.50497499
        self._B2_e = 0.09076162
        self._A3_e = 2.4904862
        self._B3_e = 23.771995   # see Note: the Handbook prints 12.771995

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
