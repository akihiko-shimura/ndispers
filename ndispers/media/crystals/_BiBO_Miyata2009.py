import sympy
from ndispers._baseclass import wl, phi, theta, T, pi
from ndispers.groups import Biax_2
from ndispers.helper import vars2

class BiBO(Biax_2):
    """
    BiBO (BiB₃O₆, bismuth triborate) crystal

    - Point group : 2  (C2)
    - Crystal system : Monoclinic
    - Dielectric principal axes: x // b (the two-fold axis); y and z lie in the
      a-c plane at ∠(a, z) = 31.1°, ∠(c, y) = 46.72° (632.8 nm) and rotate about
      b by about ±1° across the transparency range - neglected here
    - Positive biaxial (n_x < n_y < n_z), optic axes in the xz plane
    - Transparency range : 0.286 to 2.6 µm

    Sellmeier equation
    ------------------
        n(wl)**2 = A_i + B_i / (wl**2 - C_i) - D_i * wl**2   for i = x, y, z
        (B, C in µm**2, D in µm**-2)

    Validity range
    ---------------
    0.326 to 3.083 µm

    Note
    ----
    The Sellmeier equation has no temperature term: T_degC is accepted for
    signature uniformity and ignored, and dndT returns 0. Petrov et al. 2010
    (Table 3) list dn/dT at 1064 nm as +3.5e-6 (x), -5.6e-6 (y), -6.8e-6 (z) /K.

    The nonlinear tensor components are those of the dielectric xyz frame
    with x the two-fold axis: measured by Hellwig et al. 1999 (Maker fringes,
    1079.5 nm), transformed into the xyz frame by Ghotbi & Ebrahim-Zadeh 2004
    and averaged under Kleinman symmetry by Tzankov & Petrov 2005, as quoted
    in Petrov et al. 2010: d11 = 2.53, d12 = 3.2, d13 = -1.76, d14 = 1.66 pm/V.
    d_eff in the principal planes reproduces Table 1 of Petrov et al. 2010
    (walk-off added); in the xz plane the interaction type changes across the
    optic-axis angle Omega, and in the yz plane the three terms carry relative
    signs that matter.

    Two conventions for naming BiBO's planes are in circulation. This class
    follows the dielectric frame with the two-fold axis along x, in which the
    zx plane is the x-z plane of Petrov et al. 2010 (its type-I ooe locus runs
    1181.5 to 2272.8 nm and the eeo branch starts at 614.3 nm, both reproduced
    here). The cut usually quoted for 1064 nm SHG - "xz, theta = 168.5 deg,
    d_eff = 3.5 pm/V" - is written in the other convention (two-fold axis along
    y) and is this class's **yz** plane at theta = 168.8 deg, where d_eff comes
    out 3.62 pm/V.

    Because the group has no mirror plane through the z axis, theta and
    180 deg - theta phase-match alike (the index is the same) but give
    different d_eff: for 800 nm type-I SHG in the yz plane, theta = 28.8 deg
    gives 1.1 pm/V and 151.2 deg gives 3.9 pm/V. pmAngles_sfg returns the
    0-90 deg solution; evaluate deff_sfg at pi - theta as well.

    In the current version, biaxial crystals are limited to the principal
    dielectric planes, xy, yz or zx. Correspondence between principal plane,
    polarization orientations of o-wave and e-wave, polar (theta) and
    azimuthal (phi) angles of a wavevector with respect to z and x principal
    axes, respectively, are shown in the table below ('var' marks the angle
    passed as the method argument).

    | plane | o-wave | e-wave | theta | phi  |
    |-------|--------|--------|-------|------|
    | xy    | z      | xy     | pi/2  | var  |
    | yz    | x      | yz     | var   | pi/2 |
    | zx    | y      | zx     | var   | 0    |

    Ref
    ---
    Sellmeier equation (Table 2 of the review, from Miyata et al.):
      Miyata, K., Umemura, N., & Kato, K. (2009). Phase-matched pure χ(3) third-harmonic generation in noncentrosymmetric BiB3O6. Optics Letters, 34(4), 500-502. https://doi.org/10.1364/ol.34.000500
      Petrov, V., Ghotbi, M., Kokabee, O., et al. (2010). Femtosecond nonlinear frequency conversion based on BiB3O6. Laser & Photonics Reviews, 4(1), 53-98. https://doi.org/10.1002/lpor.200810075
    Nonlinear optical coefficients:
      Hellwig, H., Liebertz, J., & Bohatý, L. (1999). Exceptional large nonlinear optical coefficients in the monoclinic bismuth borate BiB3O6 (BIBO). Solid State Communications, 109(4), 249-251. https://doi.org/10.1016/s0038-1098(98)00538-9
    Transformation of the tensor to the xyz frame:
      Ghotbi, M., & Ebrahim-Zadeh, M. (2004). Optical second harmonic generation properties of BiB3O6. Optics Express, 12(24), 6002-6019. https://doi.org/10.1364/opex.12.006002
    d_eff expressions:
      Tzankov, P., & Petrov, V. (2005). Effective second-order nonlinearity in acentric optical crystals with low symmetry. Applied Optics, 44(32), 6971-6985. https://doi.org/10.1364/ao.44.006971
    """

    __slots__ = ["_A_x", "_B_x", "_C_x", "_D_x",
                 "_A_y", "_B_y", "_C_y", "_D_y",
                 "_A_z", "_B_z", "_C_z", "_D_z"]

    _d_ref = {"d11": (2.53, 1.0795, 1.0795),
              "d12": (3.2, 1.0795, 1.0795),
              "d13": (-1.76, 1.0795, 1.0795),
              "d14": (1.66, 1.0795, 1.0795)}
    _d_note = ("Hellwig et al. 1999 (Maker fringes at 1079.5 nm), transformed to the dielectric "
               "xyz frame (Ghotbi & Ebrahim-Zadeh 2004) and Kleinman-averaged (Tzankov & Petrov "
               "2005), as quoted by Petrov et al. 2010: d14 = d25 = d36, d12 = d26, d13 = d35, "
               "deviations from Kleinman symmetry within the 10% experimental error. Relative "
               "signs as given there. theta and pi - theta give different d_eff. Miller scaling "
               "untested for BiBO.")

    def __init__(self):
        super().__init__()

        """ Constants of dispersion formula, Petrov et al. 2010 Table 2 (from Miyata et al. 2009) """
        # for x-axis
        self._A_x = 3.0759
        self._B_x = 0.03169
        self._C_x = 0.03323
        self._D_x = 0.01402
        # for y-axis
        self._A_y = 3.1698
        self._B_y = 0.03666
        self._C_y = 0.03599
        self._D_y = 0.01819
        # for z-axis
        self._A_z = 3.6546
        self._B_z = 0.05116
        self._C_z = 0.03713
        self._D_z = 0.02299

    def n_x_expr(self):
        """ sympy expresssion, dispersion formula of x-axis (principal dielectric axis) """
        return sympy.sqrt(self._A_x + self._B_x / (wl**2 - self._C_x) - self._D_x * wl**2)

    def n_y_expr(self):
        """ sympy expresssion, dispersion formula of y-axis (principal dielectric axis) """
        return sympy.sqrt(self._A_y + self._B_y / (wl**2 - self._C_y) - self._D_y * wl**2)

    def n_z_expr(self):
        """ sympy expresssion, dispersion formula of z-axis (principal dielectric axis) """
        return sympy.sqrt(self._A_z + self._B_z / (wl**2 - self._C_z) - self._D_z * wl**2)

class BiBO_xy(BiBO):
    __slots__ = []

    def __init__(self):
        super().__init__()
        self._plane = 'xy'
        self._theta_rad = 0.5*pi
        self._phi_rad = 'var'


    def n_o_expr(self):
        """ sympy expresssion,
        dispersion formula for o-wave polarization for a given principal plane """
        return super().n_z_expr()

    def n_e_expr(self):
        """ sympy expresssion,
        dispersion formula for e-wave polarization for a given principal plane """
        return super().n_x_expr() * super().n_y_expr() / sympy.sqrt( super().n_x_expr()**2 * sympy.cos(phi)**2 + super().n_y_expr()**2 * sympy.sin(phi)**2 )

    def n_expr(self, pol):
        """ sympy expresssion,
        dispersion formula for a given polarization """
        if pol == 'o':
            return self.n_o_expr()
        elif pol == 'e':
            return self.n_e_expr()
        else:
            raise ValueError("pol = '%s' must be 'o' or 'e'" % pol)


class BiBO_yz(BiBO):
    __slots__ = []

    def __init__(self):
        super().__init__()
        self._plane = 'yz'
        self._phi_rad = 0.5*pi
        self._theta_rad = 'var'


    def n_o_expr(self):
        """ sympy expresssion,
        dispersion formula for o-wave polarization for yx principal plane """
        return super().n_x_expr()

    def n_e_expr(self):
        """ sympy expresssion,
        dispersion formula for e-wave polarization for yz principal plane """
        return super().n_y_expr() * super().n_z_expr() / sympy.sqrt( super().n_y_expr()**2 * sympy.sin(theta)**2 + super().n_z_expr()**2 * sympy.cos(theta)**2 )

    def n_expr(self, pol):
        """ sympy expresssion,
        dispersion formula for a given polarization """
        if pol == 'o':
            return self.n_o_expr()
        elif pol == 'e':
            return self.n_e_expr()
        else:
            raise ValueError("pol = '%s' must be 'o' or 'e'" % pol)


class BiBO_zx(BiBO):
    __slots__ = []

    def __init__(self):
        super().__init__()
        self._plane = 'zx'
        self._theta_rad = 'var'
        self._phi_rad = 0.0


    def n_o_expr(self):
        """ sympy expresssion,
        dispersion formula for o-wave polarization for zx principal plane """
        return super().n_y_expr()

    def n_e_expr(self):
        """ sympy expresssion,
        dispersion formula for e-wave polarization for zx principal plane """
        return super().n_z_expr() * super().n_x_expr() / sympy.sqrt( super().n_z_expr()**2 * sympy.cos(theta)**2 + super().n_x_expr()**2 * sympy.sin(theta)**2 )

    def n_expr(self, pol):
        """ sympy expresssion,
        dispersion formula for a given polarization """
        if pol == 'o':
            return self.n_o_expr()
        elif pol == 'e':
            return self.n_e_expr()
        else:
            raise ValueError("pol = '%s' must be 'o' or 'e'" % pol)
