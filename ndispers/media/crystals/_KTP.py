from ndispers._sym import sympy
from ndispers._baseclass import wl, phi, theta, T, pi
from ndispers.groups import Biax_mm2
from ndispers.helper import vars2

class KTP(Biax_mm2):
    """
    KTP (KTiOPO₄, potassium titanyl phosphate) crystal

    - Point group : mm2  (C2v)
    - Crystal system : Orthorhombic
    - Dielectric principal axes, x // a, y // b, z // c
    - Biaxial, with two optic axes in xz plane, symmetric with respect to z-axis
    - Transparency range : 0.35 to 4.5 µm

    Sellmeier equation
    ------------------
        n(wl) = sqrt(A_i + B_i/(wl**2 - C_i) + D_i/(wl**2 - E_i))  for i = x, y, z

    Thermo-optic coefficient
    -------------------------
        dn/dT = (At_i/wl**3 + Bt_i/wl**2 + Ct_i/wl + Dt_i)*1e-5 for i = x,y,x

    The temperature correction dn/dT * (T - 20) vanishes at 20 degC.

    Validity range
    ---------------
        dn/dT : 0.43 to 1.58 µm

    Note
    ----
    In the current version, biaxial crystals are limited to the principal dielectric planes,
    xy, yz or zx planes. In other words, a wavevector of light must be within any one of
    the three planes. Correspondence between principal plane, polarization orientations of
    o-wave and e-wave, polar (theta) and azimuthal (phi) angles of a wavevector with respect
    to z and x principal axes, respectively, are shown in the table below
    ('var' marks the angle passed as the method argument).

    | plane | o-wave | e-wave | theta | phi  |
    |-------|--------|--------|-------|------|
    | xy    | z      | xy     | pi/2  | var  |
    | yz    | x      | yz     | var   | pi/2 |
    | zx    | y      | zx     | var   | 0    |

    Ref
    ---
    Kato, K., & Takaoka, E. (2002). Sellmeier and thermo-optic dispersion formulas for KTP. Applied Optics, 41(24), 5040-5044. https://doi.org/10.1364/ao.41.005040
    """

    _wl_range = (0.35, 4.5)  # um, Sellmeier validity (see docstring)

    __slots__ = ["_A_x", "_B_x", "_C_x", "_D_x", "_E_x",
                 "_A_y", "_B_y", "_C_y", "_D_y", "_E_y",
                 "_A_z", "_B_z", "_C_z", "_D_z", "_E_z"]

    # mm2 tensor in the crystallographic frame (a, b, c) with c polar;
    # dielectric x // a, y // b, z // c
    _mm2_axes = ("x", "y", "z")
    _d_ref = {"d31": (2.2, 1.064, 1.064),
              "d32": (3.7, 1.064, 1.064),
              "d33": (14.6, 1.064, 1.064)}
    _d_note = ("Shoji et al. 1997, absolute, 1.064 um SHG, converted to the conventional "
               "x // a, y // b frame (their Table 10 follows Roberts's a/b interchange: "
               "their d31, d15 are d32, d24 here). Shoji et al. found Miller's delta of "
               "KTP to change by 35% between 1.06 and 1.31 um, so treat Miller scaling "
               "for KTP as rough; Alford & Smith 2001 judge it acceptable once earlier "
               "data are reinterpreted.")

    def __init__(self):
        super().__init__()

        # for x-axis
        self._A_x = 3.29100
        self._B_x = 0.04140
        self._C_x = 0.03978
        self._D_x = 9.35522
        self._E_x = 31.45571
        # for y-axis
        self._A_y = 3.45018
        self._B_y = 0.04341
        self._C_y = 0.04597
        self._D_y = 16.98825
        self._E_y = 39.43799
        # z-axis
        self._A_z = 4.59423
        self._B_z = 0.06206
        self._C_z = 0.04763
        self._D_z = 110.80672
        self._E_z = 86.12171
        #dn/dT


    @property
    def _dndT_x(self):
        # a sympy expression, built on demand so that constructing
        # the medium does not need sympy
        return (0.1717/wl**3 - 0.5353/wl**2 + 0.8416/wl + 0.1627)*1e-5  # 1/K

    @property
    def _dndT_y(self):
        # a sympy expression, built on demand so that constructing
        # the medium does not need sympy
        return (0.1997/wl**3 - 0.4063/wl**2 + 0.5154/wl + 0.5425)*1e-5  # 1/K

    @property
    def _dndT_z(self):
        # a sympy expression, built on demand so that constructing
        # the medium does not need sympy
        return (0.9221/wl**3 - 2.9220/wl**2 + 3.6677/wl - 0.1897)*1e-5  # 1/K

    def n_x_expr(self):
        """ sympy expresssion, dispersion formula of x-axis (principal dielectric axis) """
        return sympy.sqrt(self._A_x + self._B_x/(wl**2 - self._C_x) + self._D_x/(wl**2 - self._E_x)) + self._dndT_x * (T - 20)

    def n_y_expr(self):
        """ sympy expresssion, dispersion formula of y-axis (principal dielectric axis) """
        return sympy.sqrt(self._A_y + self._B_y/(wl**2 - self._C_y) + self._D_y/(wl**2 - self._E_y)) + self._dndT_y * (T - 20)

    def n_z_expr(self):
        """ sympy expresssion, dispersion formula of z-axis (principal dielectric axis) """
        return sympy.sqrt(self._A_z + self._B_z/(wl**2 - self._C_z) + self._D_z/(wl**2 - self._E_z)) + self._dndT_z * (T - 20)

class KTP_xy(KTP):
    __slots__ = []

    def __init__(self):
        super().__init__()
        self._plane = 'xy'
        self._theta_rad = 0.5*pi
        self._phi_rad = 'var'


    def n_o_expr(self):
        """ sympy expresssion,
        dispersion formula for o-wave polarization for a given principal plane
        """
        return super().n_z_expr()

    def n_e_expr(self):
        """ sympy expresssion,
        dispersion formula for e-wave polarization for a given principal plane """
        return super().n_x_expr() * super().n_y_expr() / sympy.sqrt( super().n_x_expr()**2 * sympy.cos(phi)**2 + super().n_y_expr()**2 * sympy.sin(phi)**2 )

    def n_expr(self, pol):
        """ sympy expresssion,
        dispersion formula for a given polarization
        """
        if pol == 'o':
            return self.n_o_expr()
        elif pol == 'e':
            return self.n_e_expr()
        else:
            raise ValueError("pol = '%s' must be 'o' or 'e'" % pol)


class KTP_yz(KTP):
    __slots__ = []

    def __init__(self):
        super().__init__()
        self._plane = 'yz'
        self._phi_rad = 0.5*pi
        self._theta_rad = 'var'


    def n_o_expr(self):
        """ sympy expresssion,
        dispersion formula for o-wave polarization for yx principal plane
        """
        return super().n_x_expr()

    def n_e_expr(self):
        """ sympy expresssion,
        dispersion formula for e-wave polarization for yz principal plane
        """
        return super().n_y_expr() * super().n_z_expr() / sympy.sqrt( super().n_y_expr()**2 * sympy.sin(theta)**2 + super().n_z_expr()**2 * sympy.cos(theta)**2 )

    def n_expr(self, pol):
        """ sympy expresssion,
        dispersion formula for a given polarization
        """
        if pol == 'o':
            return self.n_o_expr()
        elif pol == 'e':
            return self.n_e_expr()
        else:
            raise ValueError("pol = '%s' must be 'o' or 'e'" % pol)


class KTP_zx(KTP):
    __slots__ = []

    def __init__(self):
        super().__init__()
        self._plane = 'zx'
        self._theta_rad = 'var'
        self._phi_rad = 0.0


    def n_o_expr(self):
        """ sympy expresssion,
        dispersion formula for o-wave polarization for zx principal plane
        """
        return super().n_y_expr()

    def n_e_expr(self):
        """ sympy expresssion,
        dispersion formula for e-wave polarization for zx principal plane
        """
        return super().n_z_expr() * super().n_x_expr() / sympy.sqrt( super().n_z_expr()**2 * sympy.cos(theta)**2 + super().n_x_expr()**2 * sympy.sin(theta)**2 )

    def n_expr(self, pol):
        """ sympy expresssion,
        dispersion formula for a given polarization
        """
        if pol == 'o':
            return self.n_o_expr()
        elif pol == 'e':
            return self.n_e_expr()
        else:
            raise ValueError("pol = '%s' must be 'o' or 'e'" % pol)
