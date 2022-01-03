import sympy
from ndispers._baseclass import Medium, wl, phi, theta, T, pi
from ndispers.helper import vars2

class LBO(Medium):
    """
    LBO (Li B_3 O_5) crystal

    - Point group : mm2
    - Crystal system : orthorhombic 
    - Dielectric principal axes, x // a, y // -c, z // b
    - Biaxial, with two optic axes in xz plane, symmetric with respect to z-axis

    Dispersion formula for refractive index
    ---------------------------------------
    n(wl) = sqrt(A_i + B_i/(wl**2 - C_i) - D_i * wl**2 + E_i * wl**4 - F_i * wl**6)    for i = x,y,z

    Thermo-optic coefficient
    -------------------------
    dn/dT = F_i + G_i * wl   for i = x,y,z
    
    Validity range
    --------------

    Ref
    ---
    Kato, K., "Temperature-tuned 90-deg phase-matching properties of LiB3O5," IEEE Journal of Quantum Electronics 30, 2950-2952 (1994).
    Kato, K., S. G. Grechin, and N. Umemura. "New thermo-optic dispersion formula for LiB3O5." Laser Physics 28.9 (2018): 095403.

    Input
    -----
    plane  :  Principal dielectric plane which includes wave vector of light ("xy", "yz" or "xz")
    
    If plane == "xy", 
        o-ray polarization // z-axis, e-ray polarization in xy-plane, phi is variable and theta = 90 deg.
    If plane == "yz", 
        o-ray polarization // x-axis, e-ray polarization in yz-plane, phi = 90 deg and theta is variable.
    If plane == "xz", 
        o-ray polarization // y-axis, e-ray polarization in xz-plane, phi = 0 deg and theta is variable.

    Example
    -------
    >>> lbo_xy = ndispers.media.crystals.LBO_xy()
    >>> lbo_xy.n(0.6, 0.3*pi, 40, pol='e') # for xy plane, 2nd argument is phi_rad. theta_rad is fixed at 0.5*pi.
    
    """

    __slots__ = ["_A_x", "_B_x", "_C_x", "_D_x", "_E_x", "_F_x",
                 "_A_y", "_B_y", "_C_y", "_D_y", "_E_y", "_F_y",
                 "_A_z", "_B_z", "_C_z", "_D_z", "_E_z", "_F_z",
                 "_F_x", "_F_y", "_F_z", 
                 "_G_x", "_G_y", "_G_z"]
    
    def __init__(self):
        super().__init__()

        # for x-axis
        self._A_x = 2.4542
        self._B_x = 0.01125
        self._C_x = 0.01135
        self._D_x = 0.01388
        self._E_x = 0.0
        self._F_x = 0.0
        # for y-axis
        self._A_y = 2.5390
        self._B_y = 0.01277
        self._C_y = 0.01189
        self._D_y = 0.01849
        self._E_y = 4.3025e-5
        self._F_y = 2.9131e-5
        # z-axis
        self._A_z = 2.5865
        self._B_z = 0.01310
        self._C_z = 0.01223
        self._D_z = 0.01862
        self._E_z = 4.5778e-5
        self._F_z = 3.2526e-5
        # dn/dT
        self._F_x = +0.2300e-5
        self._G_x = -0.3760e-5
        self._F_y = -1.9318e-5
        self._G_y = +0.5779e-5
        self._F_z = -1.1569e-5
        self._G_z = +0.4073e-5
    
    @property
    def constants(self):
        print(vars2(self))
    
    @property
    def symbols(self):
        return [wl, theta, phi, T]
    
    def _n_T20_x_expr(self):
        """ Sympy expression, dispersion formula for x-axis (principal dielectric axis) at 20degC """
        return sympy.sqrt(self._A_x + self._B_x/(wl**2 - self._C_x) - self._D_x * wl**2 + self._E_x * wl**4 - self._F_x * wl**6)
    
    def _n_T20_y_expr(self):
        """ Sympy expression, dispersion formula for y-axis (principal dielectric axis) at 20degC """
        return sympy.sqrt(self._A_y + self._B_y/(wl**2 - self._C_y) - self._D_y * wl**2 + self._E_y * wl**4 - self._F_y * wl**6)
    
    def _n_T20_z_expr(self):
        """ Sympy expression, dispersion formula for x-axis (principal dielectric axis) at 20degC """
        return sympy.sqrt(self._A_z + self._B_z/(wl**2 - self._C_z) - self._D_z * wl**2 + self._E_z * wl**4 - self._F_z * wl**6)
    
    def dndT_x_expr(self):
        return self._F_x + self._G_x * wl
    
    def dndT_y_expr(self):
        return self._F_y + self._G_y * wl
    
    def dndT_z_expr(self):
        return self._F_z + self._G_z * wl

    def n_x_expr(self):
        """ sympy expresssion, dispersion formula of x-axis (principal dielectric axis) """
        return self._n_T20_x_expr() + self.dndT_x_expr() * (T - 20)
    
    def n_y_expr(self):
        """ sympy expresssion, dispersion formula of y-axis (principal dielectric axis) """
        return self._n_T20_y_expr() + self.dndT_y_expr() * (T - 20)

    def n_z_expr(self):
        """ sympy expresssion, dispersion formula of z-axis (principal dielectric axis) """
        return self._n_T20_z_expr() + self.dndT_z_expr() * (T - 20)


class LBO_xy(LBO):
    __slots__ = ["_LBO_xy__plane", "_LBO_xy__theta_rad", "_LBO_xy__phi_rad"]

    def __init__(self):
        super().__init__()
        self._LBO_xy__plane = 'xy'
        self._LBO_xy__theta_rad = 0.5*pi
        self._LBO_xy__phi_rad = 'var'
    
    @property
    def help(self):
        print(super().__doc__)

    @property
    def plane(self):
        return self._LBO_xy__plane

    @property
    def theta_rad(self):
        return self._LBO_xy__theta_rad

    @property
    def phi_rad(self):
        return self._LBO_xy__phi_rad

    @property
    def constants(self):
        print({**vars2(super()), **vars2(self)})

    def n_o_expr(self):
        """ sympy expresssion, 
        dispersion formula for o-ray polarization for a given principal plane """
        return super().n_z_expr()
    
    def n_e_expr(self):
        """ sympy expresssion, 
        dispersion formula for e-ray polarization for a given principal plane """
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

    def n(self, wl_um, phi_rad, T_degC, pol='o'):
        """
        Refractive index in xy plane.

        input
        ------
        wl_um     :  float or array_like, wavelength in um
        phi_rad   :  float or array_like, polar angle in radians
        T_degC    :  float or array_like, temperature of crystal in degree C.
        (Note: theta_rad is fixed at 0.5*pi in xy principal plane.)

        return
        -------
        Refractive index, float or array_like
        """
        return super().n(wl_um, 0.5*pi, phi_rad, T_degC, pol=pol)

    def dn_wl(self, wl_um, phi_rad, T_degC, pol='o'):
        return super().dn_wl(wl_um, 0.5*pi, phi_rad, T_degC, pol=pol)
    
    def d2n_wl(self, wl_um, phi_rad, T_degC, pol='o'):
        return super().d2n_wl(wl_um, 0.5*pi, phi_rad, T_degC, pol=pol)

    def d3n_wl(self, wl_um, phi_rad, T_degC, pol='o'):
        return super().d3n_wl(wl_um, 0.5*pi, phi_rad, T_degC, pol=pol)

    def GD(self, wl_um, phi_rad, T_degC, pol='o'):
        """Group Delay [fs/mm]"""
        return super().GD(wl_um, 0.5*pi, phi_rad, T_degC, pol=pol)
    
    def GV(self, wl_um, phi_rad, T_degC, pol='o'):
        """Group Velocity [um/fs]"""
        return super().GV(wl_um, 0.5*pi, phi_rad, T_degC, pol=pol)
    
    def ng(self, wl_um, phi_rad, T_degC, pol='o'):
        """Group index, c/Group velocity"""
        return super().ng(wl_um, 0.5*pi, phi_rad, T_degC, pol=pol)

    def GVD(self, wl_um, phi_rad, T_degC, pol='o'):
        """Group Delay Dispersion [fs^2/mm]"""
        return super().GVD(wl_um, 0.5*pi, phi_rad, T_degC, pol=pol)

    def TOD(self, wl_um, phi_rad, T_degC, pol='o'):
        """Third Order Dispersion [fs^3/mm]"""
        return super().TOD(wl_um, 0.5*pi, phi_rad, T_degC, pol=pol)
    
    def woa_theta(self, wl_um, phi_rad, T_degC, pol='e'):
        return super().woa_theta(wl_um, 0.5*pi, phi_rad, T_degC, pol=pol)
    
    def woa_phi(self, wl_um, phi_rad, T_degC, pol='e'):
        return super().woa_phi(wl_um, 0.5*pi, phi_rad, T_degC, pol=pol)
    
    def dndT(self, wl_um, phi_rad, T_degC, pol='o'):
        return super().dndT(wl_um, 0.5*pi, phi_rad, T_degC, pol=pol)


class LBO_yz(LBO):
    __slots__ = ["_LBO_yz__plane", "_LBO_yz__theta_rad", "_LBO_yz__phi_rad"]

    def __init__(self):
        super().__init__()
        self._LBO_yz__plane = 'yz'
        self._LBO_yz__phi_rad = 0.5*pi
        self._LBO_yz__theta_rad = 'var'
    
    @property
    def help(self):
        print(super().__doc__)

    @property
    def plane(self):
        return self._LBO_yz__plane

    @property
    def theta_rad(self):
        return self._LBO_yz__theta_rad

    @property
    def phi_rad(self):
        return self._LBO_yz__phi_rad

    @property
    def constants(self):
        print({**vars2(super()), **vars2(self)})

    def n_o_expr(self):
        """ sympy expresssion, 
        dispersion formula for o-ray polarization for yx principal plane """
        return super().n_x_expr()
    
    def n_e_expr(self):
        """ sympy expresssion, 
        dispersion formula for e-ray polarization for yz principal plane """
        return super().n_y_expr() * super().n_z_expr() / sympy.sqrt( super().n_y_expr()**2 * sympy.cos(phi)**2 + super().n_z_expr()**2 * sympy.sin(phi)**2 )

    def n_expr(self, pol):
        """ sympy expresssion, 
        dispersion formula for a given polarization """
        if pol == 'o':
            return self.n_o_expr()
        elif pol == 'e':
            return self.n_e_expr()
        else:
            raise ValueError("pol = '%s' must be 'o' or 'e'" % pol)

    def n(self, wl_um, theta_rad, T_degC, pol='o'):
        """
        Refractive index in yz plane.

        input
        ------
        wl_um     :  float or array_like, wavelength in um
        theta_rad   :  float or array_like, azimuthal angle in radians
        T_degC    :  float or array_like, temperature of crystal in degree C.
        (Note: phi_rad is fixed at 0.5*pi in xy principal plane.)

        return
        -------
        Refractive index, float or array_like
        """
        return super().n(wl_um, theta_rad, 0.5*pi, T_degC, pol=pol)

    def dn_wl(self, wl_um, theta_rad, T_degC, pol='o'):
        return super().dn_wl(wl_um, theta_rad, 0.5*pi, T_degC, pol=pol)
    
    def d2n_wl(self, wl_um, theta_rad, T_degC, pol='o'):
        return super().d2n_wl(wl_um, theta_rad, 0.5*pi, T_degC, pol=pol)

    def d3n_wl(self, wl_um, theta_rad, T_degC, pol='o'):
        return super().d3n_wl(wl_um, theta_rad, 0.5*pi, pol=pol)

    def GD(self, wl_um, theta_rad, T_degC, pol='o'):
        """Group Delay [fs/mm]"""
        return super().GD(wl_um, theta_rad, 0.5*pi, T_degC, pol=pol)
    
    def GV(self, wl_um, theta_rad, T_degC, pol='o'):
        """Group Velocity [um/fs]"""
        return super().GV(wl_um, theta_rad, 0.5*pi, T_degC, pol=pol)
    
    def ng(self, wl_um, theta_rad, T_degC, pol='o'):
        """Group index, c/Group velocity"""
        return super().ng(wl_um, theta_rad, 0.5*pi, T_degC, pol=pol)

    def GVD(self, wl_um, theta_rad, T_degC, pol='o'):
        """Group Delay Dispersion [fs^2/mm]"""
        return super().GVD(wl_um, theta_rad, 0.5*pi, T_degC, pol=pol)

    def TOD(self, wl_um, theta_rad, T_degC, pol='o'):
        """Third Order Dispersion [fs^3/mm]"""
        return super().TOD(wl_um, theta_rad, 0.5*pi, T_degC, pol=pol)
    
    def woa_theta(self, wl_um, theta_rad, T_degC, pol='e'):
        return super().woa_theta(wl_um, theta_rad, 0.5*pi, T_degC, pol=pol)
    
    def woa_phi(self, wl_um, theta_rad, T_degC, pol='e'):
        return super().woa_phi(wl_um, theta_rad, 0.5*pi, T_degC, pol=pol)
    
    def dndT(self, wl_um, theta_rad, T_degC, pol='o'):
        return super().dndT(wl_um, theta_rad, 0.5*pi, T_degC, pol=pol)


class LBO_zx(LBO):
    __slots__ = ["_LBO_zx__plane", "_LBO_zx__theta_rad", "_LBO_zx__phi_rad"]

    def __init__(self):
        super().__init__()
        self._LBO_zx__plane = 'zx'
        self._LBO_zx__theta_rad = 'var'
        self._LBO_zx__phi_rad = 0.5*pi
    
    @property
    def help(self):
        print(super().__doc__)

    @property
    def plane(self):
        return self._LBO_zx__plane

    @property
    def theta_rad(self):
        return self._LBO_zx__theta_rad

    @property
    def phi_rad(self):
        return self._LBO_zx__phi_rad

    @property
    def constants(self):
        print({**vars2(super()), **vars2(self)})

    def n_o_expr(self):
        """ sympy expresssion, 
        dispersion formula for o-ray polarization for zx principal plane """
        return super().n_y_expr()
    
    def n_e_expr(self):
        """ sympy expresssion, 
        dispersion formula for e-ray polarization for zx principal plane """
        return super().n_z_expr() * super().n_z_expr() / sympy.sqrt( super().n_z_expr()**2 * sympy.cos(phi)**2 + super().n_x_expr()**2 * sympy.sin(phi)**2 )

    def n_expr(self, pol):
        """ sympy expresssion, 
        dispersion formula for a given polarization """
        if pol == 'o':
            return self.n_o_expr()
        elif pol == 'e':
            return self.n_e_expr()
        else:
            raise ValueError("pol = '%s' must be 'o' or 'e'" % pol)

    def n(self, wl_um, theta_rad, T_degC, pol='o'):
        """
        Refractive index in yz plane.

        input
        ------
        wl_um     :  float or array_like, wavelength in um
        theta_rad   :  float or array_like, azimuthal angle in radians
        T_degC    :  float or array_like, temperature of crystal in degree C.
        (Note: phi_rad is fixed at 0.5*pi in xy principal plane.)

        return
        -------
        Refractive index, float or array_like
        """
        return super().n(wl_um, theta_rad, 0.5*pi, T_degC, pol=pol)

    def dn_wl(self, wl_um, theta_rad, T_degC, pol='o'):
        return super().dn_wl(wl_um, theta_rad, 0.5*pi, T_degC, pol=pol)
    
    def d2n_wl(self, wl_um, theta_rad, T_degC, pol='o'):
        return super().d2n_wl(wl_um, theta_rad, 0.5*pi, T_degC, pol=pol)

    def d3n_wl(self, wl_um, theta_rad, T_degC, pol='o'):
        return super().d3n_wl(wl_um, theta_rad, 0.5*pi, pol=pol)

    def GD(self, wl_um, theta_rad, T_degC, pol='o'):
        """Group Delay [fs/mm]"""
        return super().GD(wl_um, theta_rad, 0.5*pi, T_degC, pol=pol)
    
    def GV(self, wl_um, theta_rad, T_degC, pol='o'):
        """Group Velocity [um/fs]"""
        return super().GV(wl_um, theta_rad, 0.5*pi, T_degC, pol=pol)
    
    def ng(self, wl_um, theta_rad, T_degC, pol='o'):
        """Group index, c/Group velocity"""
        return super().ng(wl_um, theta_rad, 0.5*pi, T_degC, pol=pol)

    def GVD(self, wl_um, theta_rad, T_degC, pol='o'):
        """Group Delay Dispersion [fs^2/mm]"""
        return super().GVD(wl_um, theta_rad, 0.5*pi, T_degC, pol=pol)

    def TOD(self, wl_um, theta_rad, T_degC, pol='o'):
        """Third Order Dispersion [fs^3/mm]"""
        return super().TOD(wl_um, theta_rad, 0.5*pi, T_degC, pol=pol)
    
    def woa_theta(self, wl_um, theta_rad, T_degC, pol='e'):
        return super().woa_theta(wl_um, theta_rad, 0.5*pi, T_degC, pol=pol)
    
    def woa_phi(self, wl_um, theta_rad, T_degC, pol='e'):
        return super().woa_phi(wl_um, theta_rad, 0.5*pi, T_degC, pol=pol)
    
    def dndT(self, wl_um, theta_rad, T_degC, pol='o'):
        return super().dndT(wl_um, theta_rad, 0.5*pi, T_degC, pol=pol)
    