import sympy
from ndispers._baseclass import Medium, wl, phi, theta, T, pi
from ndispers.helper import vars2

class LBO(Medium):
    """
    LBO (Li B_3 O_5) crystal

    - Point group : 2mm  (C_{2v})
    - Crystal system : orthorhombic 
    - Dielectric principal axes, x // a, y // -c, z // b
    - Biaxial, with two optic axes in xz plane, symmetric with respect to z-axis

    Sellmeier equation
    ------------------
    n(wl) = sqrt( A_i + B_i/(1 - C_i/wl**2) + D_i/(1 - E_i/wl**2) ) + dn/dT * (T -20)  for i = x,y,z

    Thermo-optic coefficient
    ------------------------
    dn/dT = (G_i * R_i + H_i * R_i**2) / (2*n_i)  for i = x,y,z
    R_i = wl**2/(wl**2 - wl0_i**2), wl0_i = 0.0530, 0.0327, 0.0435 for i=x,y,z
    
    Validity range
    --------------

    Ref
    ---
    Ghosh, Gorachand. "Temperature dispersion of refractive indices in β-BaB2O4 and LiB3O5 crystals for nonlinear optical devices." Journal of applied physics 78.11 (1995): 6752-6760.

    Note
    ----
    In the current version, biaxial crystals are limited to the principal dielectric planes, 
    xy, yz or zx planes. In other words, a wavevector of light must be within any one of 
    the three planes. Correspondence between principal plane, polarization orientations of 
    o-wave and e-wave, polar (theta) and azimuthal (phi) angles of a wavevector with respect 
    to z and x principal axes, respectively, are shown in the table below.

    plane  |  o-wave  |  e-wave  |  theta  |   phi   |
    ================================================
    xy     |    z    |    xy   |   pi/2  |   var   |
    yz     |    x    |    yz   |   var   |   pi/2  |
    zx     |    y    |    zx   |   var   |    0    |
    ------------------------------------------------
    pi = 3.14159...
    var : variable

    Example
    -------
    >>> lbo_xy = ndispers.media.crystals.LBO_KK2018_xy()
    >>> lbo_xy.n(0.6, 0.3*pi, 40, pol='e') # for xy plane, 2nd argument is phi_rad. theta_rad is fixed at 0.5*pi.

    """

    __slots__ = ["_A_x", "_B_x", "_C_x", "_D_x", "_E_x",
                 "_A_y", "_B_y", "_C_y", "_D_y", "_E_y",
                 "_A_z", "_B_z", "_C_z", "_D_z", "_E_z",
                 "_G_x", "_H_x", "_R_x",
                 "_G_y", "_H_y", "_R_y",
                 "_G_z", "_H_z", "_R_z"]
    
    def __init__(self):
        super().__init__()

        # for x-axis
        self._A_x = 1.4426279
        self._B_x = 1.0109932
        self._C_x = 1.1210197e-2
        self._D_x = 1.2363218
        self._E_x = 91
        # for y-axis
        self._A_y = 1.5014015
        self._B_y = 1.0388217
        self._C_y = 1.2157100e-2
        self._D_y = 1.7567133
        self._E_y = 91
        # z-axis
        self._A_z = 1.4489240
        self._B_z = 1.1365228
        self._C_z = 1.1676746e-2
        self._D_z = 1.5830069
        self._E_z = 91
        # dn/dT
        self._G_x = -127.70167e-6
        self._G_y = (372.170 - 2.199e-1 * T + 1.1748e-3 * T**2 - 2.05077e-6 * T**3)*1e-6
        self._G_z = -446.95031e-6
        self._H_x = 122.13435e-6
        self._H_y = -415.10435e-6
        self._H_z = (410.66123 + 1.667e-1 * T - 5.1887e-4 * T**2 + 5.56251e-7 * T**3)*1e-6
        self._R_x = wl**2/(wl**2 - 0.0530**2)
        self._R_y = wl**2/(wl**2 - 0.0327**2)
        self._R_z = wl**2/(wl**2 - 0.0435**2)
    
    @property
    def constants(self):
        print(vars2(self))
    
    @property
    def symbols(self):
        return [wl, theta, phi, T]
    
    def _n_T20_x_expr(self):
        """ Sympy expression, dispersion formula for x-axis (principal dielectric axis) at 20degC"""
        return sympy.sqrt(self._A_x + self._B_x / (1 - self._C_x/wl**2) + self._D_x/(1 - self._E_x/wl**2))
    
    def _n_T20_y_expr(self):
        """ Sympy expression, dispersion formula for y-axis (principal dielectric axis) at 20degC"""
        return sympy.sqrt(self._A_y + self._B_y / (1 - self._C_y/wl**2) + self._D_y/(1 - self._E_y/wl**2))
    
    def _n_T20_z_expr(self):
        """ Sympy expression, dispersion formula for x-axis (principal dielectric axis) at 20degC"""
        return sympy.sqrt(self._A_z + self._B_z / (1 - self._C_z/wl**2) + self._D_z/(1 - self._E_z/wl**2))
    
    def dndT_x_expr(self):
        return (self._G_x * self._R_x + self._H_x * self._R_x**2) / (2*self._n_T20_x_expr())
    
    def dndT_y_expr(self):
        return (self._G_y * self._R_y + self._H_y * self._R_y**2) / (2*self._n_T20_y_expr())
    
    def dndT_z_expr(self):
        return (self._G_z * self._R_z + self._H_z * self._R_z**2) / (2*self._n_T20_z_expr())

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
        dispersion formula for o-wave polarization for a given principal plane
        """
        return super().n_z_expr()
    
    def n_e_expr(self):
        """ sympy expresssion, 
        dispersion formula for e-wave polarization for a given principal plane
        """
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

    def n(self, wl_um, phi_rad, T_degC, pol='o'):
        """
        Refractive index in xy plane.

        input
        ------
        wl_um     :  float or array_like, wavelength in µm
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
        """Group Velocity [µm/fs]"""
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

    def n(self, wl_um, theta_rad, T_degC, pol='o'):
        """
        Refractive index in yz plane.

        input
        ------
        wl_um     :  float or array_like, wavelength in µm
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
        """Group Velocity [µm/fs]"""
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

    def n(self, wl_um, theta_rad, T_degC, pol='o'):
        """
        Refractive index in yz plane.

        input
        ------
        wl_um     :  float or array_like, wavelength in µm
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
        """Group Velocity [µm/fs]"""
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
    