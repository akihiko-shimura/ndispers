import sympy
from ndispers._baseclass import Medium, wl, phi, theta, T, pi
from ndispers.helper import vars2

class KTP(Medium):
    """
    KTP (K Ti O P O_4, Potassium Titanyl Phosphate) crystal

    - Point group : 2mm  (C_{2v})
    - Crystal system : Orthorhombic
    - Dielectric principal axes, x // a, y // b, z // c
    - Biaxial, with two optic axes in xz plane, symmetric with respect to z-axis
    - Tranparency range : 0.35 to 4.5 µm

    Sellmeier equation
    ------------------
    n(wl) = sqrt(A_i + B_i/(wl**2 - C_i) - D_i/(wl**2 - E_i))  for i = x, y, z
    
    Thermo-optic coefficient
    -------------------------
    dn/dT = (At_i/wl**3 + Bt_i/wl**2 + Ct_i/wl + Dt_i)*1e-5 for i = x,y,x

    Validity range
    ---------------
    dn/dT : 0.43 to 1.58 µm

    Ref
    ----
    Kato, Kiyoshi, and Eiko Takaoka. "Sellmeier and thermo-optic dispersion formulas for KTP." Applied optics 41.24 (2002): 5040-5044.
    
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

    Usage
    ------
    >>> ktp_xy = ndispers.media.crystals.KTP_xy()
    >>> ktp_xy.n(0.6, 0.3*pi, 40, pol='e') # for xy plane, 2nd argument is phi_rad. theta_rad is fixed at 0.5*pi.
    
    """

    __slots__ = ["_A_x", "_B_x", "_C_x", "_D_x", "_E_x",
                 "_A_y", "_B_y", "_C_y", "_D_y", "_E_y",
                 "_A_z", "_B_z", "_C_z", "_D_z", "_E_z",
                 "_dndT_x", "_dndT_y", "_dndT_z"]
    
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
        self._dndT_x = (0.1717/wl**3 - 0.5353/wl**2 + 0.8416/wl + 0.1627)*1e-5 #1/K
        self._dndT_y = (0.1997/wl**3 - 0.4063/wl**2 + 0.5154/wl + 0.5425)*1e-5 #1/K
        self._dndT_y = (0.9221/wl**3 - 2.9220/wl**2 + 3.6677/wl - 0.1897)*1e-5 #1/K
    
    @property
    def symbols(self):
        return [wl, theta, phi, T]

    @property
    def constants(self):
        print(vars2(self))
    
    def n_x_expr(self):
        """ sympy expresssion, dispersion formula of x-axis (principal dielectric axis) """
        return sympy.sqrt(self._A_x + self._B_x/(wl**2 - self._C_x) - self._D_x/(wl**2 - self._E_x)) + self._dndT_x * (T - 20)
    
    def n_y_expr(self):
        """ sympy expresssion, dispersion formula of y-axis (principal dielectric axis) """
        return sympy.sqrt(self._A_y + self._B_y/(wl**2 - self._C_y) - self._D_y/(wl**2 - self._E_y)) + self._dndT_y * (T - 20)

    def n_z_expr(self):
        """ sympy expresssion, dispersion formula of z-axis (principal dielectric axis) """
        return sympy.sqrt(self._A_z + self._B_z/(wl**2 - self._C_z) - self._D_z/(wl**2 - self._E_z)) + self._dndT_z * (T - 20)


class KTP_xy(KTP):
    __slots__ = ["_KTP_xy__plane", "_KTP_xy__theta_rad", "_KTP_xy__phi_rad"]

    def __init__(self):
        super().__init__()
        self._KTP_xy__plane = 'xy'
        self._KTP_xy__theta_rad = 0.5*pi
        self._KTP_xy__phi_rad = 'var'
    
    @property
    def help(self):
        print(super().__doc__)

    @property
    def plane(self):
        return self._KTP_xy__plane

    @property
    def theta_rad(self):
        return self._KTP_xy__theta_rad

    @property
    def phi_rad(self):
        return self._KTP_xy__phi_rad

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

    def n(self, wl_um, phi_rad, T_degC, pol='o'):
        """
        Refractive index in xy plane.

        input
        ------
        wl_um     :  float, wavelength in µm
        pol       :  str, 'o' or 'e', polarization of light
        phi_rad   :  float, 0 to 2pi radians
        (Note: theta_rad is fixed at 0.5*pi in xy principal plane.)

        return
        -------
        Refractive index, float

        """
        return super().n(wl_um, 0.5*pi, phi_rad, T_degC, pol=pol)

    def dn_wl(self, wl_um, phi_rad, T_degC, pol='o'):
        return super().dn_wl(wl_um, pol, 0.5*pi, phi_rad,  T_degC, pol=pol)
    
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
        return super().ng(wl_um, .5*pi, phi_rad, T_degC, pol=pol)

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


class KTP_yz(KTP):
    __slots__ = ["_KTP_yz__plane", "_KTP_yz__theta_rad", "_KTP_yz__phi_rad"]

    def __init__(self):
        super().__init__()
        self._KTP_yz__plane = 'yz'
        self._KTP_yz__phi_rad = 0.5*pi
        self._KTP_yz__theta_rad = 'var'
    
    @property
    def help(self):
        print(super().__doc__)
    
    @property
    def plane(self):
        return self._KTP_yz__plane

    @property
    def theta_rad(self):
        return self._KTP_yz__theta_rad

    @property
    def phi_rad(self):
        return self._KTP_yz__phi_rad

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
        wl_um     :  float, wavelength in µm
        pol       :  str, 'o' or 'e', polarization of light
        theta_rad :  float, 0 to 2pi radians
        (Note: phi_rad is fixed at 0.5*pi in yz principal plane.)

        return
        -------
        Refractive index, float

        """
        return super().n(wl_um, theta_rad, 0.5*pi, T_degC, pol=pol)

    def dn_wl(self, wl_um, theta_rad, T_degC, pol='o'):
        return super().dn_wl(wl_um, theta_rad, 0.5*pi, T_degC, pol=pol)
    
    def d2n_wl(self, wl_um, theta_rad, T_degC, pol='o'):
        return super().d2n_wl(wl_um, theta_rad, 0.5*pi, T_degC, pol=pol)

    def d3n_wl(self, wl_um, theta_rad, T_degC, pol='o'):
        return super().d3n_wl(wl_um, theta_rad, 0.5*pi, T_degC, pol=pol)

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


class KTP_zx(KTP):
    __slots__ = ["_KTP_zx__plane", "_KTP_zx__theta_rad", "_KTP_zx__phi_rad"]

    def __init__(self):
        super().__init__()
        self._KTP_zx__plane = 'zx'
        self._KTP_zx__theta_rad = 'arb'
        self._KTP_zx__phi_rad = 0.5*pi
    
    @property
    def help(self):
        print(super().__doc__)
    
    @property
    def plane(self):
        return self._KTP_zx__plane

    @property
    def theta_rad(self):
        return self._KTP_zx__theta_rad

    @property
    def phi_rad(self):
        return self._KTP_zx__phi_rad
    
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
        Refractive index in zx plane.

        input
        ------
        wl_um     :  float, wavelength in µm
        pol       :  str, 'o' or 'e', polarization of light
        theta_rad :  float, 0 to 2pi radians
        (Note: phi_rad is fixed at 0.5*pi in zx principal plane.)

        return
        -------
        Refractive index, float
        
        """
        return super().n(wl_um, theta_rad, 0.5*pi, T_degC, pol=pol)

    def dn_wl(self, wl_um, theta_rad, T_degC, pol='o'):
        return super().dn_wl(wl_um, theta_rad, 0.5*pi, T_degC, pol=pol)
    
    def d2n_wl(self, wl_um, theta_rad, T_degC, pol='o'):
        return super().d2n_wl(wl_um, theta_rad, 0.5*pi, T_degC, pol=pol)

    def d3n_wl(self, wl_um, theta_rad, T_degC, pol='o'):
        return super().d3n_wl(wl_um, theta_rad, 0.5*pi, T_degC, pol=pol)

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