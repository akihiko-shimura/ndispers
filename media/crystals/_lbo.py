import sympy
from ndispers._baseclass import Medium, wl, phi, theta, pi

class LBO(Medium):
    """
    LBO (Li B_3 O_5) crystal

    - Point group : mm2
    - Crystal ststem : orthorhombic 
    - Dielecric principal axes, x // a, y // -c, z // b
    - Biaxial, with two optic axes in xz plane, symmetric with respect to z-axis

    Dispersion formula of refractive index
    ---------------------------------------
        n(wl_um) = sqrt(A_i + B_i/(wl_um**2 - C_i) - D_i * wl_um**2)  for i = x, y, z
    
    Validity range
    ---------------
    0.22 - 1.32 um

    Ref
    ----
    Kato, K. "Tunable UV generation to 0.2325 mu m in LiB/sub 3/O/sub 5." IEEE journal of quantum electronics 26.7 (1990): 1173-1175.
    
    Note
    -----
    The constants of dispersion formula are at temperature T = 20 degC.

    Input
    ------
    plane  :  Principal dielectric plane which includes wave vector of light ("xy", "yz" or "xz")
    
    If plane == "xy", 
        o-ray polarization // z-axis, e-ray polarization in xy-plane, phi is a variable and theta = 90 deg.
    If plane == "yz", 
        o-ray polarization // x-axis, e-ray polarization in yz-plane, phi = 90 deg and theta is a variable.
    If plane == "xz", 
        o-ray polarization // y-axis, e-ray polarization in xz-plane, phi = 0 deg and theta is a variable.

    Usage
    ------
    # create an instance of LBO object for wave vector in 'xy' principal dielectric plane.
    lbo_xy = ndispers.media.crystals.LBO_xy()
    # Get a refractive index for e-ray as a function of wavelength (um) and phi angle
    lbo_xy.n(0.6, 0.23*pi, pol='e')
    # Note: theta_rad is fixed at 0.5*pi value for xy plane and the third argument is phi_rad.
    # Derivative dispersion quantities are also easily obtained.
    lbo_xy.GVD(0.6, 0.23*pi, pol='e')

    @author: Akihiko Shimura
    """

    __slots__ = ["_A_x", "_B_x", "_C_x", "_D_x",
                 "_A_y", "_B_y", "_C_y", "_D_y",
                 "_A_z", "_B_z", "_C_z", "_D_z",]
    
    def __init__(self):
        super().__init__()

        # for x-axis
        self._A_x = 2.4542
        self._B_x = 0.01125
        self._C_x = 0.01135
        self._D_x = 0.01388
        # for y-axis
        self._A_y = 2.5390
        self._B_y = 0.01277
        self._C_y = 0.01189
        self._D_y = 0.01848
        # z-axis
        self._A_z = 2.5865
        self._B_z = 0.01310
        self._C_z = 0.01223
        self._D_z = 0.01861
    
    @property
    def property(self):
        msg =  ["A_x = %g" % self._A_x]
        msg += ["B_x = %g" % self._B_x]
        msg += ["C_x = %g" % self._C_x]
        msg += ["D_x = %g" % self._D_x]
        msg += ["A_y = %g" % self._A_y]
        msg += ["B_y = %g" % self._B_y]
        msg += ["C_y = %g" % self._C_y]
        msg += ["D_y = %g" % self._D_y]
        msg += ["A_z = %g" % self._A_z]
        msg += ["B_z = %g" % self._B_z]
        msg += ["C_z = %g" % self._C_z]
        msg += ["D_z = %g" % self._D_z]
        print("\n".join(msg))
    
    def n_x_expr(self):
        """ sympy expresssion, dispersion formula of x-axis (principal dielectric axis) """
        return sympy.sqrt(self._A_x + self._B_x/(wl**2 - self._C_x) - self._D_x * wl**2)
    
    def n_y_expr(self):
        """ sympy expresssion, dispersion formula of y-axis (principal dielectric axis) """
        return sympy.sqrt(self._A_y + self._B_y/(wl**2 - self._C_y) - self._D_y * wl**2)

    def n_z_expr(self):
        """ sympy expresssion, dispersion formula of z-axis (principal dielectric axis) """
        return sympy.sqrt(self._A_z + self._B_z/(wl**2 - self._C_z) - self._D_z * wl**2)


class LBO_xy(LBO):
    __slots__ = ["_LBO_xy__plane", "_LBO_xy__theta_rad", "_LBO_xy__phi_rad"]

    def __init__(self):
        super().__init__()
        # self.__doc__ = super().__doc__
        self._LBO_xy__plane = 'xy'
        self._LBO_xy__theta_rad = 0.5*pi
        self._LBO_xy__phi_rad = 'arb'
    
    @property
    def help(self):
        print(super().__doc__)
    
    @property
    def angles(self):
        msg =  ["plane = %s" % self._LBO_xy__plane]
        msg += ["theta_rad = %s" % self._LBO_xy__theta_rad]
        msg += ["phi_rad = %s" % self._LBO_xy__phi_rad]
        print("\n".join(msg))

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

    def n(self, wl_um, phi_rad, pol='o'):
        """
        Refractive index in xy plane.

        input
        ------
        wl_um     :  float, wavelength in um
        pol       :  str, 'o' or 'e', polarization of light
        phi_rad   :  float, 0 to 2pi radians
        (Note: theta_rad is fixed at 0.5*pi in xy principal plane.)

        return
        -------
        Refractive index, float
        """
        return super().n(wl_um, 0.5*pi, phi_rad, pol=pol)

    def dn_wl(self, wl_um, phi_rad, pol='o'):
        return super().dn_wl(wl_um, 0.5*pi, phi_rad, pol=pol)
    
    def d2n_wl(self, wl_um, phi_rad, pol='o'):
        return super().d2n_wl(wl_um, 0.5*pi, phi_rad, pol=pol)

    def d3n_wl(self, wl_um, phi_rad, pol='o'):
        return super().d3n_wl(wl_um, 0.5*pi, phi_rad, pol=pol)

    def GD(self, wl_um, phi_rad, pol='o'):
        """Group Delay [fs/mm]"""
        return super().GD(wl_um, 0.5*pi, phi_rad, pol=pol)
    
    def GV(self, wl_um, phi_rad, pol='o'):
        """Group Velocity [um/fs]"""
        return super().GV(wl_um, 0.5*pi, phi_rad, pol=pol)
    
    def ng(self, wl_um, phi_rad, pol='o'):
        """Group index, c/Group velocity"""
        return super().ng(wl_um, 0.5*pi, phi_rad, pol=pol)

    def GVD(self, wl_um, phi_rad, pol='o'):
        """Group Delay Dispersion [fs^2/mm]"""
        return super().GVD(wl_um, 0.5*pi, phi_rad, pol=pol)

    def TOD(self, wl_um, phi_rad, pol='o'):
        """Third Order Dispersion [fs^3/mm]"""
        return super().TOD(wl_um, 0.5*pi, phi_rad, pol=pol)


class LBO_yz(LBO):
    __slots__ = ["_LBO_yz__plane", "_LBO_yz__theta_rad", "_LBO_yz__phi_rad"]

    def __init__(self):
        super().__init__()
        # self.__doc__ = super().__doc__
        self._LBO_yz__plane = 'yz'
        self._LBO_yz__phi_rad = 0.5*pi
        self._LBO_yz__theta_rad = 'arb'
    
    @property
    def help(self):
        print(super().__doc__)
    
    @property
    def angles(self):
        msg =  ["plane = %s" % self._LBO_yz__plane]
        msg += ["theta_rad = %s" % self._LBO_yz__theta_rad]
        msg += ["phi_rad = %s" % self._LBO_yz__phi_rad]
        print("\n".join(msg))

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

    def n(self, wl_um, theta_rad, pol='o'):
        """
        Refractive index in yz plane.

        input
        ------
        wl_um     :  float, wavelength in um
        pol       :  str, 'o' or 'e', polarization of light
        theta_rad :  float, 0 to 2pi radians
        (Note: phi_rad is fixed at 0.5*pi in yz principal plane.)

        return
        -------
        Refractive index, float
        """
        return super().n(wl_um, theta_rad, 0.5*pi, pol=pol)

    def dn_wl(self, wl_um, theta_rad, pol='o'):
        return super().dn_wl(wl_um, theta_rad, 0.5*pi, pol=pol)
    
    def d2n_wl(self, wl_um, theta_rad, pol='o'):
        return super().d2n_wl(wl_um, theta_rad, 0.5*pi, pol=pol)

    def d3n_wl(self, wl_um, theta_rad, pol='o'):
        return super().d3n_wl(wl_um, theta_rad, 0.5*pi, pol=pol)

    def GD(self, wl_um, theta_rad, pol='o'):
        """Group Delay [fs/mm]"""
        return super().GD(wl_um, theta_rad, 0.5*pi, pol=pol)
    
    def GV(self, wl_um, theta_rad, pol='o'):
        """Group Velocity [um/fs]"""
        return super().GV(wl_um, theta_rad, 0.5*pi, pol=pol)
    
    def ng(self, wl_um, theta_rad, pol='o'):
        """Group index, c/Group velocity"""
        return super().ng(wl_um, theta_rad, 0.5*pi, pol=pol)

    def GVD(self, wl_um, theta_rad, pol='o'):
        """Group Delay Dispersion [fs^2/mm]"""
        return super().GVD(wl_um, theta_rad, 0.5*pi, pol=pol)

    def TOD(self, wl_um, theta_rad, pol='o'):
        """Third Order Dispersion [fs^3/mm]"""
        return super().TOD(wl_um, theta_rad, 0.5*pi, pol=pol)


class LBO_zx(LBO):
    __slots__ = ["_LBO_zx__plane", "_LBO_zx__theta_rad", "_LBO_zx__phi_rad"]

    def __init__(self):
        super().__init__()
        # self.__doc__ = super().__doc__
        self._LBO_zx__plane = 'zx'
        self._LBO_zx__theta_rad = 'arb'
        self._LBO_zx__phi_rad = 0.5*pi
    
    @property
    def help(self):
        print(super().__doc__)
    
    @property
    def angles(self):
        msg =  ["plane = %s" % self._LBO_zx__plane]
        msg += ["theta_rad = %s" % self._LBO_zx__theta_rad]
        msg += ["phi_rad = %s" % self._LBO_zx__phi_rad]
        print("\n".join(msg))

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

    def n(self, wl_um, theta_rad, pol='o'):
        """
        Refractive index in zx plane.

        input
        ------
        wl_um     :  float, wavelength in um
        pol       :  str, 'o' or 'e', polarization of light
        theta_rad :  float, 0 to 2pi radians
        (Note: phi_rad is fixed at 0.5*pi in zx principal plane.)

        return
        -------
        Refractive index, float
        """
        return super().n(wl_um, theta_rad, 0.5*pi, pol=pol)

    def dn_wl(self, wl_um, theta_rad, pol='o'):
        return super().dn_wl(wl_um, theta_rad, 0.5*pi, pol=pol)
    
    def d2n_wl(self, wl_um, theta_rad, pol='o'):
        return super().d2n_wl(wl_um, theta_rad, 0.5*pi, pol=pol)

    def d3n_wl(self, wl_um, theta_rad, pol='o'):
        return super().d3n_wl(wl_um, theta_rad, 0.5*pi, pol=pol)

    def GD(self, wl_um, theta_rad, pol='o'):
        """Group Delay [fs/mm]"""
        return super().GD(wl_um, theta_rad, 0.5*pi, pol=pol)
    
    def GV(self, wl_um, theta_rad, pol='o'):
        """Group Velocity [um/fs]"""
        return super().GV(wl_um, theta_rad, 0.5*pi, pol=pol)
    
    def ng(self, wl_um, theta_rad, pol='o'):
        """Group index, c/Group velocity"""
        return super().ng(wl_um, theta_rad, 0.5*pi, pol=pol)

    def GVD(self, wl_um, theta_rad, pol='o'):
        """Group Delay Dispersion [fs^2/mm]"""
        return super().GVD(wl_um, theta_rad, 0.5*pi, pol=pol)

    def TOD(self, wl_um, theta_rad, pol='o'):
        """Third Order Dispersion [fs^3/mm]"""
        return super().TOD(wl_um, theta_rad, 0.5*pi, pol=pol)