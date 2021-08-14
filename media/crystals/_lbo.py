import sympy
from sympy.utilities import lambdify
from ndispers._baseclass import Medium, wl, phi, theta, pi
from math import isclose
from functools import cached_property

class LBO(Medium):
    """
    LBO (Li B_3 O_5) crystal

    - Point group : mm2
    - Crystal class : orthorhombic 
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
    lbo = ndispers.media.crystals.LBO(plane='xy')
    # Get a refractive index for e-ray as a function of wavelength (um) and theta angle
    lbo.n(0.6, pol='e', theta_rad=0.5*pi, phi_rad=0.2*pi)
    # Note: theta_rad is actually fixed at 0.5*pi value for xy plane. So if input any values, the output is the same.
    # Derivative dispersion quantities are also easily obtained.
    lbo.GVD(0.6, pol='e', theta_rad=0.5*pi, phi_rad=0.2*pi)

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
        msg =  ["A_x = %.4f" % self._A_x]
        msg += ["B_x = %.5f" % self._B_x]
        msg += ["C_x = %.5f" % self._C_x]
        msg += ["D_x = %.5f" % self._D_x]
        msg += ["A_y = %.4f" % self._A_y]
        msg += ["B_y = %.5f" % self._B_y]
        msg += ["C_y = %.5f" % self._C_y]
        msg += ["D_y = %.5f" % self._D_y]
        msg += ["A_z = %.4f" % self._A_z]
        msg += ["B_z = %.5f" % self._B_z]
        msg += ["C_z = %.5f" % self._C_z]
        msg += ["D_z = %.5f" % self._D_z]
        print("\n".join(msg))
    
    # sympy expressions for each case
    @cached_property
    def n_x_expr(self):
        """ sympy expresssion, dispersion formula of x-axis (principal dielectric axis) """
        return sympy.sqrt(self._A_x + self._B_x/(wl**2 - self._C_x) - self._D_x * wl**2)
    
    @cached_property
    def n_y_expr(self):
        """ sympy expresssion, dispersion formula of y-axis (principal dielectric axis) """
        return sympy.sqrt(self._A_y + self._B_y/(wl**2 - self._C_y) - self._D_y * wl**2)

    @cached_property
    def n_z_expr(self):
        """ sympy expresssion, dispersion formula of z-axis (principal dielectric axis) """
        return sympy.sqrt(self._A_z + self._B_z/(wl**2 - self._C_z) - self._D_z * wl**2)


class LBO_xy(LBO):
    __slots__ = ["_LBO_xy__plane", "_LBO_xy__theta_rad", "_LBO_xy__phi_rad"]

    def __init__(self):
        super().__init__()
        self._LBO_xy__plane = 'xy'
        self._LBO_xy__theta_rad = 0.5*pi
        self._LBO_xy__phi_rad = 'arb'
    
    @property
    def angles(self):
        msg =  ["plane = %s" % self._LBO_xy__plane]
        msg += ["theta_rad = %s" % self._LBO_xy__theta_rad]
        msg += ["phi_rad = %s" % self._LBO_xy__phi_rad]
        print("\n".join(msg))

    @cached_property
    def n_o_expr(self):
        """ sympy expresssion, 
        dispersion formula for o-ray polarization for a given principal plane """
        return super().n_z_expr()
    
    @cached_property
    def n_e_expr(self):
        """ sympy expresssion, 
        dispersion formula for e-ray polarization for a given principal plane """
        return super().n_x_expr() * super().n_y_expr() / sympy.sqrt( super().n_x_expr()**2 * sympy.cos(phi)**2 + super().n_y_expr()**2 * sympy.sin(phi)**2 )

    @cached_property
    def n_expr(self, pol):
        """ sympy expresssion, 
        dispersion formula for a given polarization """
        if pol == 'o':
            return self.n_o_expr()
        elif pol == 'e':
            return self.n_e_expr()
        else:
            raise ValueError("pol = '%s' must be 'o' or 'e'" % pol)

    def n(self, wl_um, pol, phi_rad):
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
        return super().n(wl_um, pol, 0.5*pi, phi_rad)

    def dn(self, wl_um, pol, phi_rad):
        return super().dn(wl_um, pol, 0.5*pi, phi_rad)
    
    def dn2(self, wl_um, pol, phi_rad):
        return super().dn2(wl_um, pol, 0.5*pi, phi_rad)

    def dn3(self, wl_um, pol, phi_rad):
        return super().dn3(wl_um, pol, 0.5*pi, phi_rad)

    def GD(self, wl_um, pol, phi_rad):
        """Group Delay [fs/mm]"""
        return super().GD(wl_um, pol, 0.5*pi, phi_rad)
    
    def GV(self, wl_um, pol, phi_rad):
        """Group Velocity [um/fs]"""
        return super().GV(wl_um, pol, 0.5*pi, phi_rad)
    
    def ng(self, wl_um, pol, phi_rad):
        """Group index, c/Group velocity"""
        return super().ng(wl_um, pol, 0.5*pi, phi_rad)

    def GVD(self, wl_um, pol, phi_rad):
        """Group Delay Dispersion [fs^2/mm]"""
        return super().GVD(wl_um, pol, 0.5*pi, phi_rad)

    def TOD(self, wl_um, pol, phi_rad):
        """Third Order Dispersion [fs^3/mm]"""
        return super().TOD(wl_um, pol, 0.5*pi, phi_rad)


class LBO_yz(LBO):
    __slots__ = ["_LBO_yz__plane", "_LBO_yz__theta_rad", "_LBO_yz__phi_rad"]

    def __init__(self):
        super().__init__()
        self._LBO_yz__plane = 'yz'
        self._LBO_yz__phi_rad = 0.5*pi
        self._LBO_yz__theta_rad = 'arb'
    
    @property
    def angles(self):
        msg =  ["plane = %s" % self._LBO_yz__plane]
        msg += ["theta_rad = %s" % self._LBO_yz__theta_rad]
        msg += ["phi_rad = %s" % self._LBO_yz__phi_rad]
        print("\n".join(msg))

    @cached_property
    def n_o_expr(self):
        """ sympy expresssion, 
        dispersion formula for o-ray polarization for yx principal plane """
        return super().n_x_expr()
    
    @cached_property
    def n_e_expr(self):
        """ sympy expresssion, 
        dispersion formula for e-ray polarization for yz principal plane """
        return super().n_y_expr() * super().n_z_expr() / sympy.sqrt( super().n_y_expr()**2 * sympy.cos(phi)**2 + super().n_z_expr()**2 * sympy.sin(phi)**2 )

    @cached_property
    def n_expr(self, pol):
        """ sympy expresssion, 
        dispersion formula for a given polarization """
        if pol == 'o':
            return self.n_o_expr()
        elif pol == 'e':
            return self.n_e_expr()
        else:
            raise ValueError("pol = '%s' must be 'o' or 'e'" % pol)

    def n(self, wl_um, pol, theta_rad):
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
        return super().n(wl_um, pol, theta_rad, 0.5*pi)

    def dn(self, wl_um, pol, theta_rad):
        return super().dn(wl_um, pol, theta_rad, 0.5*pi)
    
    def dn2(self, wl_um, pol, theta_rad):
        return super().dn2(wl_um, pol, theta_rad, 0.5*pi)

    def dn3(self, wl_um, pol, theta_rad):
        return super().dn3(wl_um, pol, theta_rad, 0.5*pi)

    def GD(self, wl_um, pol, theta_rad):
        """Group Delay [fs/mm]"""
        return super().GD(wl_um, pol, theta_rad, 0.5*pi)
    
    def GV(self, wl_um, pol, theta_rad):
        """Group Velocity [um/fs]"""
        return super().GV(wl_um, pol, theta_rad, 0.5*pi)
    
    def ng(self, wl_um, pol, theta_rad):
        """Group index, c/Group velocity"""
        return super().ng(wl_um, pol, theta_rad, 0.5*pi)

    def GVD(self, wl_um, pol, theta_rad):
        """Group Delay Dispersion [fs^2/mm]"""
        return super().GVD(wl_um, pol, theta_rad, 0.5*pi)

    def TOD(self, wl_um, pol, theta_rad):
        """Third Order Dispersion [fs^3/mm]"""
        return super().TOD(wl_um, pol, theta_rad, 0.5*pi)


class LBO_zx(LBO):
    __slots__ = ["_LBO_zx__plane", "_LBO_zx__theta_rad", "_LBO_zx__phi_rad"]

    def __init__(self):
        super().__init__()
        self._LBO_zx__plane = 'zx'
        self._LBO_zx__theta_rad = 'arb'
        self._LBO_zx__phi_rad = 0.5*pi
    
    @property
    def angles(self):
        msg =  ["plane = %s" % self._LBO_zx__plane]
        msg += ["theta_rad = %s" % self._LBO_zx__theta_rad]
        msg += ["phi_rad = %s" % self._LBO_zx__phi_rad]
        print("\n".join(msg))

    @cached_property
    def n_o_expr(self):
        """ sympy expresssion, 
        dispersion formula for o-ray polarization for zx principal plane """
        return super().n_y_expr()
    
    @cached_property
    def n_e_expr(self):
        """ sympy expresssion, 
        dispersion formula for e-ray polarization for zx principal plane """
        return super().n_z_expr() * super().n_z_expr() / sympy.sqrt( super().n_z_expr()**2 * sympy.cos(phi)**2 + super().n_x_expr()**2 * sympy.sin(phi)**2 )

    @cached_property
    def n_expr(self, pol):
        """ sympy expresssion, 
        dispersion formula for a given polarization """
        if pol == 'o':
            return self.n_o_expr()
        elif pol == 'e':
            return self.n_e_expr()
        else:
            raise ValueError("pol = '%s' must be 'o' or 'e'" % pol)

    def n(self, wl_um, pol, theta_rad):
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
        return super().n(wl_um, pol, theta_rad, 0.5*pi)

    def dn(self, wl_um, pol, theta_rad):
        return super().dn(wl_um, pol, theta_rad, 0.5*pi)
    
    def dn2(self, wl_um, pol, theta_rad):
        return super().dn2(wl_um, pol, theta_rad, 0.5*pi)

    def dn3(self, wl_um, pol, theta_rad):
        return super().dn3(wl_um, pol, theta_rad, 0.5*pi)

    def GD(self, wl_um, pol, theta_rad):
        """Group Delay [fs/mm]"""
        return super().GD(wl_um, pol, theta_rad, 0.5*pi)
    
    def GV(self, wl_um, pol, theta_rad):
        """Group Velocity [um/fs]"""
        return super().GV(wl_um, pol, theta_rad, 0.5*pi)
    
    def ng(self, wl_um, pol, theta_rad):
        """Group index, c/Group velocity"""
        return super().ng(wl_um, pol, theta_rad, 0.5*pi)

    def GVD(self, wl_um, pol, theta_rad):
        """Group Delay Dispersion [fs^2/mm]"""
        return super().GVD(wl_um, pol, theta_rad, 0.5*pi)

    def TOD(self, wl_um, pol, theta_rad):
        """Third Order Dispersion [fs^3/mm]"""
        return super().TOD(wl_um, pol, theta_rad, 0.5*pi)