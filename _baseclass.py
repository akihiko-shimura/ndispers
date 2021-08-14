"""
base class for medium object - _baseclass.py
"""

import sympy
from sympy.utilities import lambdify
wl = sympy.Symbol('wl')
phi = sympy.Symbol('phi')
theta = sympy.Symbol('theta')
import physcon
pi = physcon.pi
c = physcon.c * 1e-9  # um/fs, speed of light in vacuum

from functools import cached_property

class Medium:
    """
    Medium object base class. Vacuum with n=1.0.

    @author: Akihiko Shimura
    """
    __slots__ = ["__plane"]
    def __init__(self):
        self.__plane = 'arb'

    @property
    def help(self):
        print(self.__doc__)

    def n_expr(self, pol):
        return 1.0

    """ Derivative expressions """
    @cached_property
    def dn_expr(self, pol):
        """ Sympy expression for first derivative of n with respect to wl """
        return sympy.diff(self.n_expr(pol=pol), wl)

    @cached_property
    def dn2_expr(self, pol):
        """ Sympy expression for second derivative of n with respect to wl """
        return sympy.diff(self.dn_expr(pol=pol), wl)

    @cached_property
    def dn3_expr(self, pol):
        """ Sympy expression for third derivative of n with respect to wl """
        return sympy.diff(self.dn2_expr(pol=pol), wl)

    @cached_property
    def GD_expr(self, pol):
        """ Sympy expression for group delay """
        return (self.n_expr(pol=pol) - wl * self.dn_expr(pol=pol)) * 1e3 / c
    
    @cached_property
    def GV_expr(self, pol):
        """ Sympy expression for group velocity """
        return (c/self.n_expr(pol=pol)) / (1 - (wl/self.n_expr(pol=pol)) * self.dn_expr(pol=pol))
    
    @cached_property
    def ng_expr(self, pol):
        """ Sympy expression for group index """
        n_expr = self.n_expr(pol=pol)
        return n_expr * (1 - wl/n_expr * self.dn_expr(pol=pol))
    
    @cached_property
    def GVD_expr(self, pol):
        """ Sympy expression for Group Delay Dispersion """
        return wl**3/(2*pi*c**2) * self.dn2_expr(pol=pol) * 1e3
    
    @cached_property
    def TOD_expr(self, pol):
        """ Sympy expression for Third Order Dispersion """
        return - wl**4/(4*pi**2*c**3) * (3*self.dn2_expr(pol=pol) + wl * self.dn3_expr(pol=pol)) * 1e3

    """ Functions """
    def n(self, wl_um, pol, theta_rad, phi_rad):
        return lambdify([wl, theta, phi], self.n_expr(pol=pol), 'numpy')(wl_um, theta_rad, phi_rad)

    def dn(self, wl_um, pol, theta_rad, phi_rad):
        return lambdify([wl, theta, phi], self.dn_expr(pol=pol), 'numpy')(wl_um, theta_rad, phi_rad)

    def dn2(self, wl_um, pol, theta_rad, phi_rad):
        return lambdify([wl, theta, phi], self.dn2_expr(pol=pol), 'numpy')(wl_um, theta_rad, phi_rad)

    def dn3(self, wl_um, pol, theta_rad, phi_rad):
        return lambdify([wl, theta, phi], self.dn3_expr(pol=pol), 'numpy')(wl_um, theta_rad, phi_rad)

    def GD(self, wl_um, pol, theta_rad, phi_rad):
        """ Group Delay [fs/mm] """
        return lambdify([wl, theta, phi], self.GD_expr(pol=pol), 'numpy')(wl_um, theta_rad, phi_rad)
    
    def GV(self, wl_um, pol, theta_rad, phi_rad):
        """ Group Velocity [um/fs] """
        return lambdify([wl, theta, phi], self.GV_expr(pol=pol), 'numpy')(wl_um, theta_rad, phi_rad)
    
    def ng(self, wl_um, pol, theta_rad, phi_rad):
        """ Group index, c/Group velocity """
        return lambdify([wl, theta, phi], self.ng_expr(pol=pol), 'numpy')(wl_um, theta_rad, phi_rad)

    def GVD(self, wl_um, pol, theta_rad, phi_rad):
        """ Group Delay Dispersion [fs^2/mm] """
        return lambdify([wl, theta, phi], self.GVD_expr(pol=pol), 'numpy')(wl_um, theta_rad, phi_rad)

    def TOD(self, wl_um, pol, theta_rad, phi_rad):
        """ Third Order Dispersion [fs^3/mm] """
        return lambdify([wl, theta, phi], self.TOD_expr(pol=pol), 'numpy')(wl_um, theta_rad, phi_rad)