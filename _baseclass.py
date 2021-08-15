"""
base class for medium object - _baseclass.py
"""

import sympy
from sympy.utilities import lambdify
wl = sympy.Symbol('wl')
phi = sympy.Symbol('phi')
theta = sympy.Symbol('theta')
from math import pi
c_ms = 2.99792458e8 #(m/s) speed of light in vacuum
c_umfs = c_ms * 1e-9  #(um/fs)
import numpy

class Medium:
    """
    Medium object base class. Vacuum with n=1.0.

    @author: Akihiko Shimura
    """
    __slots__ = ["__plane", 
                 "_cached_n_func", "_cached_dn_wl_func", "_cached_d2n_wl_func", "_cached_d3n_wl_func",
                 "_cached_GD_func", "_cached_ng_func", "_cached_ng_func", "_cached_GV_func", "_cached_GVD_func", "_cached_TOD_func"]
    def __init__(self):
        self.__plane = 'arb'
        self._cached_n_func = {'o': 0, 'e': 0}
        self._cached_dn_wl_func = {'o': 0, 'e': 0}
        self._cached_d2n_wl_func = {'o': 0, 'e': 0}
        self._cached_d3n_wl_func = {'o': 0, 'e': 0}
        self._cached_GD_func = {'o': 0, 'e': 0}
        self._cached_ng_func = {'o': 0, 'e': 0}
        self._cached_GV_func = {'o': 0, 'e': 0}
        self._cached_GVD_func = {'o': 0, 'e': 0}
        self._cached_TOD_func = {'o': 0, 'e': 0}
    
    def clear(self):
        """ clear cached """
        self.__init__()

    @property
    def help(self):
        print(self.__doc__)

    def n_expr(self, pol):
        return 1.0

    """ Derivative expressions """
    def dn_wl_expr(self, pol):
        """ Sympy expression for first derivative of n with respect to wl """
        return sympy.diff(self.n_expr(pol), wl)
    
    def d2n_wl_expr(self, pol):
        """ Sympy expression for second derivative of n with respect to wl """
        return sympy.diff(self.dn_wl_expr(pol), wl)

    def d3n_wl_expr(self, pol):
        """ Sympy expression for third derivative of n with respect to wl """
        return sympy.diff(self.d2n_wl_expr(pol), wl)

    def GD_expr(self, pol):
        """ Sympy expression for group delay """
        return (self.n_expr(pol) - wl * self.dn_expr(pol)) * 1e3 / c_umfs
    
    def GV_expr(self, pol):
        """ Sympy expression for group velocity """
        return (c_umfs/self.n_expr(pol)) / (1 - (wl/self.n_expr(pol)) * self.dn_wl_expr(pol))
    
    def ng_expr(self, pol):
        """ Sympy expression for group index """
        n_expr = self.n_expr(pol)
        return n_expr * (1 - wl/n_expr * self.dn_wl_expr(pol))
    
    def GVD_expr(self, pol):
        """ Sympy expression for Group Delay Dispersion """
        return wl**3/(2*pi*c_umfs**2) * self.d2n_wl_expr(pol) * 1e3
    
    def TOD_expr(self, pol):
        """ Sympy expression for Third Order Dispersion """
        return - wl**4/(4*pi**2*c_umfs**3) * (3*self.d2n_wl_expr(pol) + wl * self.d3n_wl_expr(pol)) * 1e3

    """ Functions """
    def n(self, wl_um, pol, theta_rad, phi_rad):
        """ n, refractive index of a medium """
        if self._cached_n_func[pol]:
            return self._cached_n_func[pol](wl_um, theta_rad, phi_rad)
        else:
            self._cached_n_func[pol] = lambdify([wl, theta, phi], self.n_expr(pol), 'numpy')
            return self._cached_n_func[pol](wl_um, theta_rad, phi_rad)
    
    def dn_wl(self, wl_um, pol, theta_rad, phi_rad):
        if self._cached_dn_wl_func[pol]:
            return self._cached_dn_wl_func[pol](wl_um, theta_rad, phi_rad)
        else:
            self._cached_dn_wl_func[pol] = lambdify([wl, theta, phi], self.dn_wl_expr(pol), 'numpy')
        return self._cached_dn_wl_func[pol](wl_um, theta_rad, phi_rad)

    def d2n_wl(self, wl_um, pol, theta_rad, phi_rad):
        if self._cached_d2n_wl_func[pol]:
            return self._cached_d2n_wl_func[pol](wl_um, theta_rad, phi_rad)
        else:
            self._cached_d2n_wl_func[pol] = lambdify([wl, theta, phi], self.d2n_wl_expr(pol), 'numpy')
        return self._cached_d2n_wl_func[pol](wl_um, theta_rad, phi_rad)

    def d3n_wl(self, wl_um, pol, theta_rad, phi_rad):
        if self._cached_d3n_wl_func[pol]:
            return self._cached_d3n_wl_func[pol](wl_um, theta_rad, phi_rad)
        else:
            self._cached_d3n_wl_func[pol] = lambdify([wl, theta, phi], self.d3n_wl_expr(pol), 'numpy')
        return self._cached_d3n_wl_func[pol](wl_um, theta_rad, phi_rad)

    def GD(self, wl_um, pol, theta_rad, phi_rad):
        """ Group Delay [fs/mm] """
        if self._cached_GD_func[pol]:
            return self._cached_GD_func[pol](wl_um, theta_rad, phi_rad)
        else:
            self._cached_GD_func[pol] = lambdify([wl, theta, phi], self.GD_expr(pol), 'numpy')
        return self._cached_GD_func[pol](wl_um, theta_rad, phi_rad)
    
    def GV(self, wl_um, pol, theta_rad, phi_rad):
        """ Group Velocity [um/fs] """
        if self._cached_GV_func[pol]:
            return self._cached_GV_func[pol](wl_um, theta_rad, phi_rad)
        else:
            self._cached_GV_func[pol] = lambdify([wl, theta, phi], self.GV_expr(pol), 'numpy')
        return self._cached_GV_func[pol](wl_um, theta_rad, phi_rad)
    
    def ng(self, wl_um, pol, theta_rad, phi_rad):
        """ Group index, c/Group velocity [] """
        if self._cached_ng_func[pol]:
            return self._cached_ng_func[pol](wl_um, theta_rad, phi_rad)
        else:
            self._cached_ng_func[pol] = lambdify([wl, theta, phi], self.ng_expr(pol), 'numpy')
        return self._cached_ng_func[pol](wl_um, theta_rad, phi_rad)

    def GVD(self, wl_um, pol, theta_rad, phi_rad):
        """ Group Delay Dispersion [fs^2/mm] """
        if self._cached_GVD_func[pol]:
            return self._cached_GVD_func[pol](wl_um, theta_rad, phi_rad)
        else:
            self._cached_GVD_func[pol] = lambdify([wl, theta, phi], self.GVD_expr(pol), 'numpy')
        return self._cached_GVD_func[pol](wl_um, theta_rad, phi_rad)

    def TOD(self, wl_um, pol, theta_rad, phi_rad):
        """ Third Order Dispersion [fs^3/mm] """
        if self._cached_TOD_func[pol]:
            return self._cached_TOD_func[pol](wl_um, theta_rad, phi_rad)
        else:
            self._cached_TOD_func[pol] = lambdify([wl, theta, phi], self.TOD_expr(pol), 'numpy')
            return self._cached_TOD_func[pol](wl_um, theta_rad, phi_rad)