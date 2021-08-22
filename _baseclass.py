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
from collections import defaultdict

class Medium:
    """
    Medium object base class. Vacuum with n=1.0.

    @author: Akihiko Shimura
    """
    __slots__ = ["__plane", "_cached_func_dict"]

    def __init__(self):
        self.__plane = 'arb'
        self._cached_func_dict = defaultdict(dict)
        self._cached_func_dict['n_expr'] = {'o': 0, 'e': 0}
        self._cached_func_dict['dn_wl_expr'] = {'o': 0, 'e': 0}
        self._cached_func_dict['d2n_wl_expr'] = {'o': 0, 'e': 0}
        self._cached_func_dict['d3n_wl_expr'] = {'o': 0, 'e': 0}
        self._cached_func_dict['GD_expr'] = {'o': 0, 'e': 0}
        self._cached_func_dict['GV_expr'] = {'o': 0, 'e': 0}
        self._cached_func_dict['ng_expr'] = {'o': 0, 'e': 0}
        self._cached_func_dict['GD_expr'] = {'o': 0, 'e': 0}
        self._cached_func_dict['GVD_expr'] = {'o': 0, 'e': 0}
        self._cached_func_dict['TOD_expr'] = {'o': 0, 'e': 0}
    
    def clear(self):
        """ clear cached functions """
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

    """ lambdified functions """
    def _func(self, expr, *args, pol='o'):
        array_args = map(numpy.asarray, args)
        func = self._cached_func_dict[expr.__name__][pol]
        if func:
            return func(*array_args)
        else:
            func = lambdify([wl, theta, phi], expr(pol), 'numpy')
            self._cached_func_dict[expr.__name__][pol] = func
            return func(*array_args)
    
    def n(self, *args, pol='o'):
        return self._func(self.n_expr, *args, pol=pol)
    
    def dn_wl(self, *args, pol='o'):
        return self._func(self.dn_wl_expr, *args, pol=pol)
    
    def d2n_wl(self, *args, pol='o'):
        return self._func(self.d2n_wl_expr, *args, pol=pol)

    def d3n_wl(self, *args, pol='o'):
        return self._func(self.d3n_wl_expr, *args, pol=pol)

    def GD(self, *args, pol='o'):
        return self._func(self.GD_expr, *args, pol=pol)
    
    def GV(self, *args, pol='o'):
        return self._func(self.GV_expr, *args, pol=pol)
    
    def ng(self, *args, pol='o'):
        return self._func(self.ng_expr, *args, pol=pol)
    
    def GVD(self, *args, pol='o'):
        return self._func(self.GVD_expr, *args, pol=pol)
    
    def TOD(self, *args, pol='o'):
        return self._func(self.TOD_expr, *args, pol=pol)
