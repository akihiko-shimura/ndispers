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

import numpy as np
from collections import defaultdict

def _returnShape(*args):
    return np.broadcast(*args).shape

def _arg_signchange(a):
    a_sign = np.sign(a)
    if_signflip = ((np.roll(a_sign, 1) - a_sign) != 0).astype(int)
    if_signflip[0] = 0
    arg_signflip = np.where(if_signflip == 1)
    return arg_signflip

class Medium:
    """
    Medium object base class. Vacuum with n=1.0.

    @author: Akihiko Shimura
    """
    __slots__ = ["__plane", "__theta_rad", "__phi_rad", "_cached_func_dict"]


    def __init__(self):
        self.__plane = 'arb'
        self.__theta_rad = 'arb'
        self.__phi_rad = 'arb'
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
    
    @property
    def plane(self):
        return self.__plane
    @property
    def theta_rad(self):
        return self.__theta_rad
    @property
    def phi_rad(self):
        return self.__phi_rad

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
        array_args = map(np.asarray, args)
        func = self._cached_func_dict[expr.__name__][pol]
        if func:
            return np.resize(func(*args), _returnShape(*array_args))
        else:
            func = lambdify([wl, theta, phi], expr(pol), 'numpy')
            self._cached_func_dict[expr.__name__][pol] = func
            return np.resize(func(*args), _returnShape(*array_args))
    
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

    """ Methods for nonlinear optics """
    def pmAngles_sfg(self, wl1, wl2, tol_deg=0.01, deg=False):
        """
        Phase-matching angles in sum-frequency generation (SFG).

        Parameters
        ----------
        wl1 : float, 
            pump wavelength 1 in um.
        wl2 : float, 
            pump wavelength 2 in um.
        tol_deg : float, optional
            Tolerance error of angle in degree. Defaults to 0.01 deg.
        deg : bool, optional
            If returned angle is expressed in radians (False) or degrees (True).

        Return
        ------
        d : dict,
            Phase-matching angles for types of interactions.
        """
        wl3 = 1./(1./wl1 + 1./wl2) # SFG energy conversion

        def pmAngle_for_pol(pol1, pol2, pol3):
            angle_ar = np.arange(0, 90, tol_deg) * pi/180 #deg to rad, accuracy 0.01 deg
            n1 = self.n(wl1, angle_ar, pol=pol1)
            n2 = self.n(wl2, angle_ar, pol=pol2)
            n3 = self.n(wl3, angle_ar, pol=pol3)
            dk = 2*pi * (n3/wl3 - n2/wl2 - n1/wl1)
            angle_pm = angle_ar[_arg_signchange(dk)]
            if deg:
                angle_pm *= 180/pi

            pm_angles = dict()
            if self.theta_rad == 'var':
                pm_angles['theta'] = angle_pm
                pm_angles['phi'] = None
            elif self.phi_rad == 'var':
                pm_angles['phi'] = angle_pm
                pm_angles['theta'] = None
            return pm_angles
        
        d = dict()
        # Type-I interactions
        d['ooe'] = pmAngle_for_pol('o', 'o', 'e') #negative
        d['eeo'] = pmAngle_for_pol('e', 'e', 'o') #positive
        # Type-II interactions
        d['oee'] = pmAngle_for_pol('o', 'e', 'e') #nega1
        d['eoe'] = pmAngle_for_pol('e', 'o', 'e') #nega2
        d['eoo'] = pmAngle_for_pol('e', 'o', 'o') #posi1
        d['oeo'] = pmAngle_for_pol('o', 'e', 'o') #posi2
        return d
