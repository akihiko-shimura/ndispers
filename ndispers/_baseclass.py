"""
base class for medium object - _baseclass.py
"""
import sympy
from sympy.utilities import lambdify

wl = sympy.Symbol('lambda')
phi = sympy.Symbol('phi')
theta = sympy.Symbol('theta')
T = sympy.Symbol('T')

from math import pi

c_ms = 2.99792458e8 #(m/s) speed of light in vacuum
c_umfs = c_ms * 1e-9  #(µm/fs)

import numpy as np

from scipy.optimize import brentq

from .helper import returnShape, vars2


class Medium:
    """
    Medium object base class. Vacuum with n=1.0.
    """
    __slots__ = ["_plane", "_theta_rad", "_phi_rad"]

    # Lambdified functions, keyed (class, expression name, pol). Class-level:
    # a medium's constants are fixed in __init__ and its expressions depend only
    # on the class, so two instances of the same class share every function.
    # This is what makes the second instance's first call cheap.
    _lambdified = {}

    # what pol means when the caller does not say; SLN overrides with 'e'
    # because its Sellmeier equation exists only for the extraordinary ray
    _default_pol = 'o'

    def __init__(self):
        self._plane = 'arb'
        self._theta_rad = 'arb'
        self._phi_rad = 'arb'

    def clear(self):
        """clear this medium's cached functions"""
        cls = type(self)
        for key in [k for k in Medium._lambdified if k[0] is cls]:
            del Medium._lambdified[key]

    def __reduce__(self):
        # State is fixed in __init__ and never mutated, so calling the class
        # again restores an equivalent instance.
        # ponytail: assumes no-arg __init__ (true for every medium in v0.5.2);
        # pass constructor args here if any medium ever takes them.
        return (self.__class__, ())


    @property
    def plane(self):
        return self._plane

    @property
    def theta_rad(self):
        return self._theta_rad

    @property
    def phi_rad(self):
        return self._phi_rad

    @property
    def symbols(self):
        return [wl, theta, phi, T]

    @property
    def constants(self):
        """dispersion-formula constants of this medium, with their sources' names"""
        return {k: v for k, v in vars2(self).items()
                if k not in ("_plane", "_theta_rad", "_phi_rad")}

    def __repr__(self):
        return f"{self.__class__}\n  plane: {self.plane}\n  theta_rad: {self.theta_rad}\n  phi_rad: {self.phi_rad}"

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
        return (self.n_expr(pol) - wl * self.dn_wl_expr(pol)) * 1e3 / c_umfs
    
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

    def woa_theta_expr(self, pol):
        """ Sympy expression for polar walkoff angle """
        return sympy.atan(- sympy.diff(self.n_expr(pol), theta) / self.n_expr(pol))
    
    def woa_phi_expr(self, pol):
        """ Sympy expression for azimuthal walkoff angle """
        return sympy.atan(- sympy.diff(self.n_expr(pol), phi) / self.n_expr(pol))
    
    def dndT_expr(self, pol):
        """ Sympy expression for dn/dT """
        return sympy.diff(self.n_expr(pol), T)

    def dndT2_expr(self, pol):
        """ Sympy expression for d^2n/dT^2 """
        return sympy.diff(self.dndT_expr(pol), T)

    """ lambdified functions """
    def _full_args(self, args):
        """
        Fill in the angle this medium holds fixed.

        Callers pass one value per *varying* symbol - (wl, angle, T) for a
        crystal, (wl, T) for a glass - and the fixed angle is inserted here.
        A call that already carries one value per symbol passes through, so
        internal code like dk_sfg may spell out all four.
        """
        if len(args) == len(self.symbols):
            return args
        if self.theta_rad == 'var':
            phi_rad = 0 if isinstance(self.phi_rad, str) else self.phi_rad
            return (args[0], args[1], phi_rad) + tuple(args[2:])
        if self.phi_rad == 'var':
            return (args[0], self.theta_rad, args[1]) + tuple(args[2:])
        return args

    def _func(self, expr, *args, pol=None):
        if pol is None:
            pol = self._default_pol
        args = self._full_args(args)
        array_args = map(np.asarray, args)
        key = (type(self), expr.__name__, pol)
        func = Medium._lambdified.get(key)
        if func is None:
            func = lambdify(self.symbols, expr(pol), 'numpy')
            Medium._lambdified[key] = func
        out = np.resize(func(*args), returnShape(*array_args))
        # scalar in, scalar out - not a 0-d ndarray
        return out.item() if out.ndim == 0 else out
    
    def n(self, *args, pol=None):
        """
        Refractive index for one eigen polarization of light.

        Parameters
        ----------
        wl_um : float or array_like
            Wavelength in µm.
        angle_rad : float or array_like
            Angle of the wavevector, 0 to pi radians. Which angle this is
            depends on the medium: it is theta for a uniaxial crystal and for
            a principal plane that fixes phi, and phi for a plane that fixes
            theta. Read ``theta_rad`` and ``phi_rad`` to see which one varies;
            the fixed one is supplied automatically.
        T_degC : float or array_like
            Temperature of the medium in degree C. Media whose Sellmeier
            equation carries no temperature term still take this argument and
            ignore it, so that every medium has the same signature. Their
            ``dndT`` is zero.
        pol : {'o', 'e'}, optional
            Polarization of light.

        Returns
        -------
        float or numpy.ndarray
            Refractive index. Scalar input returns a float.

        Notes
        -----
        The other dispersion methods - ``dn_wl``, ``GD``, ``GV``, ``ng``,
        ``GVD``, ``TOD``, ``woa_theta``, ``woa_phi``, ``dndT`` and ``dndT2`` -
        take the same arguments and differ only in what they return. Glasses
        take ``(wl_um, T_degC)`` and no polarization.
        """
        return self._func(self.n_expr, *args, pol=pol)
    
    def dn_wl(self, *args, pol=None):
        return self._func(self.dn_wl_expr, *args, pol=pol)
    
    def d2n_wl(self, *args, pol=None):
        return self._func(self.d2n_wl_expr, *args, pol=pol)

    def d3n_wl(self, *args, pol=None):
        return self._func(self.d3n_wl_expr, *args, pol=pol)

    def GD(self, *args, pol=None):
        return self._func(self.GD_expr, *args, pol=pol)
    
    def GV(self, *args, pol=None):
        return self._func(self.GV_expr, *args, pol=pol)
    
    def ng(self, *args, pol=None):
        return self._func(self.ng_expr, *args, pol=pol)
    
    def GVD(self, *args, pol=None):
        return self._func(self.GVD_expr, *args, pol=pol)
    
    def TOD(self, *args, pol=None):
        return self._func(self.TOD_expr, *args, pol=pol)

    def woa_theta(self, *args, pol='e'):
        """ Polar walk-off angle (rad) """
        return self._func(self.woa_theta_expr, *args, pol=pol)
    
    def woa_phi(self, *args, pol='e'):
        """ Azimuthal walk-off angle (rad) """
        return self._func(self.woa_phi_expr, *args, pol=pol)

    def dndT(self, *args, pol=None):
        """
        dn/dT function for given arguments (angle, temperature and polarization)

        NOTE
        ----
        Here, self.dndT_expr is given by sympy.diff(self.n_expr(pol)), so there is no need to give dndT_expr explicitly.
        """
        return self._func(self.dndT_expr, *args, pol=pol)
    
    def dndT2(self, *args, pol=None):
        """ d^2n/dT^2 """
        return self._func(self.dndT2_expr, *args, pol=pol)

    
    """ 
    Methods for three-wave interactions
    """
    def dk_sfg(self, wl1, wl2, angle_rad, T_degC, pol1, pol2, pol3):
        """
        Wavevector mismatch for sum-frequency generation (SFG).

        Parameters
        ----------
        wl1 : float or array_like
            1st pump wavelength in µm.
        wl2 : float or array_like
            2nd pump wavelength in µm.
        angle_rad : float or array_like
            theta or phi angles in radians.
        T_degC  :  float or array_like
            Crystal temperature in degC.
        pol1: {'o', 'e'}
            Polarization of 1st pump wave.
        pol2: {'o', 'e'}
            Polarization of 2nd pump wave.
        pol3: {'o', 'e'}
            Polarization of sum-frequency wave.

        Return
        ------
        float or array_like
            Wavevector mismatch for SFG (in rad/µm)
        """
        wl3 = 1./(1./wl1 + 1./wl2)
        n1 = self.n(wl1, angle_rad, T_degC, pol=pol1)
        n2 = self.n(wl2, angle_rad, T_degC, pol=pol2)
        n3 = self.n(wl3, angle_rad, T_degC, pol=pol3)
        dk_sfg = 2*pi * (n3/wl3 - n2/wl2 - n1/wl1)
        return dk_sfg

    def qpm_period_sfg(self, wl1, wl2, angle_rad, T_degC, pol1, pol2, pol3, order=1):
        """
        Quasi-phase-matching period for sum-frequency generation.

        Lambda = 2 pi m / |dk| in µm, with m the QPM order and dk the collinear
        wavevector mismatch of ``dk_sfg``. For the usual PPLN/PPKTP geometry -
        propagation along a principal axis, all three waves polarized along the
        polar axis (d33) - use the extraordinary ray at angle_rad = pi/2 and
        ``('e', 'e', 'e')``; the temperature tuning follows from the medium's
        dn/dT.

        Parameters
        ----------
        wl1, wl2 : float or array_like
            Pump wavelengths in µm.
        angle_rad : float or array_like
            theta or phi angle in radians (see ``n``).
        T_degC : float or array_like
            Crystal temperature in degC.
        pol1, pol2, pol3 : {'o', 'e'}
            Polarizations of the two pumps and of the sum-frequency wave.
        order : int, default 1
            QPM order m.

        Return
        ------
        float or array_like
            Poling period in µm. Infinite where the interaction is
            birefringently phase-matched (dk = 0).
        """
        dk = self.dk_sfg(wl1, wl2, angle_rad, T_degC, pol1, pol2, pol3)
        return 2 * pi * order / np.abs(dk)

    def pmAngles_sfg(self, wl1, wl2, T_degC, tol_deg=0.001, deg=False):
        """
        Phase-matching (PM) angles for sum-frequency generation (SFG) and sum-frequency wavelength.

        Parameters
        ----------
        wl1 : float or array_like
            1st pump wavelength in µm.
        wl2 : float or array_like
            2nd pump wavelength in µm.
        tol_deg : float, default=0.001
            Absolute tolerance of the returned angle in degrees. An endpoint of
            the 0 to 90 degree range counts as a solution when its wavevector
            mismatch is within this tolerance of zero, so noncritical phase
            matching is reported as 90 (or 0) degrees rather than as no
            solution at all.
        deg : bool, default=False
            If returned angles are expressed in radians (False) or degrees (True).

        Return
        ------
        dict,
            'wl3'  :  wavelength of SFG
            {'ooe'}  :  PM angle for negative type-I
            {'eeo'}  :  PM angle for positive type-I
            {'oee', 'eoe'}  :  PM angles for negative type-II
            {'oee', 'eoe'}  :  PM angles for positive type-II
        """
        wl3 = 1./(1./wl1 + 1./wl2)

        def pmAngle_for_pol(pol1, pol2, pol3):
            # coarse grid locates each sign change, brentq refines it
            angle_ar = np.linspace(0, 0.5*pi, 361)
            try:
                dk_ar = self.dk_sfg(wl1, wl2, angle_ar, T_degC, pol1, pol2, pol3)
            except ValueError:
                # this medium has no Sellmeier equation for one of these
                # polarizations (SLN is e-ray only), so the combination has no
                # solution rather than being an error
                empty = {'theta': [], 'phi': []}
                empty['theta' if self.theta_rad == 'var' else 'phi'] = []
                return empty
            crossings = np.nonzero(np.diff(np.signbit(dk_ar)))[0]
            dk = lambda a: float(self.dk_sfg(wl1, wl2, a, T_degC, pol1, pol2, pol3))
            angle_pm = [brentq(dk, angle_ar[i], angle_ar[i+1], xtol=tol_deg*pi/180)
                        for i in crossings]

            # Noncritical phase matching puts the root exactly on an endpoint of
            # [0, pi/2], where a sign change never happens, so the scan above
            # misses it: LBO at its 90-degree temperature, the shortest SHG
            # wavelength of a crystal, and every NCPM cut would report "no
            # solution". Accept an endpoint whose dk is within one angular
            # tolerance of zero, judged by the local slope.
            tol_rad = tol_deg * pi / 180
            for i_edge, edge in ((0, 0.0), (-1, 0.5 * pi)):
                if any(abs(a - edge) <= tol_rad for a in angle_pm):
                    continue
                dk_edge = float(dk_ar[i_edge])
                if dk_edge == 0.0:
                    angle_pm.append(edge)
                    continue
                i_in = 1 if i_edge == 0 else -2
                slope = abs(dk_ar[i_in] - dk_edge) / (angle_ar[1] - angle_ar[0])
                if slope > 0 and abs(dk_edge) <= slope * tol_rad:
                    angle_pm.append(edge)
            angle_pm.sort()
            if deg:
                angle_pm = [a * 180/pi for a in angle_pm]
            pm_angles = dict()
            if self.theta_rad == 'var':
                pm_angles['theta'] = angle_pm
                pm_angles['phi'] = None
            elif self.phi_rad == 'var':
                pm_angles['phi'] = angle_pm
                pm_angles['theta'] = None
            return pm_angles
        
        d = dict()
        d['wl3'] = wl3
        # Type-I interaction
        d['ooe'] = pmAngle_for_pol('o', 'o', 'e') #negative
        d['eeo'] = pmAngle_for_pol('e', 'e', 'o') #positive
        # Type-II interaction
        d['oee'] = pmAngle_for_pol('o', 'e', 'e') #nega1
        d['eoe'] = pmAngle_for_pol('e', 'o', 'e') #nega2
        d['eoo'] = pmAngle_for_pol('e', 'o', 'o') #posi1
        d['oeo'] = pmAngle_for_pol('o', 'e', 'o') #posi2
        return d

    def pmFactor_sfg(self, wl1, wl2, angle_rad, T_degC, pol1, pol2, pol3, L_mm):
        """
        Phase-matching factor, sin^2((0.5*dk*L)/(0.5*dk*L)), for sum-frequency generation (SFG).

        Parameters
        ----------
        wl1 : float or array_like
            1st pump wavelength in µm.
        wl2 : float or array_like
            2nd pump wavelength in µm.
        angle_rad : float or array_like
            theta or phi angles in radians.
        T_degC  :  float or array_like
            Crystal temperature in degC.
        pol1: {'o', 'e'}
            Polarization of 1st pump wave.
        pol2: {'o', 'e'}
            Polarization of 2nd pump wave.
        pol3: {'o', 'e'}
            Polarization of sum-frequency wave.
        L_mm : float
            Crystal length in mm.

        Return
        ------
        float or array_like
            Phase-matching factor for SFG.
        """
        L_um = L_mm * 1e3
        t = 0.5 * self.dk_sfg(wl1, wl2, angle_rad, T_degC, pol1, pol2, pol3) * L_um
        return (np.sin(t)/t)**2
    
 