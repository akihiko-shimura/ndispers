"""
base class for medium object - _baseclass.py
"""
from ._sym import sympy, wl, phi, theta, T

from math import pi

c_ms = 2.99792458e8 #(m/s) speed of light in vacuum
c_umfs = c_ms * 1e-9  #(µm/fs)

import numpy as np


from .helper import returnShape, vars2, brentq

import warnings


class TemperatureWarning(UserWarning):
    """This medium's thermo-optic coefficients are zero (none reported in the
    literature, or the source states none): T_degC is accepted but has no
    effect. Sweeping temperature returns constant results by design."""


class ValidityWarning(UserWarning):
    """A wavelength lies outside the medium's Sellmeier validity range.

    The returned value is an extrapolation of the fit, not a measured-range
    interpolation. Silence with warnings.simplefilter('ignore', ValidityWarning)
    if the extrapolation is deliberate."""


def _compiled_module(cls):
    """The pre-generated module for a medium class, or None."""
    import importlib
    try:
        mod = importlib.import_module(f"ndispers._compiled.{cls.__name__}")
    except ImportError:
        return None
    # a user subclass that happens to share a name must not pick up the
    # original's functions
    if mod.CLASS != f"{cls.__module__}.{cls.__qualname__}":
        return None
    return mod


def _compiled_func(cls, expr_name, pol):
    mod = _compiled_module(cls)
    return None if mod is None else mod.FUNCS.get((expr_name, pol))


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

    # (min, max) wavelength in µm of the Sellmeier validity, transcribed from
    # the class docstring (which stays authoritative - see it for per-pol
    # ranges and caveats). None = not stated in the source; no range check.
    _wl_range = None

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
        # One informative line: agents and users print media constantly.
        cls = type(self).__name__
        doc = (type(self).__doc__ or '').strip().splitlines()
        desc = doc[0].strip() if doc else ''
        bits = [desc] if desc else []
        if self.plane != 'arb':
            bits.append(f"{self.plane} plane")
        for name, val in (("theta", self.theta_rad), ("phi", self.phi_rad)):
            if val == 'var':
                bits.append(f"angle arg: {name}")
        if self._wl_range:
            bits.append(f"{self._wl_range[0]:g}-{self._wl_range[1]:g} um")
        return f"{cls}({', '.join(bits)})"

    def n_expr(self, pol):
        return 1.0

    def n_latex(self, pol=None):
        """LaTeX of the refractive-index formula for one polarization, as a string.

        Read from the pre-generated module when there is one, so that it does
        not need sympy; otherwise rendered from ``n_expr``."""
        if pol is None:
            pol = self._default_pol
        mod = _compiled_module(type(self))
        if mod is not None and pol in mod.LATEX:
            return mod.LATEX[pol]
        return sympy.latex(self.n_expr(pol))

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

    def _check_wl(self, wl_um):
        rng = self._wl_range
        if rng is None:
            return
        w = np.asarray(wl_um)
        # 1% grace band: the fit does not fail at the stated boundary, and
        # e.g. 1.064 um sits 0.4% above BetaBBO_Eimerl1987's 1.06 um limit -
        # warning on the most common use of BBO would be noise. The target
        # here is unit mistakes and gross extrapolation.
        if np.any((w < 0.99 * rng[0]) | (w > 1.01 * rng[1])):
            hint = (" Arguments are in um - did you pass nanometers?"
                    if np.any(w > 50) else "")
            warnings.warn(
                f"{type(self).__name__}: wavelength outside the Sellmeier "
                f"validity range {rng[0]:g}-{rng[1]:g} um; the result is an "
                f"extrapolation.{hint}", ValidityWarning, stacklevel=4)

    # classes already probed for a temperature term (once per process)
    _T_checked = set()

    def _warn_if_T_ignored(self):
        cls = type(self)
        if cls in Medium._T_checked:
            return
        Medium._T_checked.add(cls)      # before the probe: dndT re-enters _func
        lo, hi = self._wl_range or (0.5, 2.0)
        wls = np.linspace(lo * 1.01, hi * 0.99, 3)
        args = (wls, 25.0) if len(self.symbols) == 2 else (wls, 0.5, 25.0)
        try:
            if np.all(np.asarray(self.dndT(*args)) == 0):
                warnings.warn(
                    f"{cls.__name__}: thermo-optic coefficients are zero (none "
                    f"implemented - see docstring); T_degC is accepted but has "
                    f"no effect, so temperature sweeps return constant results.",
                    TemperatureWarning, stacklevel=5)
        except Exception:
            pass    # a probe must never break a computation

    def _signature(self):
        """The call signature of this medium's dispersion methods, for errors."""
        cls = type(self).__name__
        if len(self.symbols) == 2:
            return f"{cls} is isotropic: methods take (wl_um, T_degC) - no angle"
        which = 'phi_rad' if self.phi_rad == 'var' else 'theta_rad'
        return f"{cls} methods take (wl_um, {which}, T_degC)"

    def _func(self, expr, *args, pol=None):
        if pol is None:
            pol = self._default_pol
        if any(isinstance(a, str) for a in args):
            raise TypeError(
                "pol is keyword-only: write pol='o' or pol='e', not a "
                "positional argument. " + self._signature())
        # lists/tuples would crash inside the numpy formulas; convert here
        args = tuple(np.asarray(a, dtype=float) if isinstance(a, (list, tuple))
                     else a for a in args)
        if args:
            self._check_wl(args[0])
        self._warn_if_T_ignored()
        n_given = len(args)
        args = self._full_args(args)
        if len(args) != len(self.symbols):
            raise TypeError(
                f"{self._signature()}; got {n_given} positional arguments "
                f"(pol is keyword-only)")
        array_args = map(np.asarray, args)
        key = (type(self), expr.__name__, pol)
        func = Medium._lambdified.get(key)
        if func is None:
            func = _compiled_func(type(self), expr.__name__, pol)
            if func is None:
                # no pre-generated module for this class (a user subclass,
                # or the generator was not run): build it with sympy
                func = sympy.lambdify([s._sympy_() for s in self.symbols],
                                      expr(pol), 'numpy')
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

        Examples
        --------
        >>> import ndispers as nd
        >>> bbo = nd.media.crystals.BetaBBO_Eimerl1987()
        >>> bbo.n(0.532, 0, 25, pol='o')
        1.674884049110459
        >>> bbo.n(1.064, 0.3994, 25, pol='e')
        1.63644345403142
        >>> print(round(bbo.pmAngles_sfg(1.064, 1.064, 25, deg=True)['ooe']['theta'][0], 6))
        22.884169
        """
        return self._func(self.n_expr, *args, pol=pol)
    
    def dn_wl(self, *args, pol=None):
        """dn/dλ (1/µm). Same arguments as ``n``."""
        return self._func(self.dn_wl_expr, *args, pol=pol)
    
    def d2n_wl(self, *args, pol=None):
        """d²n/dλ² (1/µm²). Same arguments as ``n``."""
        return self._func(self.d2n_wl_expr, *args, pol=pol)

    def d3n_wl(self, *args, pol=None):
        """d³n/dλ³ (1/µm³). Same arguments as ``n``."""
        return self._func(self.d3n_wl_expr, *args, pol=pol)

    def GD(self, *args, pol=None):
        """Group delay per unit length, fs/mm. Same arguments as ``n``."""
        return self._func(self.GD_expr, *args, pol=pol)
    
    def GV(self, *args, pol=None):
        """Group velocity, µm/fs. Same arguments as ``n``."""
        return self._func(self.GV_expr, *args, pol=pol)
    
    def ng(self, *args, pol=None):
        """Group index (dimensionless). Same arguments as ``n``."""
        return self._func(self.ng_expr, *args, pol=pol)
    
    def GVD(self, *args, pol=None):
        """Group-velocity dispersion, fs²/mm. Positive = normal dispersion. Same arguments as ``n``."""
        return self._func(self.GVD_expr, *args, pol=pol)
    
    def TOD(self, *args, pol=None):
        """Third-order dispersion, fs³/mm. Same arguments as ``n``."""
        return self._func(self.TOD_expr, *args, pol=pol)

    def GVM(self, wl1, wl2, angle_rad, T_degC, pol1='o', pol2='o'):
        """
        Group-velocity mismatch between two waves, GD(wl1) - GD(wl2), in fs/mm.

        Positive means the wl1 wave is the slower one (arrives later). This is
        the walk-off per unit length that sets e.g. an OPA gain bandwidth
        (signal vs pump) or the temporal smearing of SFG (wl1 vs wl2). For a glass pass
        ``None`` as angle_rad (it is ignored).

        Examples
        --------
        >>> import ndispers as nd
        >>> bbo = nd.media.crystals.BetaBBO_Eimerl1987()
        >>> print(round(bbo.GVM(0.532, 1.064, 0.3994, 25, pol1='e', pol2='o'), 4))
        79.0863
        """
        args1, args2 = ((wl1, T_degC), (wl2, T_degC)) if len(self.symbols) == 2 \
            else ((wl1, angle_rad, T_degC), (wl2, angle_rad, T_degC))
        return self.GD(*args1, pol=pol1) - self.GD(*args2, pol=pol2)

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
            # plain floats: np.float64 in a printed result is noise
            angle_pm = [float(a) * (180/pi if deg else 1) for a in angle_pm]
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
        Phase-matching factor sinc²(Δk·L/2) = sin²(x)/x² with x = 0.5·Δk·L, for sum-frequency generation (SFG).

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

    # ------------------------------------------------------------------
    # Difference-frequency generation / optical parametric amplification and
    # oscillation.  The interaction is the same three-wave process as SFG -
    # k_p = k_s + k_i with 1/wl_p = 1/wl_s + 1/wl_i - so every quantity is the
    # SFG one evaluated at (wl_s, wl_i); what differs is which wavelengths are
    # given.  Indices (1, 2, 3) = (signal, idler, pump), the order SNLO uses
    # (red1, red2, blue), so 'ooe' is a Type I OPA (signal and idler
    # ordinary, pump extraordinary) and 'oee' / 'eoe' are Type II.
    # ------------------------------------------------------------------

    def acceptance_sfg(self, wl1, wl2, angle_rad, T_degC, pol1, pol2, pol3,
                       L_mm, param):
        """
        Acceptance width (FWHM of sinc²(Δk·L/2)) of SFG along one variable.

        The full width over which the phase-matching factor stays >= 0.5,
        holding every other variable fixed at the given values. Evaluate at a
        phase-matched point (from ``pmAngles_sfg``) for the usual meaning.

        Parameters
        ----------
        param : {'wl', 'wl1', 'wl2', 'theta', 'T'}
            The variable to widen. 'wl1'/'wl2': that input wave alone
            (mix acceptance; the sum wavelength follows energy conservation).
            'wl': the degenerate-SHG spectral acceptance - wl1 and wl2 (equal)
            are swept together as the fundamental. 'theta': the angle argument
            (internal). 'T': temperature.
        L_mm : float
            Crystal length in mm.

        Returns
        -------
        float
            FWHM in the variable's own unit: µm for wavelengths, rad for
            'theta', K for 'T'. Returns numpy.inf if the factor never falls
            below 0.5 within a wide search (noncritical direction).

        Examples
        --------
        >>> import ndispers as nd
        >>> from math import radians
        >>> bbo = nd.media.crystals.BetaBBO_Tamosauskas2018()
        >>> th = radians(bbo.pmAngles_sfg(0.8, 0.8, 25, deg=True)['ooe']['theta'][0])
        >>> dl = bbo.acceptance_sfg(0.8, 0.8, th, 25, 'o', 'o', 'e', 1.0, 'wl')
        >>> print(round(dl * 1e3, 2))   # µm -> nm; ~5 nm: why fs SHG needs thin BBO
        4.91
        """
        if param == 'wl' and not np.isclose(wl1, wl2):
            raise ValueError("param='wl' sweeps the degenerate-SHG fundamental "
                             "and needs wl1 == wl2; use 'wl1' or 'wl2'")
        L_um = L_mm * 1e3
        base = {'wl1': wl1, 'wl2': wl2, 'wl': wl1, 'theta': angle_rad, 'T': T_degC}
        try:
            x0 = base[param]
        except KeyError:
            raise ValueError(f"param must be one of {sorted(base)}, got {param!r}")

        def s2(x):
            v = {'wl1': (x, wl2), 'wl2': (wl1, x), 'wl': (x, x)}.get(param, (wl1, wl2))
            ang = x if param == 'theta' else angle_rad
            T = x if param == 'T' else T_degC
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', ValidityWarning)
                dk = self.dk_sfg(v[0], v[1], ang, T, pol1, pol2, pol3)
            u = 0.5 * dk * L_um
            return float(np.sinc(u / pi) ** 2)

        s0 = s2(x0)
        if s0 < 0.5:
            raise ValueError(
                f"sinc2 = {s0:.3g} < 0.5 at the given point - it is not "
                f"phase-matched, so no FWHM contains it. Evaluate at an angle "
                f"from pmAngles_sfg (in radians).")

        # initial step ~ the variable's natural fine scale
        step0 = {'wl1': 1e-5, 'wl2': 1e-5, 'wl': 1e-5, 'theta': 1e-5, 'T': 1e-2}[param]

        def half_point(sign):
            # expand until sinc² < 0.5, then bisect the crossing
            a, step = x0, step0
            for _ in range(64):
                b = a + sign * step
                if s2(b) < 0.5:
                    return brentq(lambda x: s2(x) - 0.5, *sorted((a, b)))
                a, step = b, step * 2
            return np.inf
        lo, hi = half_point(-1), half_point(+1)
        return np.inf if np.isinf(lo) or np.isinf(hi) else float(hi - lo)

    def wl_idler(self, wl_p, wl_s):
        """
        Idler wavelength of DFG/OPA: 1/wl_i = 1/wl_p - 1/wl_s (µm).

        Raises ValueError unless wl_s > wl_p everywhere - the signal must carry
        less energy per photon than the pump for an idler to exist.
        """
        wl_p = np.asarray(wl_p, dtype=float)
        wl_s = np.asarray(wl_s, dtype=float)
        if np.any(wl_s <= wl_p):
            raise ValueError("DFG/OPA needs wl_s > wl_p (the idler 1/(1/wl_p - 1/wl_s) "
                             "must be positive and finite)")
        out = 1. / (1. / wl_p - 1. / wl_s)
        return out.item() if out.ndim == 0 else out

    def dk_dfg(self, wl_p, wl_s, angle_rad, T_degC, pol_s, pol_i, pol_p):
        """
        Wavevector mismatch k_p - k_s - k_i (rad/µm) for DFG / OPA / OPO.

        Identical to ``dk_sfg(wl_s, wl_i, ...)`` with wl_i = ``wl_idler``; see
        that method for the arguments. Polarizations are given in the order
        (signal, idler, pump).
        """
        return self.dk_sfg(wl_s, self.wl_idler(wl_p, wl_s), angle_rad, T_degC, pol_s, pol_i, pol_p)

    def pmFactor_dfg(self, wl_p, wl_s, angle_rad, T_degC, pol_s, pol_i, pol_p, L_mm):
        """sinc²(Δk L / 2) for DFG / OPA / OPO; see ``pmFactor_sfg``."""
        return self.pmFactor_sfg(wl_s, self.wl_idler(wl_p, wl_s), angle_rad, T_degC,
                                 pol_s, pol_i, pol_p, L_mm)

    def qpm_period_dfg(self, wl_p, wl_s, angle_rad, T_degC, pol_s, pol_i, pol_p, order=1):
        """Quasi-phase-matching period (µm) for DFG / OPA / OPO; see ``qpm_period_sfg``."""
        return self.qpm_period_sfg(wl_s, self.wl_idler(wl_p, wl_s), angle_rad, T_degC,
                                   pol_s, pol_i, pol_p, order)

    def pmAngles_dfg(self, wl_p, wl_s, T_degC, tol_deg=0.001, deg=False):
        """
        Phase-matching angles for DFG / OPA / OPO at a given pump and signal.

        The same dictionary as ``pmAngles_sfg(wl_s, wl_i, ...)`` with the key
        'wl3' replaced by 'wl_i', the idler wavelength in µm. Polarization keys read
        (signal, idler, pump): 'ooe' is Type I, 'oee' and 'eoe' Type II.
        """
        wl_i = self.wl_idler(wl_p, wl_s)
        d = self.pmAngles_sfg(wl_s, wl_i, T_degC, tol_deg=tol_deg, deg=deg)
        del d['wl3']
        return {'wl_i': wl_i, **d}

    def tuning_dfg(self, wl_p, angle_rad, T_degC, pol_s, pol_i, pol_p,
                   wl_i_max=None, tol_um=1e-6, n_grid=2001, qpm_period=None, qpm_order=1):
        """
        Signal/idler pairs that phase-match at a fixed angle and temperature -
        one point of an OPO/OPA tuning curve.

        Solves Δk(wl_s) = 0 - or, with ``qpm_period`` given, the quasi-phase-
        matching condition |Δk(wl_s)| = 2π·qpm_order/Λ - with 1/wl_i = 1/wl_p
        - 1/wl_s, for the signal on
        the short-wavelength side of degeneracy (wl_s <= wl_i). The search runs
        on an even grid in signal frequency from degeneracy (wl_s = 2 wl_p) to
        the signal whose idler is ``wl_i_max``; each sign change is refined with
        Brent's method, so every root is returned - a Type II branch or a
        retracing angle-tuning curve can give more than one. Degeneracy itself,
        where Δk touches zero without crossing, is reported when |Δk| there is
        smaller than its change over the first grid step.

        Parameters
        ----------
        wl_p : float
            Pump wavelength in µm.
        angle_rad : float
            theta or phi angle in radians (see ``n``).
        T_degC : float
            Crystal temperature in degC.
        pol_s, pol_i, pol_p : {'o', 'e'}
            Polarizations of signal, idler and pump. The branch with signal
            and idler polarizations exchanged is the triple with pol_s and
            pol_i swapped.
        wl_i_max : float, optional
            Longest idler wavelength searched, in µm. Default: the medium's
            Sellmeier validity edge (at most 20), so spurious far-extrapolation
            roots are excluded. Pass a value explicitly to extrapolate beyond
            the validity range on purpose.
        tol_um : float, default 1e-6
            Absolute tolerance of the returned signal wavelength, in µm.
        n_grid : int, default 2001
            Points of the coarse frequency grid.
        qpm_period : float, optional
            Poling period Λ in µm. When given, the condition solved is
            |Δk| = 2π m/Λ (m = ``qpm_order``) instead of Δk = 0 - the tuning of
            a periodically poled crystal at fixed Λ, e.g. PPLN's signal versus
            temperature with ``('e', 'e', 'e')`` at angle_rad = pi/2.
        qpm_order : int, default 1
            QPM order m.

        Return
        ------
        list of (wl_s, wl_i) tuples, in µm, sorted by wl_s; empty if none.
        """
        if wl_i_max is None:
            # default search stops at the Sellmeier validity edge: roots from
            # far extrapolation are spurious. Pass wl_i_max explicitly to
            # search further on purpose.
            wl_i_max = min(20.0, self._wl_range[1]) if self._wl_range else 20.0
        if wl_i_max <= 2 * wl_p:
            return []
        nu_p = 1. / wl_p
        nu_s = np.linspace(nu_p / 2, nu_p - 1. / wl_i_max, n_grid)   # even in frequency
        wl_s_ar = 1. / nu_s
        # the residual that must vanish: dk itself, or |dk| less the grating vector
        k_g = 0.0 if qpm_period is None else 2 * pi * qpm_order / qpm_period

        def resid(w):
            # this method extrapolates up to wl_i_max by design (see
            # docstring); per-point ValidityWarnings would only be noise
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', ValidityWarning)
                d = self.dk_dfg(wl_p, w, angle_rad, T_degC, pol_s, pol_i, pol_p)
            return d if qpm_period is None else np.abs(d) - k_g
        try:
            with np.errstate(invalid='ignore', divide='ignore'):
                dk_ar = np.asarray(resid(wl_s_ar), dtype=float)
        except ValueError:
            return []                         # no Sellmeier equation for a polarization
        ok = np.isfinite(dk_ar)
        dk = lambda w: float(resid(w))
        roots = []
        sign = np.signbit(dk_ar)
        for i in np.nonzero(np.diff(sign) & ok[:-1] & ok[1:])[0]:
            # brentq works in wavelength; the bracket is one grid cell
            a, b = wl_s_ar[i + 1], wl_s_ar[i]       # wl decreases along the grid
            roots.append(brentq(dk, a, b, xtol=tol_um))
        # degeneracy: Δk(nu_s) is stationary there for Type I (signal and
        # idler are the same wave), so a root at 2 wl_p never shows as a sign
        # change on the half-grid. Accept it when |Δk| is within one grid
        # step's change of zero.
        if ok[0] and ok[1] and not (roots and abs(roots[-1] - wl_s_ar[0]) <= tol_um):
            if abs(dk_ar[0]) <= abs(dk_ar[1] - dk_ar[0]):
                roots.append(float(wl_s_ar[0]))
        roots.sort()
        return [(w, float(self.wl_idler(wl_p, w))) for w in roots]
    
 