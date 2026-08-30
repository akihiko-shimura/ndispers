"""
Second-order nonlinearity shared by every point group: Miller-rule scaling of
the tensor components and the contraction that gives d_eff.

Theory and conventions: docs/dev/deff_theory.tex. In short, the tensor d_il
is a property of the crystal in its principal dielectric frame; Miller's rule
scales each component with the principal susceptibilities chi_ii = n_i^2 - 1
of the three waves, and d_eff is the contraction of that tensor with the unit
E-field vectors of the three waves, which carry all the angular dependence
(including walk-off). Point-group classes supply only the tensor pattern;
crystals supply only the reference measurements.
"""
import re

import numpy as np

from ndispers._baseclass import Medium

# contracted (Voigt) index of a pair of Cartesian indices
_L = {"xx": 0, "yy": 1, "zz": 2, "yz": 3, "zy": 3, "zx": 4, "xz": 4, "xy": 5, "yx": 5}
_AX = {"x": 0, "y": 1, "z": 2}


class NonlinearGroup(Medium):
    """
    Base class for crystals whose point group allows a second-order
    nonlinearity. Subclasses set ``_d_entries`` (the point group) and
    ``_n_axis`` (uniaxial or biaxial); crystals set ``_d_ref``.

    Class attributes
    ----------------
    _d_entries : dict
        Tensor entries of each independent component, Kleinman symmetry
        already applied: ``{"d22": [("yyy", +1), ("yxx", -1), ("xxy", -1)], ...}``.
        Each entry names d_ijk in the principal frame with its sign relative
        to the component; the first entry defines the component for Miller's
        rule. One entry per d_il cell (do not list both "xzx" and "xxz").
    _d_ref : dict
        Reference measurements set by the crystal:
        ``{"d22": (value_pm_per_V, wl1_um, wl2_um)}`` - the value of d_il for
        SFG of the two wavelengths (SHG when equal). Components of
        ``_d_entries`` absent from ``_d_ref`` are treated as zero in d_eff;
        ``deff_sfg`` raises if such a component would actually contribute.
    _d_note : str
        Where the reference values come from and how far to trust Miller
        scaling for this crystal.
    """
    __slots__ = []
    _d_entries = {}
    _d_ref = {}
    _d_note = ""

    def _n_axis(self, axis, wl_um, T_degC):
        """Principal refractive index n_x, n_y or n_z."""
        raise NotImplementedError

    def _chi(self, axis, wl_um, T_degC):
        """Principal linear susceptibility chi_ii = n_i**2 - 1."""
        return self._n_axis(axis, wl_um, T_degC)**2 - 1

    def _chi3(self, il, wl1, wl2, T_degC):
        """chi_ii(wl3) chi_jj(wl1) chi_kk(wl2) for the component's defining ijk."""
        i, j, k = self._d_entries[il][0][0]
        wl3 = 1. / (1. / wl1 + 1. / wl2)
        return (self._chi(i, wl3, T_degC) * self._chi(j, wl1, T_degC)
                * self._chi(k, wl2, T_degC))

    def miller_delta(self, il, T_degC):
        """
        Miller's delta of a tensor component, in pm/V (epsilon_0 absorbed).

        Parameters
        ----------
        il : str
            Component name, e.g. ``"d22"``; one of ``_d_ref``'s keys.
        T_degC : float
            Temperature at which the susceptibilities are evaluated.
        """
        d0, wl1, wl2 = self._d_ref[il]
        return d0 / self._chi3(il, wl1, wl2, T_degC)

    def d_sfg(self, il, wl1, wl2, T_degC):
        """
        Tensor component d_il (pm/V) for SFG of wl1 and wl2, scaled from the
        reference measurement by Miller's rule. Exact at the reference pair.

        Parameters
        ----------
        il : str
            Component name, e.g. ``"d22"``; one of ``_d_ref``'s keys.
        wl1, wl2 : float
            Pump wavelengths in µm.
        T_degC : float
            Crystal temperature in degC.
        """
        return self.miller_delta(il, T_degC) * self._chi3(il, wl1, wl2, T_degC)

    def __getattr__(self, name):
        # d22_sfg(wl1, wl2, T) etc., one per reference component
        m = re.fullmatch(r"(d\d\d)_sfg", name)
        if m and m.group(1) in self._d_ref:
            il = m.group(1)
            return lambda wl1, wl2, T_degC: self.d_sfg(il, wl1, wl2, T_degC)
        raise AttributeError(f"{type(self).__name__!r} object has no attribute {name!r}")

    def _fixed_angles(self, theta_rad, phi_rad):
        """Fill in, or check, the angle a principal plane holds fixed."""
        out = []
        for name, given in (("theta_rad", theta_rad), ("phi_rad", phi_rad)):
            fixed = getattr(self, name)
            if isinstance(fixed, str):          # 'var' or 'arb': the caller's value stands
                if given is None:
                    raise ValueError(f"{name} is required for {type(self).__name__}")
                out.append(given)
            elif given is None:
                out.append(fixed)
            elif np.allclose(given, fixed):
                out.append(fixed)
            else:
                raise ValueError(f"{name} is fixed at {fixed:.6g} rad in the "
                                 f"{self.plane} plane of {type(self).__name__}")
        return tuple(out)

    def _evec(self, pol, wl_um, theta_rad, phi_rad, T_degC):
        """
        Unit E-field vector of an 'o' or 'e' wave in the principal frame.

        With alpha the angle that varies for this medium (theta, or phi in the
        xy plane) and k(alpha) the wavevector direction, the e-wave field is
        dk/dalpha evaluated at alpha + rho (rho the walk-off angle, so that E
        is perpendicular to the Poynting vector) and the o-wave field is
        k x dk/dalpha. This reproduces the textbook vectors for uniaxial
        crystals and for the three principal planes of a biaxial crystal.
        """
        th, ph = np.asarray(theta_rad, float), np.asarray(phi_rad, float)
        if self.phi_rad == 'var':               # principal plane xy, theta = pi/2
            k = lambda p: np.stack(np.broadcast_arrays(np.cos(p), np.sin(p), 0 * p))
            dk = lambda p: np.stack(np.broadcast_arrays(-np.sin(p), np.cos(p), 0 * p))
            a = ph
            rho = self.woa_phi(wl_um, ph, T_degC, pol='e') if pol == 'e' else 0.
        else:                                   # theta varies
            k = lambda t: np.stack(np.broadcast_arrays(
                np.sin(t) * np.cos(ph), np.sin(t) * np.sin(ph), np.cos(t)))
            dk = lambda t: np.stack(np.broadcast_arrays(
                np.cos(t) * np.cos(ph), np.cos(t) * np.sin(ph), -np.sin(t)))
            a = th
            rho = self.woa_theta(wl_um, th, T_degC, pol='e') if pol == 'e' else 0.
        if pol == 'e':
            return dk(a + rho)
        if pol == 'o':
            return np.cross(k(a), dk(a), axis=0)
        raise ValueError(f"pol = {pol!r} must be 'o' or 'e'")

    def eta_sfg(self, beam1, beam2, theta_rad, phi_rad, T_degC, pol3, L_mm,
                deff_pmV=None, dk_offset=0.0, n_grid=50, details=False):
        """
        Semi-analytical SFG conversion efficiency with pump depletion.

        Armstrong's plane-wave elliptic solution averaged over the Gaussian
        profiles of two collimated beams; neglects walk-off, group-velocity
        mismatch, diffraction and GVD (a ``ModelValidityWarning`` fires when
        those are not small). Theory: docs/dev/efficiency_theory.pdf.
        Requires scipy (``pip install ndispers[eff]``).

        Parameters
        ----------
        beam1, beam2 : ndispers.PulsedBeam or ndispers.CWBeam
            The two input beams; wavelength, energy/power, radius, duration
            and polarization live on the beam. Both must be the same type -
            it decides pulsed (energy efficiency) vs cw (power efficiency).
        theta_rad, phi_rad : float or None
            Propagation angles, as in ``deff_sfg``: n and dk use the
            varying angle of this medium, d_eff uses both. Pass None for
            the angle a principal plane fixes; with ``deff_pmV`` given,
            the azimuth may be None.
        T_degC : float
        pol3 : {'o', 'e'}
            Polarization of the sum-frequency wave.
        L_mm : float
            Crystal length.
        deff_pmV : float, optional
            Override the Miller-scaled ``deff_sfg`` with a measured value.
        dk_offset : float, optional
            Residual phase mismatch added to ``dk_sfg`` (rad/um).
        n_grid : int, optional
            Simpson mesh per integration axis. The default resolves equal
            beams to ~1e-4 relative; strongly unequal beam sizes or
            durations may need more (the quadrature is normalized to
            beam 1, so swapping very unequal beams changes the result by
            the quadrature error, ~3e-3 at the default).
        details : bool, optional
            False (default): return the efficiency eta = out/(in1+in2) as a
            float. True: return a dict with eta, per-wave photon conversion,
            output energy/power, peak intensities, deff, dk, K1, and the
            neglected-effect ratios under 'model_ratios'.

        Examples
        --------
        Third-harmonic stage of a Q-switched Nd:YAG laser (1064 + 532 ->
        355 nm, type I, photon-balanced inputs, d_eff Miller-scaled
        automatically):

        >>> import numpy as np
        >>> import ndispers as nd
        >>> bbo = nd.media.crystals.BetaBBO_Tamosauskas2018()
        >>> ooe = bbo.pmAngles_sfg(1.064, 0.532, 25)['ooe']
        >>> theta_pm, = ooe['theta']
        >>> b1 = nd.PulsedBeam(wl_um=1.064, E_uJ=5e3, w_um=1500.0,
        ...                    t_fs=8e6, pol='o')
        >>> b2 = nd.PulsedBeam(wl_um=0.532, E_uJ=10e3, w_um=1500.0,
        ...                    t_fs=8e6, pol='o')
        >>> eta = bbo.eta_sfg(b1, b2, theta_pm, np.pi/2, 25, 'e', 5.0)
        >>> print(round(eta, 4))
        0.1657
        """
        from ndispers import _efficiency
        return _efficiency.eta_sfg(self, beam1, beam2, theta_rad, phi_rad,
                                   T_degC, pol3, L_mm, deff_pmV=deff_pmV,
                                   dk_offset=dk_offset, n_grid=n_grid,
                                   details=details)

    def eta_shg(self, beam, theta_rad, phi_rad, T_degC, pol3, L_mm,
                deff_pmV=None, dk_offset=0.0, n_grid=50, details=False):
        """
        Semi-analytical type-I SHG conversion efficiency with depletion.

        The single-field Armstrong solution (tanh-squared at dk = 0)
        averaged over one Gaussian beam/pulse; the degeneracy factor is
        handled here, so this is NOT ``eta_sfg`` called with the same beam
        twice. Type-II SHG is degenerate SFG of the o- and e-projections of
        one beam - use ``eta_sfg`` for it. Arguments as in ``eta_sfg`` with
        a single beam whose ``pol`` is the fundamental polarization; the
        efficiency is out/in of that beam. Theory:
        docs/dev/efficiency_theory.pdf. Requires scipy
        (``pip install ndispers[eff]``).

        Examples
        --------
        Doubling a Q-switched Nd:YAG laser (30 mJ, 8 ns, w = 1.5 mm,
        7 mm crystal - deep depletion):

        >>> import numpy as np
        >>> import ndispers as nd
        >>> bbo = nd.media.crystals.BetaBBO_Tamosauskas2018()
        >>> ooe = bbo.pmAngles_sfg(1.064, 1.064, 25)['ooe']
        >>> theta_pm, = ooe['theta']
        >>> b = nd.PulsedBeam(wl_um=1.064, E_uJ=30e3, w_um=1500.0,
        ...                   t_fs=8e6, pol='o')
        >>> print(round(bbo.eta_shg(b, theta_pm, np.pi/2, 25, 'e', 7.0), 4))
        0.2832

        Doubling 100 fs pulses at 800 nm in 1 mm of BBO raises
        ``ModelValidityWarning``: the group-velocity mismatch (~190 fs/mm)
        exceeds the pulse duration, which is why femtosecond SHG uses
        0.1-0.3 mm crystals.
        """
        from ndispers import _efficiency
        return _efficiency.eta_shg(self, beam, theta_rad, phi_rad, T_degC,
                                   pol3, L_mm, deff_pmV=deff_pmV,
                                   dk_offset=dk_offset, n_grid=n_grid,
                                   details=details)

    def deff_sfg(self, wl1, wl2, theta_rad, phi_rad, T_degC, pol1, pol2, pol3):
        """
        Effective second-order nonlinear coefficient d_eff (pm/V) for
        collinear SFG, with the tensor components scaled to the pump
        wavelengths by Miller's rule and walk-off included in the field
        directions. Overall sign is a convention; compare magnitudes.

        Parameters
        ----------
        wl1, wl2 : float
            Pump wavelengths in µm.
        theta_rad, phi_rad : float or array_like or None
            Polar and azimuthal angles of the wavevector. A principal-plane
            medium fixes one of them: pass ``None`` to use the plane's value
            (passing that value explicitly is also accepted).
        T_degC : float
            Crystal temperature in degC.
        pol1, pol2, pol3 : {'o', 'e'}
            Polarizations of the two pumps and of the sum-frequency wave.
        """
        theta_rad, phi_rad = self._fixed_angles(theta_rad, phi_rad)
        wl3 = 1. / (1. / wl1 + 1. / wl2)
        e1 = self._evec(pol1, wl1, theta_rad, phi_rad, T_degC)
        e2 = self._evec(pol2, wl2, theta_rad, phi_rad, T_degC)
        e3 = self._evec(pol3, wl3, theta_rad, phi_rad, T_degC)
        u = [e1[0] * e2[0], e1[1] * e2[1], e1[2] * e2[2],
             e1[1] * e2[2] + e1[2] * e2[1],
             e1[2] * e2[0] + e1[0] * e2[2],
             e1[0] * e2[1] + e1[1] * e2[0]]

        def contract(d):
            D = [[0.] * 6 for _ in range(3)]
            for il, entries in self._d_entries.items():
                for ijk, sign in entries:
                    D[_AX[ijk[0]]][_L[ijk[1:]]] += sign * d.get(il, 0.)
            return sum(e3[i] * sum(D[i][l] * u[l] for l in range(6)) for i in range(3))

        # a component this crystal has no reference value for must not be
        # silently zero if the geometry would actually use it
        missing = {il: 1. for il in self._d_entries if il not in self._d_ref}
        if missing and np.any(np.abs(contract(missing)) > 1e-12):
            raise ValueError(
                f"{type(self).__name__} has no reference value for "
                f"{', '.join(missing)}, which {pol1}{pol2}{pol3} needs. {self._d_note}")

        out = contract({il: self.d_sfg(il, wl1, wl2, T_degC) for il in self._d_ref})
        return out.item() if np.ndim(out) == 0 else out


    # DFG / OPA / OPO: the same tensor contraction. Under Kleinman symmetry d is
    # invariant to permuting the three waves, and Miller's chi(w3) chi(w1) chi(w2)
    # is symmetric too, so d_eff for the parametric process equals the SFG one
    # evaluated at (signal, idler). Indices (1, 2, 3) = (signal, idler, pump).

    def d_dfg(self, il, wl_p, wl_s, T_degC):
        """Tensor component d_il (pm/V) for DFG/OPA of pump wl_p and signal wl_s;
        see ``d_sfg``. Equals ``d_sfg(il, wl_s, wl_i, T)``."""
        return self.d_sfg(il, wl_s, self.wl_idler(wl_p, wl_s), T_degC)

    def deff_dfg(self, wl_p, wl_s, theta_rad, phi_rad, T_degC, pol_s, pol_i, pol_p):
        """Effective nonlinear coefficient d_eff (pm/V) for collinear DFG / OPA /
        OPO; see ``deff_sfg``. Polarizations read (signal, idler, pump)."""
        return self.deff_sfg(wl_s, self.wl_idler(wl_p, wl_s), theta_rad, phi_rad, T_degC,
                             pol_s, pol_i, pol_p)


class UniaxialGroup(NonlinearGroup):
    """Uniaxial: chi_xx = chi_yy from n_o, chi_zz from the principal n_e."""
    __slots__ = []

    def _n_axis(self, axis, wl_um, T_degC):
        if axis == 'z':
            return self.n(wl_um, 0.5 * np.pi, T_degC, pol='e')
        return self.n(wl_um, 0., T_degC, pol='o')


class BiaxialGroup(NonlinearGroup):
    """Biaxial: principal indices from the crystal's n_x_expr, n_y_expr, n_z_expr."""
    __slots__ = []

    def _n_axis(self, axis, wl_um, T_degC):
        def expr(pol):
            return getattr(self, f"n_{axis}_expr")()
        expr.__name__ = f"n_{axis}"             # the lambdify cache key
        return self._func(expr, wl_um, 0., 0., T_degC)
