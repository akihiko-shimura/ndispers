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
