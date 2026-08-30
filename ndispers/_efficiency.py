"""
Semi-analytical conversion efficiency of pulsed and cw SFG / type-I SHG.

Armstrong's Jacobi-elliptic plane-wave solution with pump depletion and
phase mismatch, averaged over Gaussian spatial (and, for pulses, temporal)
profiles. Neglects spatial walk-off, group-velocity mismatch, diffraction
and GVD; ModelValidityWarning fires when those are not small. Theory,
conventions and the validation of every formula against Armstrong (1962)
and Boyd (2020): docs/dev/efficiency_theory.tex.

Requires scipy (Jacobi sn): pip install ndispers[eff].
"""
from dataclasses import dataclass

import numpy as np

c_ms = 2.99792458e8
eps0 = 8.8541878128e-12
_SQRT2 = np.sqrt(2.0)
_FWHM_TO_1E = 2.0 * np.sqrt(np.log(2.0))
_R_MAX = 5.0     # integration bound in 1/e radii / half-durations
_RATIO_WARN = 0.3


class ModelValidityWarning(UserWarning):
    """A physical effect this semi-analytical model neglects (walk-off,
    group-velocity mismatch, diffraction or GVD) is not small for the given
    parameters; the returned efficiency is an overestimate. The ratios are
    in the details=True dict under 'model_ratios'."""


@dataclass(frozen=True)
class PulsedBeam:
    """A collimated Gaussian pulse train's single pulse.

    wl_um: vacuum wavelength (µm). E_uJ: pulse energy (µJ).
    w_um: 1/e² intensity radius (µm), the Gaussian w0 with I ∝ exp(−2r²/w²).
    t_fs: FWHM intensity duration (fs). pol: 'o' or 'e'.
    """
    wl_um: float
    E_uJ: float
    w_um: float
    t_fs: float
    pol: str

    def __post_init__(self):
        _check_beam(self, ('E_uJ', 't_fs'))

    @property
    def I_peak_Wm2(self):
        """Peak intensity (W/m²): E/(π r0² √π t0), 1/e radius and half-duration."""
        r0 = self.w_um * 1e-6 / _SQRT2
        t0 = self.t_fs * 1e-15 / _FWHM_TO_1E
        return self.E_uJ * 1e-6 / (np.pi * r0**2 * np.sqrt(np.pi) * t0)


@dataclass(frozen=True)
class CWBeam:
    """A collimated Gaussian cw beam.

    wl_um: vacuum wavelength (µm). P_W: power (W).
    w_um: 1/e² intensity radius (µm). pol: 'o' or 'e'.
    """
    wl_um: float
    P_W: float
    w_um: float
    pol: str

    def __post_init__(self):
        _check_beam(self, ('P_W',))

    @property
    def I_peak_Wm2(self):
        """Axial intensity (W/m²): P/(π r0²), 1/e radius."""
        r0 = self.w_um * 1e-6 / _SQRT2
        return self.P_W / (np.pi * r0**2)


def _check_beam(b, extra_fields):
    if b.pol not in ('o', 'e'):
        raise ValueError(f"pol must be 'o' or 'e', got {b.pol!r}")
    for f in ('wl_um', 'w_um') + extra_fields:
        v = getattr(b, f)
        if not (isinstance(v, (int, float)) and v > 0):
            raise ValueError(f"{type(b).__name__}.{f} must be a positive "
                             f"number, got {v!r}")


def _ellipj_sn(u, m):
    try:
        from scipy.special import ellipj
    except ImportError:
        raise ImportError(
            "the efficiency methods need scipy for the Jacobi sn function: "
            "pip install ndispers[eff]") from None
    # m -> 1 at strong matched conversion; round-off can push it above
    return ellipj(u, np.clip(m, 0.0, 1.0))[0]


def _y_pw(alpha, beta, sigma):
    """Normalized SF intensity y(ζ=1) of the plane-wave Armstrong solution:
    (dy/dζ)² = 4y(α−y)(β−y) − σ²y², y(0)=0  →  y = b1 sn²(√b2; b1/b2).
    Verified against Armstrong (1962) eqs. (6.6)-(6.13); see the theory note.
    """
    a1 = alpha + beta + 0.25 * sigma**2
    a3 = np.sqrt(np.maximum(a1**2 - 4.0 * alpha * beta, 0.0))
    b1 = 0.5 * (a1 - a3)
    b2 = 0.5 * (a1 + a3)
    with np.errstate(divide='ignore', invalid='ignore'):
        m = np.where(b2 > 0, b1 / np.where(b2 > 0, b2, 1.0), 0.0)
    return b1 * _ellipj_sn(np.sqrt(b2), m) ** 2


def _simpson(f, x):
    from scipy.integrate import simpson
    return simpson(f, x=x)


def _characteristic_intensity(n1, n2, n3, wl1_um, wl2_um, deff_pmV, L_mm):
    """I_a = c ε0 n1 n2 n3 λ1 λ2 / (8π² deff² L²), all SI, in W/m²."""
    return (c_ms * eps0 * n1 * n2 * n3 * wl1_um * 1e-6 * wl2_um * 1e-6
            / (8.0 * np.pi**2 * (deff_pmV * 1e-12)**2 * (L_mm * 1e-3)**2))


def _angles_for(medium, theta_rad, phi_rad):
    """The varying angle (for n, dk) from the (theta, phi) pair."""
    if medium.theta_rad == 'var':
        return theta_rad
    if medium.phi_rad == 'var':
        return phi_rad
    return theta_rad


def _model_ratios(medium, angle, T_degC, beams, pol3, wl3, L_mm):
    """How large the neglected effects are; None where not computable."""
    L_um = L_mm * 1e3
    out = {}
    try:      # spatial walk-off vs beam radius (e-waves walk off)
        woas = []
        for wl, pol in [(b.wl_um, b.pol) for b in beams] + [(wl3, pol3)]:
            for f in (medium.woa_theta, medium.woa_phi):
                try:
                    woas.append(abs(float(f(wl, angle, T_degC, pol=pol))))
                except Exception:
                    pass
        out['walkoff'] = max(woas) * L_um / min(b.w_um for b in beams) \
            if woas else None
    except Exception:
        out['walkoff'] = None
    try:      # diffraction: L / Rayleigh range
        zr = min(np.pi * float(medium.n(b.wl_um, angle, T_degC, pol=b.pol))
                 * b.w_um**2 / b.wl_um for b in beams)
        out['diffraction'] = L_um / zr
    except Exception:
        out['diffraction'] = None
    pulsed = hasattr(beams[0], 't_fs')
    try:      # temporal walk-off between any pair of waves
        if pulsed:
            t_min = min(b.t_fs for b in beams)
            ws = [(b.wl_um, b.pol) for b in beams] + [(wl3, pol3)]
            gvms = [abs(float(medium.GVM(wa, wb, angle, T_degC,
                                         pol1=pa, pol2=pb)))
                    for i, (wa, pa) in enumerate(ws)
                    for (wb, pb) in ws[i + 1:]]
            out['gvm'] = max(gvms) * L_mm / t_min
        else:
            out['gvm'] = None
    except Exception:
        out['gvm'] = None
    try:      # pulse spreading by GVD
        if pulsed:
            out['gvd'] = max(abs(float(medium.GVD(b.wl_um, angle, T_degC,
                                                  pol=b.pol))) * L_mm
                             / min(b.t_fs for b in beams)**2 for b in beams)
        else:
            out['gvd'] = None
    except Exception:
        out['gvd'] = None
    return out


def _warn_ratios(ratios):
    import warnings
    bad = {k: v for k, v in ratios.items()
           if v is not None and v > _RATIO_WARN}
    if bad:
        worst = max(bad, key=bad.get)
        warnings.warn(
            f"neglected effect(s) not small: "
            + ", ".join(f"{k} ratio {v:.2g}" for k, v in bad.items())
            + f"; the semi-analytical efficiency is an overestimate "
            f"(worst: {worst}). See docs/dev/efficiency_theory.tex.",
            ModelValidityWarning, stacklevel=4)


def eta_sfg(medium, beam1, beam2, theta_rad, phi_rad, T_degC, pol3, L_mm,
            deff_pmV=None, dk_offset=0.0, n_grid=50, details=False):
    """See NonlinearGroup.eta_sfg."""
    if type(beam1) is not type(beam2):
        raise TypeError(
            "beam1 and beam2 must both be PulsedBeam or both CWBeam "
            f"(got {type(beam1).__name__} and {type(beam2).__name__}); this "
            "model has no mixed pulsed/cw solution")
    pulsed = isinstance(beam1, PulsedBeam)
    if not pulsed and not isinstance(beam1, CWBeam):
        raise TypeError(f"beams must be PulsedBeam or CWBeam, got "
                        f"{type(beam1).__name__}")
    wl1, wl2 = beam1.wl_um, beam2.wl_um
    wl3 = 1.0 / (1.0 / wl1 + 1.0 / wl2)
    angle = _angles_for(medium, theta_rad, phi_rad)

    n1 = float(medium.n(wl1, angle, T_degC, pol=beam1.pol))
    n2 = float(medium.n(wl2, angle, T_degC, pol=beam2.pol))
    n3 = float(medium.n(wl3, angle, T_degC, pol=pol3))
    if deff_pmV is None:
        deff_pmV = float(medium.deff_sfg(wl1, wl2, theta_rad, phi_rad,
                                         T_degC, beam1.pol, beam2.pol, pol3))
    dk = float(medium.dk_sfg(wl1, wl2, angle, T_degC,
                             beam1.pol, beam2.pol, pol3)) + dk_offset
    sigma = dk * L_mm * 1e3          # rad/µm × µm

    I_a = _characteristic_intensity(n1, n2, n3, wl1, wl2, abs(deff_pmV), L_mm)
    I10, I20 = beam1.I_peak_Wm2, beam2.I_peak_Wm2
    K1, p1 = I10 / I_a, I20 / I10
    d = wl1 / wl3
    d1 = wl2 / wl3
    rho = beam1.w_um / beam2.w_um

    x = np.linspace(0.0, _R_MAX, n_grid)
    if pulsed:
        tau = beam1.t_fs / beam2.t_fs
        r1, t1 = np.meshgrid(x, x, indexing='ij')
        G1 = np.exp(-r1**2) * np.exp(-t1**2)
        G2 = np.exp(-rho**2 * r1**2) * np.exp(-tau**2 * t1**2)
        y = _y_pw(d * K1 * G1, d1 * p1 * K1 * G2, sigma)
        integral = _simpson(np.array([_simpson(y[i] * x[i], x)
                                      for i in range(n_grid)]), x)
        eta = 4.0 / (np.sqrt(np.pi) * K1
                     * (1.0 + p1 / (rho**2 * tau))) * integral
    else:
        G1 = np.exp(-x**2)
        G2 = np.exp(-rho**2 * x**2)
        y = _y_pw(d * K1 * G1, d1 * p1 * K1 * G2, sigma)
        integral = _simpson(y * x, x)
        eta = 2.0 / (K1 * (1.0 + p1 / rho**2)) * integral
    eta = float(eta)

    ratios = _model_ratios(medium, angle, T_degC, (beam1, beam2),
                           pol3, wl3, L_mm)
    _warn_ratios(ratios)
    if not details:
        return eta
    in1 = beam1.E_uJ if pulsed else beam1.P_W
    in2 = beam2.E_uJ if pulsed else beam2.P_W
    out3 = eta * (in1 + in2)
    return {
        'eta': eta,
        'eta_photon_1': out3 * wl3 / (in1 * wl1),
        'eta_photon_2': out3 * wl3 / (in2 * wl2),
        ('E3_uJ' if pulsed else 'P3_W'): out3,
        'I1_peak_MWcm2': I10 * 1e-10, 'I2_peak_MWcm2': I20 * 1e-10,
        'deff_pmV': deff_pmV, 'dk_radum': dk, 'K1': K1,
        'model_ratios': ratios,
    }


def eta_shg(medium, beam, theta_rad, phi_rad, T_degC, pol3, L_mm,
            deff_pmV=None, dk_offset=0.0, n_grid=50, details=False):
    """See NonlinearGroup.eta_shg."""
    if not isinstance(beam, (PulsedBeam, CWBeam)):
        raise TypeError(f"beam must be PulsedBeam or CWBeam, got "
                        f"{type(beam).__name__}")
    pulsed = isinstance(beam, PulsedBeam)
    wl1 = beam.wl_um
    wl3 = wl1 / 2.0
    angle = _angles_for(medium, theta_rad, phi_rad)

    n1 = float(medium.n(wl1, angle, T_degC, pol=beam.pol))
    n3 = float(medium.n(wl3, angle, T_degC, pol=pol3))
    if deff_pmV is None:
        deff_pmV = float(medium.deff_sfg(wl1, wl1, theta_rad, phi_rad,
                                         T_degC, beam.pol, beam.pol, pol3))
    dk = float(medium.dk_sfg(wl1, wl1, angle, T_degC,
                             beam.pol, beam.pol, pol3)) + dk_offset
    sigma = dk * L_mm * 1e3

    # I_a^SHG: the SFG normalization under the degenerate substitution,
    # no extra factor (pinned against Boyd (2.7.36)/(2.7.38); theory note §6)
    I_a = _characteristic_intensity(n1, n1, n3, wl1, wl1, abs(deff_pmV), L_mm)
    K = beam.I_peak_Wm2 / I_a

    x = np.linspace(0.0, _R_MAX, n_grid)
    if pulsed:
        r1, t1 = np.meshgrid(x, x, indexing='ij')
        G = np.exp(-r1**2) * np.exp(-t1**2)
        y = _y_pw(K * G, K * G, sigma)
        integral = _simpson(np.array([_simpson(y[i] * x[i], x)
                                      for i in range(n_grid)]), x)
        eta = 4.0 / (np.sqrt(np.pi) * K) * integral
    else:
        G = np.exp(-x**2)
        y = _y_pw(K * G, K * G, sigma)
        eta = 2.0 / K * _simpson(y * x, x)
    eta = float(eta)

    ratios = _model_ratios(medium, angle, T_degC, (beam,), pol3, wl3, L_mm)
    _warn_ratios(ratios)
    if not details:
        return eta
    in1 = beam.E_uJ if pulsed else beam.P_W
    return {
        'eta': eta,
        'eta_photon_1': eta,      # SHG: photon and energy conversion coincide
        ('E3_uJ' if pulsed else 'P3_W'): eta * in1,
        'I1_peak_MWcm2': beam.I_peak_Wm2 * 1e-10,
        'deff_pmV': deff_pmV, 'dk_radum': dk, 'K1': K,
        'model_ratios': ratios,
    }
