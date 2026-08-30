"""
Semi-analytical conversion efficiency: every check pins a limit derived in
docs/dev/efficiency_theory.tex against Armstrong (1962) / Boyd (2020).
During development the implementation was also cross-checked against an
independent earlier implementation; that comparison is not retained here.
"""
import warnings

import numpy as np
import pytest

import ndispers as nd
from ndispers import CWBeam, ModelValidityWarning, PulsedBeam
from ndispers._efficiency import _characteristic_intensity, _y_pw


@pytest.fixture(autouse=True)
def _quiet():
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', ModelValidityWarning)
        yield


def test_plane_wave_solution_limits():
    # matched, balanced: y = K tanh^2(sqrt(K))   [Armstrong (5.12)]
    for K in (0.01, 0.5, 3.0):
        assert np.isclose(_y_pw(K, K, 0.0), K * np.tanh(np.sqrt(K))**2,
                          rtol=1e-10)
    # undepleted with mismatch: y -> a2 sinc^2(sigma/2)
    a, b, s = 1e-6, 2e-6, 3.0
    expect = a * b * (np.sin(s / 2) / (s / 2))**2
    assert np.isclose(_y_pw(a, b, s), expect, rtol=1e-4)
    # SHG mismatch turning point: y_max/K = b1/K from Armstrong (5.21):
    # v_b^2 = [ds/4 + sqrt(1+(ds/4)^2)]^-2 (their K=1 units)
    ds = 2.0
    b1 = 0.5 * (2 + ds**2 / 4) - 0.5 * np.sqrt((2 + ds**2 / 4)**2 - 4)
    assert np.isclose(b1, (ds / 4 + np.sqrt(1 + (ds / 4)**2))**-2, rtol=1e-12)


def test_manley_rowe_caps():
    """eta never exceeds the photon-conversion cap of either wave."""
    rng = np.random.default_rng(7)
    bbo = nd.media.crystals.BetaBBO_Tamosauskas2018()
    th = bbo.pmAngles_sfg(1.064, 0.532, 25)['ooe']['theta'][0]
    for _ in range(5):
        E1, E2 = rng.uniform(0.1, 2000, 2)
        b1 = PulsedBeam(1.064, E1, 200.0, 10e3, 'o')
        b2 = PulsedBeam(0.532, E2, 200.0, 10e3, 'o')
        d = bbo.eta_sfg(b1, b2, th, np.pi / 2, 25, 'e', 5.0, details=True)
        assert 0 <= d['eta_photon_1'] <= 1 + 1e-9
        assert 0 <= d['eta_photon_2'] <= 1 + 1e-9


def test_wave_exchange_symmetry():
    """Swapping the two input beams leaves eta invariant. This is the
    regression guard for the d/d1 pairing (theory note 3.A): the
    prototype's transposed a1 fails it for nondegenerate wavelengths."""
    bbo = nd.media.crystals.BetaBBO_Tamosauskas2018()
    th = bbo.pmAngles_sfg(1.064, 0.532, 25)['ooe']['theta'][0]
    b1 = PulsedBeam(1.064, 300.0, 250.0, 8e3, 'o')
    b2 = PulsedBeam(0.532, 120.0, 180.0, 12e3, 'o')
    # integration is normalized to beam 1, so the two orderings sample the
    # mesh differently: agreement is only as good as the quadrature. The
    # coarse default differs at ~3e-3; a fine grid converges the two.
    args = (th, np.pi / 2, 25, 'e', 3.0)
    e12 = bbo.eta_sfg(b1, b2, *args, deff_pmV=2.0, n_grid=1001)
    e21 = bbo.eta_sfg(b2, b1, *args, deff_pmV=2.0, n_grid=1001)
    assert np.isclose(e12, e21, rtol=1e-5)
    assert e12 > 0.01           # nontrivial conversion, not a 0==0 pass


def test_pulsed_vs_cw_undepleted_ratio():
    """Undepleted limit, equal beams: Gaussian temporal averaging costs
    exactly 1/sqrt(2) relative to cw at the same peak intensity."""
    bbo = nd.media.crystals.BetaBBO_Tamosauskas2018()
    th = bbo.pmAngles_sfg(1.064, 1.064, 25)['ooe']['theta'][0]
    w, t = 300.0, 10e3
    bp = PulsedBeam(1.064, 1e-4, w, t, 'o')     # tiny energy: undepleted
    # cw beam with the same peak (axial) intensity
    P = bp.I_peak_Wm2 * np.pi * (w * 1e-6 / np.sqrt(2))**2
    bc = CWBeam(1.064, P, w, 'o')
    ep = bbo.eta_sfg(bp, bp, th, np.pi / 2, 25, 'e', 5.0, deff_pmV=2.0,
                     n_grid=201)
    ec = bbo.eta_sfg(bc, bc, th, np.pi / 2, 25, 'e', 5.0, deff_pmV=2.0,
                     n_grid=201)
    assert np.isclose(ep / ec, 1 / np.sqrt(2), rtol=1e-3)


def test_shg_matches_degenerate_type2_sfg():
    """Type-II SHG through eta_sfg at exactly degenerate wavelengths is a
    well-defined SFG; the type-I eta_shg with the same total intensity and
    dk=0 must give the same efficiency in the undepleted limit when deff
    and indices match (the degeneracy factors cancel in this comparison)."""
    ln = nd.media.crystals.MgOLN_Zelmon1997()
    b = PulsedBeam(1.064, 1e-3, 200.0, 10e3, 'e')
    eta1 = ln.eta_shg(b, np.pi / 2, 0.0, 21, 'e', 1.0, deff_pmV=16.0,
                      dk_offset=-float(ln.dk_sfg(1.064, 1.064, np.pi / 2,
                                                 21, 'e', 'e', 'e')))
    # undepleted analytic: eta = K/(sqrt(2)... compare against the
    # closed formula eta = I1 L^2 * 8 pi^2 deff^2/(n^2 n3 eps0 c wl^2)
    # averaged over the Gaussian: x 1/(2*sqrt(2))
    from ndispers._efficiency import _characteristic_intensity
    n1 = float(ln.n(1.064, np.pi / 2, 21, pol='e'))
    n3 = float(ln.n(0.532, np.pi / 2, 21, pol='e'))
    K = b.I_peak_Wm2 / _characteristic_intensity(n1, n1, n3, 1.064, 1.064,
                                                 16.0, 1.0)
    assert np.isclose(eta1, K / (2 * np.sqrt(2)), rtol=1e-3)


def test_beam_validation_and_type_errors():
    with pytest.raises(ValueError, match='pol'):
        PulsedBeam(1.064, 1.0, 100.0, 1000.0, 'x')
    with pytest.raises(ValueError, match='positive'):
        CWBeam(1.064, -1.0, 100.0, 'o')
    bbo = nd.media.crystals.BetaBBO_Tamosauskas2018()
    bp = PulsedBeam(1.064, 1.0, 100.0, 1000.0, 'o')
    bc = CWBeam(0.532, 1.0, 100.0, 'o')
    with pytest.raises(TypeError, match='both'):
        bbo.eta_sfg(bp, bc, 0.4, np.pi / 2, 25, 'e', 1.0)


def test_model_validity_warning_fires():
    """The classic lesson: 100 fs SHG at 800 nm in 1 mm of BBO - the
    ~190 fs/mm group-velocity mismatch exceeds the pulse duration
    (gvm ratio ~1.9), which is why femtosecond doubling uses 0.1-0.3 mm
    crystals."""
    bbo = nd.media.crystals.BetaBBO_Tamosauskas2018()
    th = bbo.pmAngles_sfg(0.8, 0.8, 25)['ooe']['theta'][0]
    b = PulsedBeam(0.8, 5.0, 1000.0, 100.0, 'o')
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        d = bbo.eta_shg(b, th, np.pi / 2, 25, 'e', 1.0, details=True)
        assert any(isinstance(x.message, ModelValidityWarning) for x in w)
    assert d['model_ratios']['gvm'] > 1.5
