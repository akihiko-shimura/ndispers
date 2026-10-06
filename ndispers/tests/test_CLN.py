"""Undoped congruent LiNbO3 (Edwards & Lawrence 1984, Jundt 1997, Zelmon et
al. 1997): the transcriptions pinned to the sources and to each other."""
import warnings

import numpy as np
import pytest

import ndispers
import ndispers.media.crystals as C


@pytest.fixture(scope="module")
def ed():
    return C.CLN_Edwards1984()


@pytest.fixture(scope="module")
def ju():
    return C.CLN_Jundt1997()


@pytest.fixture(scope="module")
def ze():
    return C.CLN_Zelmon1997()


def n_o(x, wl, T):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return x.n(wl, 0.0, T, pol='o')


def n_e(x, wl, T):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return x.n(wl, np.pi / 2, T, pol='e')


# Edwards & Lawrence 1984, Table II: phase-matching temperature calculated
# from their Eq. (1) for noncritical type-I difference-frequency mixing of the
# 488 nm argon line (e) with a dye laser (o), by infrared wavelength (o).
@pytest.mark.parametrize("wl_ir,T_calc", [
    (2.159, 176.4), (2.249, 200.1), (2.337, 221.3), (2.432, 242.2), (2.541, 264.0),
    (2.643, 282.6), (2.771, 303.7), (2.901, 323.0), (3.049, 342.6), (3.235, 364.5)])
def test_edwards_reproduces_table_II(ed, wl_ir, T_calc):
    """One degree is 1e-4 in n_e at 488 nm, so this pins all fourteen
    coefficients, temperature terms included, for both rays. The residual
    0.2 to 0.5 degC is the pump wavelength (0.4880 here; with the line's
    vacuum wavelength 0.48813 it is 0.1 degC)."""
    wl_p = 0.4880
    wl_s = 1 / (1 / wl_p - 1 / wl_ir)
    dk = lambda T: n_e(ed, wl_p, T) / wl_p - n_o(ed, wl_s, T) / wl_s - n_o(ed, wl_ir, T) / wl_ir
    lo, hi = 0.0, 600.0
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if dk(lo) * dk(mid) <= 0 else (mid, hi)
    assert 0.5 * (lo + hi) == pytest.approx(T_calc, abs=1.0)


def test_zelmon_table_1_columns(ze):
    """Table 1 is labelled correctly (Table 2, the MgO-doped crystal of the
    same paper, is not - see MgOLN_Zelmon1997): n_o > n_e, and undoped
    LiNbO3 has the higher indices."""
    assert (n_o(ze, 1.064, 21), n_e(ze, 1.064, 21)) == pytest.approx((2.2321, 2.1555), abs=5e-5)
    mg = C.MgOLN_Zelmon1997()
    assert n_o(ze, 1.064, 21) > n_o(mg, 1.064, 21) and n_e(ze, 1.064, 21) > n_e(mg, 1.064, 21)


@pytest.mark.parametrize("wl", [0.532, 0.633, 1.064, 1.55, 2.0])
def test_three_sources_agree(ed, ju, ze, wl):
    """Edwards and Jundt fit the same room-temperature data (Nelson &
    Mikulyak, 24.5 degC); Zelmon measured another crystal at 21 degC. Observed
    spread 2e-4 (5e-4 for n_e between Jundt and Zelmon at 0.532 um)."""
    assert n_o(ed, wl, 21) == pytest.approx(n_o(ze, wl, 21), abs=3e-4)
    assert n_e(ed, wl, 24.5) == pytest.approx(n_e(ju, wl, 24.5), abs=3e-4)
    assert n_e(ju, wl, 21) == pytest.approx(n_e(ze, wl, 21), abs=6e-4)


def test_edwards_leaves_the_others_in_the_mid_infrared(ed, ju, ze):
    """Why its range stops at 3.4 um: one wl**2 term stands for the infrared
    absorption, and by 5 um it is 6e-3 above both infrared-corrected fits,
    which agree with each other to 1e-3 all the way."""
    assert n_e(ed, 5.0, 24.5) - n_e(ju, 5.0, 24.5) == pytest.approx(6.1e-3, abs=5e-4)
    assert n_o(ed, 5.0, 21) - n_o(ze, 5.0, 21) == pytest.approx(6.1e-3, abs=5e-4)
    for wl in (3.0, 4.0, 5.0):
        assert n_e(ju, wl, 21) == pytest.approx(n_e(ze, wl, 21), abs=1e-3)
    with pytest.warns(ndispers.ValidityWarning):
        ed.n(4.0, 0.0, 25, pol='o')


def test_temperature_dependence(ed, ju, ze):
    """Edwards and Jundt both fit Smith et al.'s temperature data at 0.633 um
    and must agree there; the ordinary index moves five times less."""
    args = (0.633, np.pi / 2, 25)
    assert ed.dndT(*args, pol='e') == pytest.approx(ju.dndT(*args, pol='e'), rel=0.02)
    assert ed.dndT(*args, pol='e') == pytest.approx(4.9e-5, rel=0.03)
    assert 0 < ed.dndT(0.633, 0.0, 25, pol='o') < 0.25 * ed.dndT(*args, pol='e')
    # both equations are referenced to 24.5 degC: the temperature terms vanish there
    a = ju.constants
    n2 = a['_a1_e'] + a['_a2_e'] / (1.064**2 - a['_a3_e']**2) + a['_a4_e'] / (1.064**2 - a['_a5_e']**2) - a['_a6_e'] * 1.064**2
    assert n_e(ju, 1.064, 24.5) == pytest.approx(np.sqrt(n2), rel=1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert ze.dndT(1.064, 0.3, 25, pol='e') == 0


def test_jundt_is_extraordinary_only(ju):
    assert ju.n(1.064, 0.3, 25) == ju.n(1.064, 0.3, 25, pol='e')      # default pol is 'e'
    with pytest.raises(ValueError):
        ju.n(1.064, 0.3, 25, pol='o')
    assert ju.d_sfg('d33', 1.064, 1.064, 25) == pytest.approx(25.2)


def test_nonlinear_coefficients_and_qpm(ed, ju, ze):
    """Shoji et al. 1997, Table 6 (congruent LiNbO3, 1.064 um SHG): d33 = 25.2,
    d31 = 4.6 pm/V. The 1.064 um SHG period for d33 agrees among the sources."""
    for x in (ed, ze):
        assert [x.d_sfg(il, 1.064, 1.064, 25) for il in ('d33', 'd31', 'd22')] == pytest.approx([25.2, 4.6, -2.1])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        periods = [x.qpm_period_sfg(1.064, 1.064, np.pi / 2, 25, 'e', 'e', 'e') for x in (ju, ed, ze)]
    assert periods == pytest.approx([6.8] * 3, abs=0.03)
