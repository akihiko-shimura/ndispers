"""
Second-order nonlinearity: Miller-rule scaling and d_eff (ndispers.groups).

Theory and conventions are in docs/dev/deff_theory.tex. These tests pin
(1) the round trip of Miller's delta at the reference point, (2) that every
crystal with nonlinear data evaluates d_eff for every polarization triple
without raising, (3) that the generic tensor contraction reproduces the
published closed forms for each point group, and (4) bookkeeping between a
crystal's docstring and the group it inherits.
"""
import re

import numpy as np
import pytest

import ndispers.media.crystals as C
from ndispers.groups import NonlinearGroup, Uniax_3m, Uniax_32, Uniax_42m, Uniax_4mm

NAMES = sorted(n for n in dir(C) if isinstance(getattr(C, n), type)
               and issubclass(getattr(C, n), NonlinearGroup))
T = 25.0


def can_o(x):
    try:
        x.n_expr('o')
    except ValueError:
        return False
    return True


# ---------------------------------------------------------------------------
# 1. Miller's rule reproduces the reference value at the reference wavelengths
@pytest.mark.parametrize("name", NAMES)
def test_reference_round_trip(name):
    x = getattr(C, name)()
    for il, (d0, wl1, wl2) in x._d_ref.items():
        assert x.d_sfg(il, wl1, wl2, T) == pytest.approx(d0, rel=1e-12), (name, il)
        # the generated alias
        assert getattr(x, f"{il}_sfg")(wl1, wl2, T) == pytest.approx(d0, rel=1e-12)
    assert x._d_note, f"{name}: _d_note should say where the values come from"


def test_unknown_component_raises():
    with pytest.raises(KeyError):
        C.KDP().d_sfg("d22", 1.064, 1.064, T)
    with pytest.raises(AttributeError):
        C.KDP().d22_sfg


# ---------------------------------------------------------------------------
# 2. d_eff evaluates for every crystal and polarization triple
POLS = ["ooe", "oee", "eoe", "eeo", "oeo", "eoo", "ooo", "eee"]


@pytest.mark.parametrize("name", NAMES)
@pytest.mark.parametrize("pols", POLS)
def test_deff_smoke(name, pols):
    x = getattr(C, name)()
    if 'o' in pols and not can_o(x):
        pytest.skip("no o-ray Sellmeier equation")
    th = 0.6 if x.theta_rad == 'var' else None
    ph = 0.4 if x.phi_rad in ('var', 'arb') else None
    for wl1, wl2 in [(1.064, 1.064), (0.8, 0.4)]:
        try:
            v = x.deff_sfg(wl1, wl2, th, ph, T, *pols)
        except ValueError as e:
            # a component the crystal has no reference value for
            assert "no reference value" in str(e), (name, pols, e)
            continue
        assert np.isfinite(v), (name, pols, v)
    # array angle in, array out, matching the scalar values
    if x.theta_rad == 'var':
        ths = np.array([0.3, 0.6, 0.9])
        try:
            arr = x.deff_sfg(1.064, 1.064, ths, ph, T, *pols)
        except ValueError:
            return
        assert arr.shape == ths.shape
        assert arr[1] == pytest.approx(x.deff_sfg(1.064, 1.064, 0.6, ph, T, *pols))


def test_missing_component_is_loud():
    # LB4 holds d31 only; eee would need d33
    with pytest.raises(ValueError, match="d33"):
        C.LB4().deff_sfg(1.064, 1.064, 0.5, 0.3, T, 'e', 'e', 'e')
    # ...while ooe, which does not use d33, is fine
    assert np.isfinite(C.LB4().deff_sfg(1.064, 1.064, 0.5, 0.3, T, 'o', 'o', 'e'))


# ---------------------------------------------------------------------------
# 3. The generic contraction reproduces the published closed forms
# (Dmitriev, Gurzadyan & Nikogosyan 1999, Sec. 2.1; Roberts 1992; the LBO
# forms also appear on the Castech data sheet). Theta_m = theta + rho_m with
# rho_m the walk-off of wave m; overall sign is a convention, so |.| compared.
WL1, WL2 = 1.064, 0.8
WL3 = 1 / (1 / WL1 + 1 / WL2)


def Th(x, wl, th, pol):
    return th + (x.woa_theta(wl, th, T, pol='e') if pol == 'e' else 0.)


def Ph(x, wl, ph, pol):
    return ph + (x.woa_phi(wl, ph, T, pol='e') if pol == 'e' else 0.)


def d(x, il):
    return x.d_sfg(il, WL1, WL2, T)


@pytest.mark.parametrize("th,ph", [(0.4, 0.3), (1.1, 1.9)])
def test_closed_form_3m(th, ph):
    x = C.BetaBBO_Eimerl1987()
    d22, d31 = d(x, "d22"), d(x, "d31")
    T1, T2, T3 = Th(x, WL1, th, 'e'), Th(x, WL2, th, 'e'), Th(x, WL3, th, 'e')
    ref = {
        "ooe": d31 * np.sin(T3) - d22 * np.cos(T3) * np.sin(3 * ph),
        "oee": d22 * np.cos(T2) * np.cos(T3) * np.cos(3 * ph),
        "eoe": d22 * np.cos(T1) * np.cos(T3) * np.cos(3 * ph),
    }
    for pols, r in ref.items():
        assert abs(x.deff_sfg(WL1, WL2, th, ph, T, *pols)) == pytest.approx(abs(r), rel=1e-12), pols


@pytest.mark.parametrize("th,ph", [(0.4, 0.3), (1.1, 1.9)])
def test_closed_form_32(th, ph):
    x = C.KBBF()
    d11 = d(x, "d11")
    T1, T2, T3 = Th(x, WL1, th, 'e'), Th(x, WL2, th, 'e'), Th(x, WL3, th, 'e')
    ref = {
        "ooe": d11 * np.cos(T3) * np.cos(3 * ph),
        "oee": d11 * np.cos(T2) * np.cos(T3) * np.sin(3 * ph),
        "eoe": d11 * np.cos(T1) * np.cos(T3) * np.sin(3 * ph),
    }
    for pols, r in ref.items():
        assert abs(x.deff_sfg(WL1, WL2, th, ph, T, *pols)) == pytest.approx(abs(r), rel=1e-12), pols


@pytest.mark.parametrize("th,ph", [(0.4, 0.3), (1.1, 1.9)])
def test_closed_form_42m(th, ph):
    x = C.KDP()
    d36 = d(x, "d36")
    T1, T2, T3 = Th(x, WL1, th, 'e'), Th(x, WL2, th, 'e'), Th(x, WL3, th, 'e')
    ref = {
        "ooe": d36 * np.sin(T3) * np.sin(2 * ph),
        "oee": d36 * np.sin(T2 + T3) * np.cos(2 * ph),
        "eoe": d36 * np.sin(T1 + T3) * np.cos(2 * ph),
    }
    for pols, r in ref.items():
        assert abs(x.deff_sfg(WL1, WL2, th, ph, T, *pols)) == pytest.approx(abs(r), rel=1e-12), pols


@pytest.mark.parametrize("th,ph", [(0.4, 0.3), (1.1, 1.9)])
def test_closed_form_4mm(th, ph):
    x = C.LB4()
    d31 = d(x, "d31")
    T3 = Th(x, WL3, th, 'e')
    assert abs(x.deff_sfg(WL1, WL2, th, ph, T, 'o', 'o', 'e')) == pytest.approx(abs(d31 * np.sin(T3)), rel=1e-12)
    # type II vanishes identically under Kleinman symmetry
    assert x.deff_sfg(WL1, WL2, th, ph, T, 'o', 'e', 'e') == pytest.approx(0., abs=1e-15)


# ---------------------------------------------------------------------------
# 4. Bookkeeping: the group a crystal inherits matches its docstring
GROUP_OF = {"3m": Uniax_3m, "32": Uniax_32, "4̄2m": Uniax_42m, "4mm": Uniax_4mm}


@pytest.mark.parametrize("name", NAMES)
def test_group_matches_point_group(name):
    cls = getattr(C, name)
    pg = re.search(r'Point group\s*:\s*(\S+)', cls.__doc__).group(1)
    assert issubclass(cls, GROUP_OF[pg]), f"{name}: point group {pg} but inherits {cls.__mro__[1:3]}"


