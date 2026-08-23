"""DFG / OPA / OPO: the SFG interaction with the pump given instead of found."""
import numpy as np
import pytest

import ndispers.media.crystals as C

NL = [C.BetaBBO_Eimerl1987, C.LBO_KK2018_xy, C.LBO_KK2018_zx, C.KTP_xy, C.KTP_zx,
      C.CLBO, C.KDP, C.BiBO_Miyata2009_yz, C.MgOLN_Zelmon1997, C.AGS_Kato1996, C.ZGP_Das2003]


def test_idler_energy_conservation_and_domain():
    b = C.BetaBBO_Eimerl1987()
    assert abs(1 / b.wl_idler(0.355, 0.5) - (1 / 0.355 - 1 / 0.5)) < 1e-15
    np.testing.assert_allclose(b.wl_idler(1.064, np.array([1.5, 2.0])), [3.66055, 2.27350], rtol=1e-5)
    for wl_s in (0.355, 0.3):                       # signal at or above the pump energy
        with pytest.raises(ValueError):
            b.wl_idler(0.355, wl_s)


@pytest.mark.parametrize("cls", NL, ids=lambda c: c.__name__)
def test_dfg_is_sfg_at_signal_and_idler(cls):
    """Every DFG method is the SFG one at (wl_s, wl_i): same angles, same dk,
    same d_eff, same QPM period. This is the whole content of stage 1."""
    x = cls()
    wl_p, wl_s, T = 1.064, 1.55, 25
    if isinstance(x, (C.AGS_Kato1996, C.ZGP_Das2003)):
        wl_p, wl_s = 2.05, 3.0                          # mid-infrared pump
    wl_i = x.wl_idler(wl_p, wl_s)
    a = x.pmAngles_dfg(wl_p, wl_s, T, deg=True)
    b = x.pmAngles_sfg(wl_s, wl_i, T, deg=True)
    assert a.pop('wl_i') == wl_i and abs(b.pop('wl3') - wl_p) < 1e-12
    assert a == b
    ang = 0.7
    for pols in ("ooe", "eeo", "oee", "eoe", "eoo", "oeo"):
        try:
            d1 = x.dk_dfg(wl_p, wl_s, ang, T, *pols)
        except ValueError:
            continue
        assert d1 == x.dk_sfg(wl_s, wl_i, ang, T, *pols)
        assert x.pmFactor_dfg(wl_p, wl_s, ang, T, *pols, 5.0) == x.pmFactor_sfg(wl_s, wl_i, ang, T, *pols, 5.0)
        assert x.qpm_period_dfg(wl_p, wl_s, ang, T, *pols) == x.qpm_period_sfg(wl_s, wl_i, ang, T, *pols)
        if hasattr(x, "deff_sfg"):
            th, ph = (ang, 0.3) if x.theta_rad == "var" else (None, ang)
            if x.theta_rad == "var" and x.phi_rad != "arb":
                ph = None
            assert x.deff_dfg(wl_p, wl_s, th, ph, T, *pols) == x.deff_sfg(wl_s, wl_i, th, ph, T, *pols)
    for il in getattr(x, "_d_ref", {}):
        assert x.d_dfg(il, wl_p, wl_s, T) == x.d_sfg(il, wl_s, wl_i, T)


@pytest.mark.parametrize("cls,wl_p,wl_s,pols", [
    (C.BetaBBO_Eimerl1987, 0.355, 0.5, "ooe"),
    (C.BetaBBO_Eimerl1987, 0.355, 0.5, "oee"),
    (C.BetaBBO_Eimerl1987, 0.355, 0.5, "eoe"),
    (C.LBO_KK2018_xy, 0.355, 0.5, "ooe"),
    (C.KTP_zx, 1.064, 1.55, "eeo"),
    (C.KTP_zx, 1.064, 1.55, "eoo"),
    (C.KTP_xy, 1.064, 1.55, "eoe"),
    (C.BiBO_Miyata2009_yz, 0.4, 0.6, "eeo"),
    (C.BiBO_Miyata2009_yz, 0.4, 0.6, "oeo"),
    (C.LBO_KK2018_zx, 0.355, 0.5, "ooe"),
], ids=str)
def test_tuning_round_trip(cls, wl_p, wl_s, pols):
    """At the angle that phase-matches (wl_p, wl_s), the tuning curve passes
    through wl_s - the wavelength root finder and the angle root finder agree."""
    x = cls()
    pm = x.pmAngles_dfg(wl_p, wl_s, 25, tol_deg=1e-6)
    key = "theta" if x.theta_rad == "var" else "phi"
    angles = pm[pols][key]
    assert angles, f"{cls.__name__} {pols}: no phase matching to round-trip"
    pairs = x.tuning_dfg(wl_p, angles[0], 25, *pols, wl_i_max=5.0)
    # the angle's residual dk maps to a wavelength offset through the (often
    # small) slope of the tuning curve; 10 pm covers every case here
    assert any(abs(ws - wl_s) < 1e-5 for ws, wi in pairs), pairs
    for ws, wi in pairs:
        assert ws <= wi + 1e-9                       # signal side of degeneracy only
        assert abs(1 / ws + 1 / wi - 1 / wl_p) < 1e-12


def test_tuning_type_ii_branches_are_the_swapped_triples():
    """Type II: the branch with the e-polarized wave as signal is 'eoe', with it
    as idler 'oee'; each triple returns only the wl_s <= wl_i side, so together
    they cover the curve once."""
    b = C.BetaBBO_Eimerl1987()
    th = np.radians(35)
    a = b.tuning_dfg(0.355, th, 25, "o", "e", "e", wl_i_max=3.0)
    c = b.tuning_dfg(0.355, th, 25, "e", "o", "e", wl_i_max=3.0)
    assert a and not c or c and not a or (a and c)  # whichever branch exists at this angle
    for ws, wi in a + c:
        assert ws <= wi


def test_tuning_degenerate_point_is_found():
    """At the Type I SHG angle of 532 -> 1064 nm, the 532-nm-pumped OPO is
    degenerate: the curve touches 1064 + 1064 nm. The flat top means a tiny
    residual in the angle moves the root a lot, so only ask for the pair to be
    within the grid resolution of degeneracy."""
    b = C.BetaBBO_Eimerl1987()
    th = b.pmAngles_sfg(1.064, 1.064, 25)["ooe"]["theta"][0]
    pairs = b.tuning_dfg(0.532, th, 25, "o", "o", "e", wl_i_max=3.0)
    assert pairs and abs(pairs[0][0] - 1.064) < 0.01


def test_tuning_no_solution_cases():
    b = C.BetaBBO_Eimerl1987()
    assert b.tuning_dfg(0.355, np.radians(5), 25, "o", "o", "e", wl_i_max=3.0) == []   # far from PM
    assert C.SLN().tuning_dfg(1.064, np.pi / 2, 25, "o", "o", "e") == []               # no o-ray
    assert b.tuning_dfg(0.355, 0.5, 25, "o", "o", "e", wl_i_max=0.5) == []             # empty range
