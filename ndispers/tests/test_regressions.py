"""Regression tests for bugs fixed in _baseclass.py and the KTP / SLN modules."""
import pickle
import re

import numpy as np

import pytest

import ndispers.media.crystals as C

PLANES = ['KTP_xy', 'KTP_yz', 'KTP_zx']


@pytest.mark.parametrize("cls", PLANES)
@pytest.mark.parametrize("pol", ['o', 'e'])
def test_ktp_all_planes_and_pols(cls, pol):
    """_dndT_z was clobbered by a duplicate _dndT_y assignment, so half of these
    raised AttributeError and the rest silently used the z coefficient for n_y."""
    n = getattr(C, cls)().n(1.064, 0.3, 100, pol=pol)
    assert 1.7 < n < 1.9


# o-wave of each principal plane is one principal index: xy->z, yz->x, zx->y
@pytest.mark.parametrize("cls,n_lit", [('KTP_xy', 1.8302), ('KTP_yz', 1.7377), ('KTP_zx', 1.7453)])
def test_ktp_principal_indices(cls, n_lit):
    """The Sellmeier D term had a flipped sign, putting every index 0.17-0.61 high."""
    assert getattr(C, cls)().n(1.064, 0.3, 20, pol='o') == pytest.approx(n_lit, abs=2e-3)


def test_ktp_type2_shg_angle():
    """Type-II SHG of 1.064 um in the xy plane is near phi = 23.3 deg."""
    phi = C.KTP_xy().pmAngles_sfg(1.064, 1.064, 20, deg=True)['oee']['phi']
    assert phi and phi[0] == pytest.approx(23.5, abs=0.5)


def test_ktp_dndT_y_is_not_dndT_z():
    ktp = C.KTP_xy()
    # y and z thermo-optic coefficients are different expressions in Kato 2002
    assert ktp._dndT_y != ktp._dndT_z


# o-wave of each plane again: xy->z, yz->x, zx->y. Kato 2002 gives dn_x < dn_y < dn_z.
# The index tests above all sit at T=20 degC, where dndT_i*(T-20) vanishes, so this
# is what actually pins the _dndT_y/_dndT_z fix.
@pytest.mark.parametrize("cls,dndT_lit", [
    ('KTP_yz', 6.2338e-6), ('KTP_zx', 8.3379e-6), ('KTP_xy', 1.4418e-5)])
def test_ktp_dndT_per_axis(cls, dndT_lit):
    assert getattr(C, cls)().dndT(1.064, 0.3, 20, pol='o') == pytest.approx(dndT_lit, rel=1e-3)


def test_ktp_zx_has_variable_theta():
    """theta_rad was 'arb', so pmAngles_sfg took neither branch and returned {}.
    Whether a given interaction has a solution is a separate question - what this
    checks is that the angle is reported under 'theta'."""
    zx = C.KTP_zx()
    assert zx.theta_rad == 'var'
    assert zx.pmAngles_sfg(1.064, 1.064, 25, deg=True)['ooe'] == {'theta': [], 'phi': None}


def test_sln_is_callable():
    """SLN assigned only _o constants while n_expr accepts only 'e'."""
    # value at the reference temperature, cross-checked against
    # refractiveindex.info Gayer-1-e (1% MgO stoichiometric LiNbO3, e-ray)
    assert C.SLN().n(1.064, 0.3, 24.5, pol='e') == pytest.approx(2.1432, abs=1e-4)


@pytest.mark.parametrize("T_degC,wl_nm", [(40, 1544), (200, 1571)])
def test_sln_qpm_matches_gayer_figure4(T_degC, wl_nm):
    """b1 was the pre-erratum 4.677e-7. Only the corrected 4.677e-6 reproduces
    Fig. 4 of Gayer 2008 (SHG in a 19.36 um period crystal); the old value drifts
    to 1614 nm at 200 degC. Room-temperature n is nearly insensitive to b1, so
    phase matching is what pins it."""
    from ndispers.helper import brentq
    sln = C.SLN()
    period = 19.36 * (1 + 1.54e-5 * (T_degC - 25))       # a-axis thermal expansion

    def mismatch(wl_f):
        n = lambda w: float(sln.n(w, 0.3, T_degC, pol='e'))
        return wl_f / (2 * (n(wl_f / 2) - n(wl_f))) - period

    assert brentq(mismatch, 1.30, 1.90) * 1000 == pytest.approx(wl_nm, abs=3)


@pytest.mark.parametrize("cls", ['CLBO', 'KTP_zx', 'SLN', 'BetaBBO_Eimerl1987'])
def test_picklable_after_use(cls):
    """lambdify caches made used instances unpicklable, breaking multiprocessing."""
    x = getattr(C, cls)()
    before = x.n(1.064, 0.3, 100, pol='e')      # populate the cache
    y = pickle.loads(pickle.dumps(x))
    assert y.n(1.064, 0.3, 100, pol='e') == before
    assert (y.plane, y.theta_rad, y.phi_rad) == (x.plane, x.theta_rad, x.phi_rad)


# ndispers ships several independent parameterisations of the same crystal, so they
# can check each other without any external reference data. This is what caught the
# _betaBBO_KK2010 dn/dT bug (_I used twice, _H never), and it is the only automated
# defence against a transcription error that returns a plausible wrong number
# instead of raising - the class of bug the KTP Sellmeier sign error belonged to.
SAME_CRYSTAL = {
    'BBO': ['BetaBBO_Eimerl1987', 'BetaBBO_Ghosh1995', 'BetaBBO_KK2010',
            'BetaBBO_Tamosauskas2018'],
    **{f'LBO_{p}': [f'LBO_{s}_{p}' for s in
                    ('Castech', 'Ghosh1995', 'KK1994', 'KK2018', 'Newlight')]
       for p in ('xy', 'yz', 'zx')},
    # infrared crystals: the sets were taken from refractiveindex.info's
    # transcription of the papers, so two sources agreeing is the check that
    # the transcription (theirs and ours) is sound - evaluated at their own
    # wavelengths below
    'ZGP': ['ZGP_Zelmon2001', 'ZGP_Das2003'],
    'AGS': ['AGS_Kato1996', 'AGS_Takaoka1999'],
}
IR_WL = {'ZGP': [2.5, 3.0, 5.0, 8.0], 'AGS': [1.064, 2.0, 5.0, 10.0]}


def _spread(names, meth, wl, T, pol):
    v = [float(getattr(getattr(C, nm)(), meth)(wl, 0.3, T, pol=pol)) for nm in names]
    return max(v) - min(v)


@pytest.mark.parametrize("group", sorted(g for g in SAME_CRYSTAL if g not in IR_WL))
@pytest.mark.parametrize("wl,T", [(0.4, 20), (0.532, 20), (1.064, 20), (1.064, 100), (1.55, 20)])
@pytest.mark.parametrize("pol", ['o', 'e'])
def test_sources_agree_on_index(group, wl, T, pol):
    """Observed worst case is 1.6e-3; anything past 5e-3 is a transcription error,
    not source scatter (the KTP sign bug was 0.17-0.61)."""
    assert _spread(SAME_CRYSTAL[group], 'n', wl, T, pol) < 5e-3


@pytest.mark.parametrize("group", sorted(IR_WL))
@pytest.mark.parametrize("pol", ['o', 'e'])
def test_ir_sources_agree_on_index(group, pol):
    """ZGP: Zelmon 2001 vs Das 2003 differ by up to 5e-3 at 8 um (different
    samples, different ranges); AGS: Kato 1996 vs Takaoka 1999 by < 1e-3."""
    for wl in IR_WL[group]:
        assert _spread(SAME_CRYSTAL[group], 'n', wl, 20, pol) < 6e-3, (group, wl)


def test_bibo_reproduces_review_index_values():
    """Petrov et al. 2010, Table 3: n at 632.8 nm and 1064 nm computed from the
    Table 2 Sellmeier set. Pins the transcription of that set."""
    x = C.BiBO_Miyata2009_xy()
    nx = lambda wl: float(x._n_axis('x', wl, 20)); ny = lambda wl: float(x._n_axis('y', wl, 20)); nz = lambda wl: float(x._n_axis('z', wl, 20))
    assert (round(nx(0.6328), 5), round(ny(0.6328), 5), round(nz(0.6328), 5)) == (1.77668, 1.80641, 1.94582)
    assert (round(nx(1.064), 5), round(nz(1.064), 5)) == (1.75752, 1.91711)


def test_isotropic_spot_values():
    """YAG n_d = 1.8328 (Zelmon 1998 fit; Hrabovsky 2021 gives 1.8325); N-BK7
    n_d = 1.5168 is the catalog value; LiF, BaF2 and MgF2 at 0.5893 um against
    Tropf's Table 20 n_infinity-style listings (1.388 / 1.4744 / 1.3734(o)) and the
    textbook MgF2 birefringence (+0.012)."""
    import ndispers.media.glasses as G
    assert G.YAG().n(0.5876, 20) == pytest.approx(1.8328, abs=1e-4)
    assert G.NBK7().n(0.5876, 20) == pytest.approx(1.5168, abs=1e-4)
    assert G.LiF().n(0.5893, 20) == pytest.approx(1.392, abs=2e-3)
    assert G.BaF2().n(0.5893, 20) == pytest.approx(1.474, abs=2e-3)
    m = C.MgF2()
    assert m.n(0.5893, 0, 20, pol='o') == pytest.approx(1.378, abs=1e-3)
    assert m.n(0.5893, np.pi / 2, 20, pol='e') - m.n(0.5893, 0, 20, pol='o') == pytest.approx(0.0119, abs=5e-4)
    # the rest of the isotropic set, against textbook values
    assert G.ZnSe().n(10.6, 20) == pytest.approx(2.403, abs=1e-3)
    assert G.ZnS().n(1.064, 20) == pytest.approx(2.288, abs=2e-3)
    assert G.Si().n(2.0, 20) == pytest.approx(3.453, abs=2e-3)
    assert G.Ge().n(4.0, 20) == pytest.approx(4.025, abs=2e-3)
    assert G.Diamond().n(0.5876, 20) == pytest.approx(2.4175, abs=5e-4)
    assert G.SF10().n(0.5876, 20) == pytest.approx(1.72825, abs=1e-4)
    assert G.SF11().n(0.5876, 20) == pytest.approx(1.78472, abs=1e-4)
    assert G.SF57().n(0.5876, 20) == pytest.approx(1.84666, abs=1e-4)
    # YVO4: strongly positive uniaxial; LiIO3 negative
    y = C.YVO4()
    assert y.n(1.064, np.pi / 2, 20, pol='e') - y.n(1.064, 0, 20, pol='o') == pytest.approx(0.208, abs=3e-3)
    li = C.LiIO3()
    assert (li.n(1.064, 0, 20, pol='o'), li.n(1.064, np.pi / 2, 20, pol='e')) == pytest.approx((1.8559, 1.7164), abs=1e-3)


def test_qpm_period():
    """First-order period for 1064 nm SHG in 5% MgO:LiNbO3 along x, eee (d33):
    Lambda = lambda / (2 (n_2w - n_w)) for the extraordinary index; and the
    relation to dk_sfg."""
    x = C.MgOLN_Zelmon1997()
    ne = lambda wl: x.n(wl, np.pi / 2, 20, pol='e')
    expect = 1.064 / (2 * (ne(0.532) - ne(1.064)))
    got = x.qpm_period_sfg(1.064, 1.064, np.pi / 2, 20, 'e', 'e', 'e')
    assert got == pytest.approx(expect, rel=1e-12)
    assert 6.0 < got < 7.5                                   # the familiar ~6.5-7 um
    assert x.qpm_period_sfg(1.064, 1.064, np.pi / 2, 20, 'e', 'e', 'e', order=3) == pytest.approx(3 * got)
    # the documented trap: at theta = 0 the e-ray degenerates to the o-ray, so
    # the same call returns a period computed from n_o
    assert x.qpm_period_sfg(1.064, 1.064, 0.0, 20, 'e', 'e', 'e') == pytest.approx(5.91, abs=0.02)


@pytest.mark.parametrize("wl", [0.532, 1.064])
@pytest.mark.parametrize("pol", ['o', 'e'])
def test_bbo_sources_agree_on_dndT(wl, pol):
    """All four BBO sources agree to 5e-7 once _H is used. Note this check is only
    meaningful for BBO: the five LBO sources genuinely disagree on dn/dT by up to
    1e-5 (vendor values vs Kato's fits), so no useful bound exists for them and a
    KK2010-sized error would hide inside that scatter."""
    assert _spread(SAME_CRYSTAL['BBO'], 'dndT', wl, 20, pol) < 2e-6


def test_betaBBO_KK2010_uses_H_coefficient():
    """dndT_o_expr had _I_o in both the 1/wl and the constant slot, leaving _H_o
    unused and dn_o/dT 1.76x too large."""
    x = C.BetaBBO_KK2010()
    assert x.dndT(1.064, 0.3, 20, pol='o') == pytest.approx(-1.61e-5, rel=5e-3)


@pytest.mark.parametrize("meth", ['n', 'dn_wl', 'd2n_wl', 'd3n_wl', 'GD', 'GV', 'ng',
                                  'GVD', 'TOD', 'dndT'])
def test_kdp_methods_callable(meth):
    """Every KDP wrapper but n and TOD declared T_deg and passed T_degC."""
    getattr(C.KDP(), meth)(0.532, 0.3, 20, pol='o')


# Zernike's measured indices, the data Ghosh's Sellmeier coefficients were
# fitted to. The previous implementation dropped the wl**2 factor of the
# lattice term and flipped its sign, putting n_o 0.023 high at 1.064 um.
@pytest.mark.parametrize("wl_um,n_o,n_e", [
    (0.4047, 1.52449, 1.47787), (0.5461, 1.51152, 1.47044),
    (0.6943, 1.50529, 1.46685), (1.0642, 1.49384, 1.46041)])
def test_kdp_indices(wl_um, n_o, n_e):
    kdp = C.KDP()
    assert kdp.n(wl_um, 0, 24.8, pol='o') == pytest.approx(n_o, abs=2e-3)
    assert kdp.n(wl_um, np.pi/2, 24.8, pol='e') == pytest.approx(n_e, abs=2e-3)


def test_kdp_thermo_optic_is_populated():
    """The coefficients were zero placeholders until the Ghosh 1992 values were
    entered; KDP is a negative thermo-optic crystal with |dn_o/dT| > |dn_e/dT|."""
    kdp = C.KDP()
    dn_o = kdp.dndT(1.0642, 0, 24.8, pol='o')
    dn_e = kdp.dndT(1.0642, np.pi/2, 24.8, pol='e')
    assert dn_o == pytest.approx(-3.72e-5, rel=2e-2)
    assert dn_e == pytest.approx(-2.32e-5, rel=2e-2)
    assert dn_o < dn_e < 0


def test_dndT2_runs_and_is_nonzero_for_e_ray():
    """dndT2 raised KeyError (unseeded cache), then TypeError (no phi injection)."""
    clbo = C.CLBO()
    assert clbo.dndT2(0.532, 0.5, 100.0, pol='o') == 0     # n_o is linear in T
    assert clbo.dndT2(0.532, 0.5, 100.0, pol='e') != 0     # n_e(theta) is not


ALL_MEDIA = [nm for nm in sorted(dir(C))
             if not nm.startswith('_') and isinstance(getattr(C, nm), type)]


# A point group belongs to exactly one crystal system, and only a
# non-centrosymmetric one can have a second-order nonlinearity at all. Both
# facts are checkable against the docstring header.
#
# Note what this does not catch: alpha-BBO carried beta-BBO's point group (3m
# where it should be 3̄m), and since both are trigonal no consistency check
# would have flagged it - that one needed reading the literature.
POINT_GROUPS = {                    # point group: (crystal system, centrosymmetric)
    '3m': ('trigonal', False), '3̄m': ('trigonal', True), '32': ('trigonal', False),
    '4̄2m': ('tetragonal', False), '4mm': ('tetragonal', False), '4/mmm': ('tetragonal', True), '4/m': ('tetragonal', True),
    '6': ('hexagonal', False), '4̄3m': ('cubic', False),
    'mm2': ('orthorhombic', False), '2': ('monoclinic', False), 'm3̄m': ('cubic', True),
}


@pytest.mark.parametrize("name", ALL_MEDIA)
def test_point_group_matches_crystal_system(name):
    import re
    doc = getattr(C, name).__doc__
    pg = re.search(r'Point group\s*:\s*(\S+)', doc)
    cs = re.search(r'Crystal system\s*:\s*(\w+)', doc)
    assert pg and cs, f"{name}: header is missing point group or crystal system"
    entry = POINT_GROUPS.get(pg.group(1))
    assert entry, f"{name}: unrecognised point group {pg.group(1)!r}"
    system, _ = entry
    assert cs.group(1).lower() == system, (
        f"{name}: crystal system {cs.group(1)!r} contradicts point group {pg.group(1)!r}")


@pytest.mark.parametrize("name", ALL_MEDIA)
def test_centrosymmetric_media_have_no_nonlinear_coefficients(name):
    """Second-order nonlinearity vanishes in a centrosymmetric crystal, so a
    medium carrying d coefficients must not be labelled with a centrosymmetric
    point group, and vice versa."""
    import re
    x = getattr(C, name)()
    pg = re.search(r'Point group\s*:\s*(\S+)', type(x).__doc__).group(1)
    _, centrosymmetric = POINT_GROUPS[pg]
    has_d = any('shg' in k for k in x.constants)
    assert not (centrosymmetric and has_d), (
        f"{name}: point group {pg} is centrosymmetric but the class defines "
        f"second-harmonic coefficients")


@pytest.mark.parametrize("name", ALL_MEDIA)
def test_references_are_resolvable(name):
    """Every reference is one line carrying a year and, unless it is vendor data
    or a work with no registered DOI, a doi.org URL. The 25 DOIs in the package
    were each resolved against Crossref when this style was introduced; this
    test only guards the shape, since resolving them needs the network."""
    import re
    doc = getattr(C, name).__doc__
    body = doc[doc.index('Ref'):]
    body = body[body.index('---') + 3:]
    lines = [l.strip() for l in body.split('\n') if l.strip() and not l.strip().endswith(':')]
    assert lines, f"{name}: Ref section is empty"
    for line in lines:
        assert re.search(r'\(\d{4}\)|accessed \d{4}-\d{2}-\d{2}', line), \
            f"{name}: reference has no year: {line[:60]}"
        assert 'https://doi.org/' in line or 'vendor data' in line \
            or 'not registered with Crossref' in line or 'Handbook of Optics' in line, \
            f"{name}: reference has no DOI: {line[:60]}"


@pytest.mark.parametrize("name", ALL_MEDIA)
def test_docstring_template(name):
    """The class docstrings are the source of the documentation's media catalog,
    so template drift is documentation drift. Mandatory sections must be present;
    the recurring copy-paste typos and the deleted Example/Usage sections must
    not come back."""
    doc = getattr(C, name).__doc__
    assert doc, f"{name} has no docstring"
    for section in ("Sellmeier equation", "Validity range", "Ref"):
        assert section in doc, f"{name}: missing '{section}' section"
    for banned in ("Dielectic", "Tranparency", "Example\n", "Usage\n"):
        assert banned not in doc, f"{name}: contains {banned.strip()!r}"
DISPERSION_METHODS = ['n', 'dn_wl', 'd2n_wl', 'd3n_wl', 'GD', 'GV', 'ng', 'GVD',
                      'TOD', 'woa_theta', 'woa_phi', 'dndT', 'dndT2']


@pytest.mark.parametrize("name", ALL_MEDIA)
def test_injected_angle_matches_hand_written_wrapper(name):
    """Every method, not just n: a wrapper that drops an argument is no longer a
    TypeError now that _func fills the angle in, because a short call is exactly
    what _func expects. The ten LBO yz/zx d3n_wl wrappers dropped T_degC and so
    silently pushed 0.5*pi into the temperature slot. Comparing each wrapper
    against the base is what catches that, and nothing else does."""
    from ndispers._baseclass import Medium
    x = getattr(C, name)()
    n_sym = len(x.symbols)
    for args in ([(0.4, 0.0, 20), (0.532, 0.3, 25), (1.064, 1.2, 150)] if n_sym == 4
                 else [(0.4, 0.0), (0.532, 0.3), (1.064, 1.2)]):
        for meth in DISPERSION_METHODS:
            for pol in ('o', 'e'):
                try:
                    by_wrapper = float(getattr(x, meth)(*args, pol=pol))
                except (ValueError, TypeError):   # pol not defined, or complex outside the range
                    continue
                by_base = float(getattr(Medium, meth)(x, *x._full_args(args), pol=pol))
                # outside a crystal's range both are nan (e.g. AgGaSe2 at 0.4 um)
                assert by_base == by_wrapper or (np.isnan(by_base) and np.isnan(by_wrapper)), (name, meth, args, pol)


@pytest.mark.parametrize("name", ALL_MEDIA)
def test_every_dispersion_method_is_wired(name):
    """Each medium used to re-implement all ~12 public methods as a wrapper that
    only injects its fixed angle, and any one of them could be forgotten: dndT2
    was missing everywhere, d3n_wl in the eight LBO yz/zx classes, and AlphaBBO
    and Calcite demanded a T_degC they have no symbol for. _func fills the angle
    in now, so a missing wrapper is no longer a missing method."""
    x = getattr(C, name)()
    args = (0.532, 0.3, 25) if len(x.symbols) == 4 else (0.532, 0.3)
    for meth in DISPERSION_METHODS:
        for pol in ('o', 'e'):
            try:
                getattr(x, meth)(*args, pol=pol)
            except ValueError:
                pass          # medium legitimately refuses this polarization


# ---------------------------------------------------------------------------
# A uniaxial crystal's docstring says whether it is negative or positive; the
# Sellmeier sets must agree. This is the check that catches swapped o/e
# columns - Zelmon 1997's Table 2 for MgO:LiNbO3 is printed with its n_e and
# n_o headings interchanged, which is exactly the transcription trap.
@pytest.mark.parametrize("name", [n for n in ALL_MEDIA if n != "SLN"])
def test_optic_sign_matches_docstring(name):
    x = getattr(C, name)()
    doc = type(x).__doc__
    m = re.search(r"(Negative|Positive) uniaxial", doc)
    if not m:
        pytest.skip("biaxial")
    no = x.n(1.0, 0.0, 25, pol="o")
    ne = x.n(1.0, 0.5 * np.pi, 25, pol="e")
    if m.group(1) == "Negative":
        assert ne < no, f"{name}: docstring says negative uniaxial but n_e = {ne:.4f} > n_o = {no:.4f}"
    else:
        assert ne > no, f"{name}: docstring says positive uniaxial but n_e = {ne:.4f} < n_o = {no:.4f}"


def test_pmAngles_skips_unsupported_polarizations():
    """SLN has no ordinary-ray Sellmeier equation. pmAngles_sfg used to raise
    ValueError out of dk_sfg for every combination containing 'o', which took
    the phase-matching app down when SLN was selected; those combinations now
    report no solution."""
    pm = C.SLN().pmAngles_sfg(1.064, 1.064, 25, deg=True)
    assert pm['wl3'] == pytest.approx(0.532)
    for key in ('ooe', 'eeo', 'oee', 'eoe', 'eoo', 'oeo'):
        assert pm[key]['theta'] == []
    # a medium that does support both rays still finds its solutions
    assert C.BetaBBO_Eimerl1987().pmAngles_sfg(1.064, 1.064, 25, deg=True)['ooe']['theta']


@pytest.mark.parametrize("name", ALL_MEDIA)
def test_docstring_states_a_wavelength_range(name):
    """The explorer sweeps each medium over the range its docstring states, so
    every medium needs one that parses and is ordered."""
    import re
    doc = getattr(C, name).__doc__ or ""
    for pattern in (r"Validity range\s*\n\s*-+\s*\n(.{0,600}?)(?:\n\s*\n|$)",
                    r"Transparency range\s*:([^\n]*)"):
        m = re.search(pattern, doc, re.S)
        if m and re.search(r"(\d+\.?\d*)\s*(?:to|-|–)\s*(\d+\.?\d*)", m.group(1)):
            got = re.search(r"(\d+\.?\d*)\s*(?:to|-|–)\s*(\d+\.?\d*)", m.group(1))
            lo, hi = float(got.group(1)), float(got.group(2))
            assert 0 < lo < hi, f"{name}: range {lo}-{hi} um is not ordered"
            return
    pytest.fail(f"{name}: no parsable validity or transparency range in the docstring")


def test_noncritical_phase_matching_is_found():
    """A noncritical solution sits exactly on an end of the scanned 0-90 degree
    range, where dk never changes sign, so the sign-change scan alone reported
    no solution for every NCPM cut. At the exact NCPM point the endpoint is now
    returned; a tenth of a degree or a tenth of a nanometre away there really is
    no solution, and none is reported.

    The three NCPM points below agree with the sources: LBO's 90-degree
    temperature for 1064 nm SHG is 147.0 C here against Kato, Grechin & Umemura
    2018 (149 C calculated, 148-151 C observed); Li2B4O7's shortest SHG is
    487.7 -> 243.9 nm against Komatsu et al. 1997 Table I (487.6 -> 243.8 nm at
    theta = 90 deg); CLBO's is 473.4 -> 236.7 nm against Umemura & Kato 1997
    (236.8 nm).
    """
    from ndispers.helper import brentq

    lbo = C.LBO_KK2018_xy()                    # xy plane: the NCPM cut is phi = 0
    T_ncpm = brentq(lambda T: float(lbo.dk_sfg(1.064, 1.064, 0.0, T, 'o', 'o', 'e')), 100, 300)
    assert T_ncpm == pytest.approx(147.0, abs=1.0)
    assert lbo.pmAngles_sfg(1.064, 1.064, T_ncpm, deg=True)['ooe']['phi'] == pytest.approx([0.0], abs=0.01)
    assert lbo.pmAngles_sfg(1.064, 1.064, T_ncpm + 0.1, deg=True)['ooe']['phi'] == []

    for x, expect_um, lo, hi in ((C.LB4(), 0.4877, 0.45, 0.55), (C.CLBO(), 0.4734, 0.45, 0.52)):
        wl = brentq(lambda w: float(x.dk_sfg(w, w, 0.5 * np.pi, 20, 'o', 'o', 'e')), lo, hi)
        assert wl == pytest.approx(expect_um, abs=5e-4)
        assert x.pmAngles_sfg(wl, wl, 20, deg=True)['ooe']['theta'] == pytest.approx([90.0], abs=0.01)
        # shorter than the NCPM wavelength there is no solution at all; longer,
        # the angle backs off the endpoint in the ordinary way
        assert x.pmAngles_sfg(wl - 1e-4, wl - 1e-4, 20, deg=True)['ooe']['theta'] == []
        assert x.pmAngles_sfg(wl + 1e-4, wl + 1e-4, 20, deg=True)['ooe']['theta'][0] == pytest.approx(88.7, abs=0.3)


def test_pmAngles_does_not_invent_endpoint_solutions():
    """The endpoint rule must not fire where dk is merely small: beta-BBO at
    1064 nm has one type-I solution at 22.88 deg and nothing at 90 deg."""
    pm = C.BetaBBO_Eimerl1987().pmAngles_sfg(1.064, 1.064, 25, deg=True)
    assert pm['ooe']['theta'] == pytest.approx([22.884], abs=0.01)
    assert pm['eeo']['theta'] == []


def test_brentq_replaces_scipy():
    """scipy was dropped in 0.15 (19 MB of the browser apps' 48 MB download,
    for one root finder). helper.brentq must hit a root to the requested
    tolerance, honour exact endpoints, and refuse an unbracketed interval."""
    from math import cos, pi
    from ndispers.helper import brentq
    assert abs(brentq(cos, 1.0, 2.0) - pi / 2) < 1e-12
    assert abs(brentq(lambda x: x ** 3 - 2 * x - 5, 2.0, 3.0) - 2.0945514815423265) < 1e-12
    assert brentq(lambda x: x, 0.0, 1.0) == 0.0                      # root on an endpoint
    with pytest.raises(ValueError):
        brentq(lambda x: x * x + 1, -1.0, 1.0)
    # the number the README and validation page quote must not move
    assert abs(C.BetaBBO_Eimerl1987().pmAngles_sfg(1.064, 1.064, 25, deg=True)
               ['ooe']['theta'][0] - 22.884169498625802) < 1e-9


# --- pre-generated functions (ndispers/_compiled, tools/compile_media.py) ----

def _all_media():
    import ndispers.media.glasses as G
    for mod in (C, G):
        for name in sorted(dir(mod)):
            obj = getattr(mod, name)
            if isinstance(obj, type) and not name.startswith("_"):
                yield obj


def test_compiled_modules_are_current():
    """Every medium ships a generated module whose hash matches the sources
    it was generated from. A coefficient edited without re-running
    tools/compile_media.py would otherwise be silently ignored at run time."""
    import importlib.util, os
    spec = importlib.util.spec_from_file_location(
        "compile_media", os.path.join(os.path.dirname(__file__), "..", "..", "tools", "compile_media.py"))
    gen = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gen)
    from ndispers._baseclass import _compiled_module
    stale = []
    for cls in _all_media():
        mod = _compiled_module(cls)
        if mod is None or mod.SOURCE_HASH != gen.source_hash(cls):
            stale.append(cls.__name__)
    assert not stale, f"run `uv run python tools/compile_media.py`; stale: {stale}"


@pytest.mark.parametrize("cls", [C.BetaBBO_Eimerl1987, C.LBO_KK2018_xy, C.SLN,
                                 C.BiBO_Miyata2009_yz, C.KTP_zx])
def test_compiled_functions_agree_with_sympy(cls):
    """The shipped functions reproduce sympy's lambdify (cse merely reorders
    the arithmetic) for every expression and polarization of a medium."""
    import warnings
    from sympy.utilities import lambdify
    from ndispers._baseclass import _compiled_module
    x = cls()
    mod = _compiled_module(cls)
    syms = [s._sympy_() for s in x.symbols]
    wl = np.array([0.4, 0.8, 1.064, 1.55, 2.5])[:, None, None]
    th = np.array([0.0, 0.4, 1.0, 1.5])[None, :, None]
    T = np.array([20., 120.])[None, None, :]
    args = (wl, th, 0.3, T) if len(syms) == 4 else (wl, T)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for (ename, pol), f in mod.FUNCS.items():
            expr = getattr(x, ename)(pol) if ename.endswith("_expr") \
                else getattr(x, ename + "_expr")()
            a = np.asarray(f(*args), float)
            b = np.asarray(lambdify(syms, expr, "numpy")(*args), float)
            m = np.isfinite(b)
            assert np.array_equal(np.isfinite(a), m), (ename, pol)
            assert np.allclose(a[m], b[m], rtol=1e-10, atol=1e-10 * np.max(np.abs(b[m]), initial=0)), (ename, pol)


def test_package_runs_without_sympy():
    """sympy is optional since 0.16: indices, dispersion, phase matching and
    d_eff come from the generated modules. Run in a subprocess with sympy
    blocked so this test cannot be fooled by the test process having it."""
    import subprocess, sys
    code = r"""
import sys
sys.modules['sympy'] = None          # makes `import sympy` raise ImportError
sys.modules['mpmath'] = None
import numpy as np
import ndispers as nd
C, G = nd.media.crystals, nd.media.glasses
b = C.BetaBBO_Eimerl1987()
assert abs(b.n(0.532, 0, 25, pol='o') - 1.674884049110459) < 1e-12
assert abs(b.pmAngles_sfg(1.064, 1.064, 25, deg=True)['ooe']['theta'][0] - 22.884169498625802) < 1e-9
b.GVD(0.8, 0.3, 25, pol='e'); b.TOD(0.8, 0.3, 25, pol='e'); b.woa_theta(0.8, 0.4, 25, pol='e')
b.deff_sfg(1.064, 1.064, np.radians(22.88), np.radians(90), 25, 'o', 'o', 'e')
l = C.LBO_KK2018_xy(); l.pmAngles_sfg(1.064, 1.064, 25); l.d_sfg('d32', 1.064, 1.064, 25)
C.SLN().qpm_period_sfg(1.064, 1.064, np.pi/2, 25, 'e', 'e', 'e')
G.FusedSilica().GVD(0.8, 20)
assert b.n_latex('o').startswith('- 1.66')
for cls in (C.BetaBBO_Eimerl1987, C.LBO_KK2018_xy, C.SLN, C.KTP_zx):
    cls().n(1.064, 0.5, 25)
# a user subclass has no generated module and needs sympy: say so clearly
class MyBBO(C.BetaBBO_Eimerl1987):
    __slots__ = []
try:
    MyBBO().n(0.5, 0, 25)
except ImportError as e:
    assert 'ndispers[sym]' in str(e)
else:
    raise AssertionError('expected ImportError without sympy')
print('OK')
"""
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert r.returncode == 0 and r.stdout.strip() == "OK", r.stderr
