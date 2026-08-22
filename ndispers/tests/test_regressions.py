"""Regression tests for bugs fixed in _baseclass.py and the KTP / SLN modules."""
import pickle

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
    from scipy.optimize import brentq
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
}


def _spread(names, meth, wl, T, pol):
    v = [float(getattr(getattr(C, nm)(), meth)(wl, 0.3, T, pol=pol)) for nm in names]
    return max(v) - min(v)


@pytest.mark.parametrize("group", sorted(SAME_CRYSTAL))
@pytest.mark.parametrize("wl,T", [(0.4, 20), (0.532, 20), (1.064, 20), (1.064, 100), (1.55, 20)])
@pytest.mark.parametrize("pol", ['o', 'e'])
def test_sources_agree_on_index(group, wl, T, pol):
    """Observed worst case is 1.6e-3; anything past 5e-3 is a transcription error,
    not source scatter (the KTP sign bug was 0.17-0.61)."""
    assert _spread(SAME_CRYSTAL[group], 'n', wl, T, pol) < 5e-3


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
    '4̄2m': ('tetragonal', False), '4mm': ('tetragonal', False),
    'mm2': ('orthorhombic', False), 'm3̄m': ('cubic', True),
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
                except ValueError:
                    continue
                by_base = float(getattr(Medium, meth)(x, *x._full_args(args), pol=pol))
                assert by_base == by_wrapper, (name, meth, args, pol)


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
