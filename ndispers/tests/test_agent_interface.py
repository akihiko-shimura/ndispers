"""
The agent-facing surface must not drift from the code.

llms.txt is the reference AI assistants act on; a stale name in it is worse
than no reference. These tests verify every identifier it mentions exists,
and run the doctests that agents will copy.
"""
import doctest
import re
from pathlib import Path

import pytest

import ndispers as nd
import ndispers._catalog
from ndispers._baseclass import Medium

LLMS = Path(__file__).resolve().parents[2] / 'docs' / 'llms.txt'

# backticked tokens in llms.txt that are not ndispers attributes
NOT_API = {
    'pip install ndispers', 'dir(nd.media.crystals)', 'dir(nd.media.glasses)',
    'help(obj)', 'obj?', 'obj.constants', 'obj.theta_rad', 'obj.phi_rad',
    'obj.dndT(...) != 0', 'nd.media.glasses', "pol='o'", "'o'", "'e'",
    "d_sfg(\"d31\", wl1, wl2, T)", 'qpm_period=', '_xy', '_yz', '_zx',
    '(wl_um, angle_rad, T_degC, pol)', '(wl_um, T_degC)', 'deg=True',
    'deff_sfg(wl1, wl2, theta_rad, phi_rad, T, pol1, pol2,\npol3)',
    'pmAngles_sfg(wl1, wl2,\nT_degC, deg=False)',
    # argument names / builtins / private-but-documented attributes
    'pol', '_d_note', 'n_expr', 'repr(obj)', 'wl_i_max=', 'qpm_period=',
    'pip install ndispers[eff]', 'details=True',
    'T_degC', 'GVM(wl1, wl2, angle_rad, T_degC, pol1, pol2)', 'dois',
}


def catalog_names():
    return {e['name'] for e in nd.catalog()}


@pytest.mark.skipif(not LLMS.exists(), reason='docs/ not present (installed package)')
def test_llms_txt_method_names_exist():
    bbo = nd.media.crystals.BetaBBO_Eimerl1987()
    # drop fenced code blocks; scan inline `...` tokens only
    text = re.sub(r'```.*?```', '', LLMS.read_text(), flags=re.S)
    for tok in re.findall(r'`([^`]+)`', text):
        if tok in NOT_API:
            continue
        name = re.match(r'\w+', tok.removeprefix('nd.').removeprefix('ndispers.'))
        if name is None:
            continue
        name = name.group()
        assert (hasattr(bbo, name) or hasattr(Medium, name)
                or name in catalog_names() or hasattr(nd, name)), \
            f'llms.txt mentions `{tok}` but {name!r} does not exist in the API'


@pytest.mark.skipif(not LLMS.exists(), reason='docs/ not present (installed package)')
def test_llms_txt_media_list_matches_package():
    """Every class in the package appears in llms.txt's Media section, after
    expanding the BetaBBO_{A, B}-style brace shorthand."""
    text = LLMS.read_text()
    sec = text.split('## Media')[1].split('## Gotchas')[0]
    listed = set()
    for stem, alts in re.findall(r'(\w+)_\{([^}]+)\}', sec):
        parts = [a.strip() for a in alts.split(',')]
        for a in parts:
            listed.add(f'{stem}_{a}')
        # nested shorthand LBO_{src}_{plane}: expand the plane suffix too
        m = re.search(re.escape('{' + ', '.join(parts) + '}') + r'_\{([^}]+)\}', sec)
        if m:
            listed |= {f'{stem}_{a}_{p.strip()}' for a in parts
                       for p in m.group(1).split(',')}
    listed |= set(re.findall(r'\b[A-Z]\w+\b', re.sub(r'\w+_\{[^}]+\}(_\{[^}]+\})?', '', sec)))
    missing = catalog_names() - listed
    assert not missing, f'media missing from llms.txt: {sorted(missing)}'


@pytest.mark.skipif(not LLMS.exists(), reason='docs/ not present (installed package)')
def test_llms_full_txt_is_current():
    """docs/llms-full.txt is generated; regenerate with
    `uv run python tools/gen_llms_full.py` after editing its sources."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        'gen_llms_full', LLMS.parents[1] / 'tools' / 'gen_llms_full.py')
    gen = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gen)
    assert (LLMS.parent / 'llms-full.txt').read_text() == gen.build()


def test_doctests():
    import warnings
    import ndispers.groups._base
    warnings.simplefilter('ignore', nd.ModelValidityWarning)
    for mod in (ndispers._catalog, nd._baseclass, ndispers.groups._base):
        r = doctest.testmod(mod, verbose=False)
        assert r.failed == 0, f'{mod.__name__}: {r.failed} doctest failures'


def test_temperature_warning_fires_once_for_T_independent_media():
    import warnings
    from ndispers._baseclass import Medium, TemperatureWarning
    Medium._T_checked.discard(nd.media.crystals.MgOLN_Zelmon1997)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        m = nd.media.crystals.MgOLN_Zelmon1997()
        m.n(1.064, 0.5, 25)
        m.n(1.064, 0.5, 100)
        assert sum(isinstance(x.message, TemperatureWarning) for x in w) == 1
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        nd.media.crystals.KTP_zx().n(1.064, 0.5, 25)   # has dn/dT: no warning
        assert not any(isinstance(x.message, TemperatureWarning) for x in w)


def test_tuning_dfg_default_stops_at_validity_edge():
    import numpy as np
    sln = nd.media.crystals.SLN()
    roots = sln.tuning_dfg(1.064, np.pi/2, 100, 'e', 'e', 'e', qpm_period=30.0)
    assert all(wl_i <= sln._wl_range[1] for _, wl_i in roots)
    more = sln.tuning_dfg(1.064, np.pi/2, 100, 'e', 'e', 'e', qpm_period=30.0,
                          wl_i_max=20)
    assert len(more) > len(roots)      # explicit extrapolation still possible


def test_gvm_is_gd_difference():
    b = nd.media.crystals.BetaBBO_Eimerl1987()
    expect = b.GD(0.532, 0.4, 25, pol='e') - b.GD(1.064, 0.4, 25, pol='o')
    assert b.GVM(0.532, 1.064, 0.4, 25, pol1='e', pol2='o') == expect


def test_acceptance_sfg_matches_direct_sinc2_scan():
    """BBO 800 nm type-I SHG, L = 1 mm: values pinned against an independent
    sinc2 scan (and the literature ballpark 0.5 nm cm, 0.31 mrad cm)."""
    from math import radians
    b = nd.media.crystals.BetaBBO_Tamosauskas2018()
    th = radians(b.pmAngles_sfg(0.8, 0.8, 25, deg=True)['ooe']['theta'][0])
    import numpy as np
    assert np.isclose(b.acceptance_sfg(0.8, 0.8, th, 25, 'o', 'o', 'e', 1.0, 'wl'),
                      4.911e-3, rtol=1e-3)
    assert np.isclose(b.acceptance_sfg(0.8, 0.8, th, 25, 'o', 'o', 'e', 1.0, 'theta'),
                      3.143e-3, rtol=1e-3)
    assert np.isclose(b.acceptance_sfg(0.8, 0.8, th, 25, 'o', 'o', 'e', 1.0, 'T'),
                      174.0, rtol=1e-2)
    with pytest.raises(ValueError, match='not\\s+phase-matched'):
        b.acceptance_sfg(0.8, 0.8, 0.1, 25, 'o', 'o', 'e', 1.0, 'wl')
    with pytest.raises(ValueError, match='wl1 == wl2'):
        b.acceptance_sfg(1.0, 0.5, th, 25, 'o', 'o', 'e', 1.0, 'wl')
