"""
Machine-readable catalog of every medium, for programmatic (and AI-agent)
discovery: filter by wavelength range, point group or optical class without
reading the docs.
"""
import re

from . import media


def _point_group(cls):
    # non-centrosymmetric crystals inherit a class from ndispers.groups
    # named for the point group (Uniax_3m, Biax_mm2, ...)
    for base in cls.__mro__:
        mod = getattr(base, '__module__', '')
        if mod.startswith('ndispers.groups'):
            return base.__name__.split('_', 1)[1]
    return None


def catalog():
    """
    One dict per medium class, for programmatic discovery.

    Keys: name (class name, states the Sellmeier source), kind ('crystal' or
    'glass'), description (docstring first line), point_group (None for media
    with no second-order nonlinearity), plane ('xy'/'yz'/'zx' for a principal-
    plane class of a biaxial crystal, else None), wl_range ((min, max) um of
    the Sellmeier validity, None if the source states none), d_components
    (nonlinear tensor component names, e.g. ['d22', 'd31', 'd33']).

    Examples
    --------
    Media usable at 3 um with a nonlinear tensor:

    >>> import ndispers as nd
    >>> [e['name'] for e in nd.catalog()
    ...  if e['wl_range'] and e['wl_range'][0] < 3 < e['wl_range'][1]
    ...  and e['d_components']]  # doctest: +NORMALIZE_WHITESPACE
    ['AGS_Kato1996', 'AGS_Takaoka1999', 'AGSe_Kato2021', 'BetaBBO_KK2010',
     'BetaBBO_Tamosauskas2018', 'BiBO_Miyata2009_xy', 'BiBO_Miyata2009_yz',
     'BiBO_Miyata2009_zx', 'KTP_xy', 'KTP_yz', 'KTP_zx', 'LiIO3',
     'MgOLN_Zelmon1997', 'RBBF', 'SLN', 'SLT', 'ZGP_Das2003', 'ZGP_Zelmon2001']
    """
    out = []
    for kind, pkg in (('crystal', media.crystals), ('glass', media.glasses)):
        for name in dir(pkg):
            cls = getattr(pkg, name)
            if not (isinstance(cls, type) and name[0].isupper()):
                continue
            doc = (cls.__doc__ or '').strip().splitlines()
            inst = cls()
            out.append({
                'name': name,
                'kind': kind,
                'description': doc[0].strip() if doc else '',
                'point_group': _point_group(cls),
                'plane': None if inst.plane == 'arb' else inst.plane,
                'wl_range': cls._wl_range,
                'd_components': sorted(getattr(cls, '_d_ref', {})),
                'dois': sorted(set(re.findall(
                    r'doi\.org/(10\.\S+?)(?=[\s)>,]|$)', cls.__doc__ or ''))),
            })
    return out
