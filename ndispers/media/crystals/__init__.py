"""
ndispers.media.crystals/__init__.py
anisotropic crystals
"""

from ._AGS_Kato1996 import AGS as AGS_Kato1996
from ._AGS_Takaoka1999 import AGS as AGS_Takaoka1999
from ._AGSe_Kato2021 import AGSe as AGSe_Kato2021
from ._alphaBBO import AlphaBBO
from ._betaBBO_Eimerl1987 import BetaBBO as BetaBBO_Eimerl1987
from ._betaBBO_Ghosh1995 import BetaBBO as BetaBBO_Ghosh1995
from ._betaBBO_KK2010 import BetaBBO as BetaBBO_KK2010
from ._betaBBO_Tamosauskas2018 import BetaBBO as BetaBBO_Tamosauskas2018
from ._BiBO_Miyata2009 import BiBO_xy as BiBO_Miyata2009_xy
from ._BiBO_Miyata2009 import BiBO_yz as BiBO_Miyata2009_yz
from ._BiBO_Miyata2009 import BiBO_zx as BiBO_Miyata2009_zx
from ._calcite import Calcite
from ._CLBO import CLBO
from ._DKDP import DKDP
from ._KBBF_Li2016 import KBBF
from ._KDP import KDP
from ._KTP import KTP_xy, KTP_yz, KTP_zx
from ._LB4 import LB4
from ._LBO_Castech import LBO_xy as LBO_Castech_xy
from ._LBO_Castech import LBO_yz as LBO_Castech_yz
from ._LBO_Castech import LBO_zx as LBO_Castech_zx
from ._LBO_Ghosh1995 import LBO_xy as LBO_Ghosh1995_xy
from ._LBO_Ghosh1995 import LBO_yz as LBO_Ghosh1995_yz
from ._LBO_Ghosh1995 import LBO_zx as LBO_Ghosh1995_zx
from ._LBO_KK1994 import LBO_xy as LBO_KK1994_xy
from ._LBO_KK1994 import LBO_yz as LBO_KK1994_yz
from ._LBO_KK1994 import LBO_zx as LBO_KK1994_zx
from ._LBO_KK2018 import LBO_xy as LBO_KK2018_xy
from ._LBO_KK2018 import LBO_yz as LBO_KK2018_yz
from ._LBO_KK2018 import LBO_zx as LBO_KK2018_zx
from ._LBO_Newlight import LBO_xy as LBO_Newlight_xy
from ._LBO_Newlight import LBO_yz as LBO_Newlight_yz
from ._LBO_Newlight import LBO_zx as LBO_Newlight_zx
from ._MgF2 import MgF2
from ._MgOLN_Zelmon1997 import MgOLN as MgOLN_Zelmon1997
from ._quartz import Quartz
from ._RBBF_Chen2009 import RBBF
from ._sapphire import Sapphire
from ._SLN_MgO_doped import SLN
from ._SLT_MgO_doped import SLT
from ._ZGP_Das2003 import ZGP as ZGP_Das2003
from ._ZGP_Zelmon2001 import ZGP as ZGP_Zelmon2001

# The classes above are exported under source-named aliases (several modules
# define a class literally named BetaBBO or LBO_xy). Stamp the public name
# onto each class so that repr, pickle-by-reference and autodoc all see the
# name users actually import. Python does not inherit class docstrings, so the
# principal-plane subclasses (KTP_xy, LBO_*_yz, ...) would otherwise show no
# material description in help() or the docs catalog - give them their base
# class's.
for _name, _cls in list(globals().items()):
    if not isinstance(_cls, type):
        continue
    if _cls.__name__ != _name:
        _cls.__name__ = _name
        _cls.__qualname__ = _name
        _cls.__module__ = __name__
    if _cls.__doc__ is None:
        _cls.__doc__ = next(
            (b.__doc__ for b in _cls.__mro__[1:] if b.__doc__), None)
del _name, _cls
