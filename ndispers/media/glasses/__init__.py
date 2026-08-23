"""
ndispers.media.glasses

Optically isotropic media: glasses, and crystals of the cubic system.

The split between this package and ndispers.media.crystals is by optical
isotropy, not by whether a medium is crystalline. An isotropic medium has a
single refractive index with no direction or polarization dependence, so its
methods take (wl_um, T_degC) and no angle or pol argument. CaF2 is a crystal
(cubic, m3-bar-m) and lives here for that reason; alpha-BBO and calcite are
centrosymmetric too but are uniaxial, so they are birefringent and belong with
the crystals.
"""

from ._fusedsilica import FusedSilica
from ._caf2 import CaF2
from ._yag import YAG
from ._nbk7 import NBK7
from ._lif import LiF
from ._baf2 import BaF2
from ._znse import ZnSe
from ._zns import ZnS
from ._si import Si
from ._ge import Ge
from ._diamond import Diamond
from ._sf10 import SF10
from ._sf11 import SF11
from ._sf57 import SF57
