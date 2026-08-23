"""
ndispers.groups - point-group classes that add the second-order
nonlinearity (Miller-rule scaled tensor components and d_eff) to Medium.
"""
from ._base import NonlinearGroup, UniaxialGroup, BiaxialGroup
from ._point_groups import Uniax_3m, Uniax_32, Uniax_42m, Uniax_4mm, Biax_mm2, Biax_2

# names used before v0.10
Uniax_neg_3m = Uniax_3m
Uniax_neg_32 = Uniax_32
