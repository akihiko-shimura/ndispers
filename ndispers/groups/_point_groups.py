"""
Point groups: the d tensor pattern of each, with Kleinman symmetry applied.

The matrices below are the standard forms (Boyd, Nonlinear Optics, Table
1.5.3; Roberts, IEEE JQE 28, 2057 (1992)) in the contracted notation
d_il, i = x,y,z and l = xx,yy,zz,yz,zx,xy. d_eff is not written out per
group: NonlinearGroup.deff_sfg contracts the tensor with the field vectors,
and tests/test_nonlinear.py checks the result against the published
closed forms.
"""
from ._base import UniaxialGroup


class Uniax_3m(UniaxialGroup):
    """
    Point group 3m (C3v): beta-BBO, LiNbO3, LiTaO3. Mirror plane
    perpendicular to x.

        d = [[   0,    0,   0,   0,  d15, -d22 ],
             [ -d22,  d22,  0,  d15,  0,    0  ],
             [  d31,  d31, d33,  0,   0,    0  ]],   d15 = d31 (Kleinman)

    Independent: d22, d31, d33. d33 does not enter d_eff of any interaction
    with an ordinary wave; it matters for eee (quasi-phase-matching).
    """
    __slots__ = []
    _d_entries = {
        "d22": [("yyy", +1), ("yxx", -1), ("xxy", -1)],
        "d31": [("zxx", +1), ("zyy", +1), ("xzx", +1), ("yzy", +1)],
        "d33": [("zzz", +1)],
    }


class Uniax_32(UniaxialGroup):
    """
    Point group 32 (D3): KBBF, RBBF, alpha-quartz.

        d = [[ d11, -d11,  0,  d14,   0,    0  ],
             [  0,    0,   0,   0,  -d14, -d11 ],
             [  0,    0,   0,   0,    0,    0  ]]

    Kleinman symmetry gives d14 = d25 = -d14, so d14 = 0 and d11 is the
    only independent component.
    """
    __slots__ = []
    _d_entries = {
        "d11": [("xxx", +1), ("xyy", -1), ("yxy", -1)],
    }


class Uniax_42m(UniaxialGroup):
    """
    Point group -42m (D2d): KDP, CLBO.

        d = [[ 0, 0, 0, d14,  0,   0  ],
             [ 0, 0, 0,  0,  d14,  0  ],
             [ 0, 0, 0,  0,   0,  d36 ]],   d14 = d36 (Kleinman)
    """
    __slots__ = []
    _d_entries = {
        "d36": [("zxy", +1), ("xyz", +1), ("yzx", +1)],
    }


class Uniax_4mm(UniaxialGroup):
    """
    Point group 4mm (C4v): Li2B4O7.

        d = [[  0,    0,   0,   0,  d15,  0 ],
             [  0,    0,   0,  d15,  0,   0 ],
             [ d31,  d31, d33,  0,   0,   0 ]],   d15 = d31 (Kleinman)

    d33 enters only eee. Type-II d_eff vanishes identically.
    """
    __slots__ = []
    _d_entries = {
        "d31": [("zxx", +1), ("zyy", +1), ("xzx", +1), ("yzy", +1)],
        "d33": [("zzz", +1)],
    }
