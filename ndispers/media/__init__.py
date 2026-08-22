"""
ndispers.media contains the collections of specific media.

They are split by optical isotropy, which is what decides a medium's call
signature:

- ndispers.media.crystals - birefringent (uniaxial or biaxial) media, whose
  methods take (wl_um, angle_rad, T_degC, pol).
- ndispers.media.glasses - optically isotropic media, glasses and cubic
  crystals alike, whose methods take (wl_um, T_degC).
"""
from . import crystals
from . import glasses
