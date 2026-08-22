Media catalog & API reference
=============================

Each crystal class below carries, in its docstring, the Sellmeier
equation it implements, its validity range and the literature the
coefficients come from. Crystals available in several
parameterisations are named after their source.

Crystals
--------

.. autoclass:: ndispers.media.crystals.AlphaBBO

.. autoclass:: ndispers.media.crystals.BetaBBO_Eimerl1987

.. autoclass:: ndispers.media.crystals.BetaBBO_Ghosh1995

.. autoclass:: ndispers.media.crystals.BetaBBO_KK2010

.. autoclass:: ndispers.media.crystals.BetaBBO_Tamosauskas2018

.. autoclass:: ndispers.media.crystals.CLBO

.. autoclass:: ndispers.media.crystals.Calcite

.. autoclass:: ndispers.media.crystals.KBBF

.. autoclass:: ndispers.media.crystals.KDP

.. autoclass:: ndispers.media.crystals.KTP_xy

.. autoclass:: ndispers.media.crystals.KTP_yz

.. autoclass:: ndispers.media.crystals.KTP_zx

.. autoclass:: ndispers.media.crystals.LB4

.. autoclass:: ndispers.media.crystals.LBO_Castech_xy

.. autoclass:: ndispers.media.crystals.LBO_Castech_yz

.. autoclass:: ndispers.media.crystals.LBO_Castech_zx

.. autoclass:: ndispers.media.crystals.LBO_Ghosh1995_xy

.. autoclass:: ndispers.media.crystals.LBO_Ghosh1995_yz

.. autoclass:: ndispers.media.crystals.LBO_Ghosh1995_zx

.. autoclass:: ndispers.media.crystals.LBO_KK1994_xy

.. autoclass:: ndispers.media.crystals.LBO_KK1994_yz

.. autoclass:: ndispers.media.crystals.LBO_KK1994_zx

.. autoclass:: ndispers.media.crystals.LBO_KK2018_xy

.. autoclass:: ndispers.media.crystals.LBO_KK2018_yz

.. autoclass:: ndispers.media.crystals.LBO_KK2018_zx

.. autoclass:: ndispers.media.crystals.LBO_Newlight_xy

.. autoclass:: ndispers.media.crystals.LBO_Newlight_yz

.. autoclass:: ndispers.media.crystals.LBO_Newlight_zx

.. autoclass:: ndispers.media.crystals.RBBF

.. autoclass:: ndispers.media.crystals.SLN

.. autoclass:: ndispers.media.crystals.SLT

Glasses
-------

.. autoclass:: ndispers.media.glasses.CaF2

.. autoclass:: ndispers.media.glasses.FusedSilica

Base class
----------

All dispersion methods live on the base class and are shared by
every medium.

.. autoclass:: ndispers._baseclass.Medium
   :members:
