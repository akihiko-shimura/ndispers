---
title: "Medium: dispersion and phase-matching methods"
description: "API reference of the ndispers Medium class: refractive index, group index, GVD, TOD, walk-off, phase mismatch, phase-matching angles, acceptance widths, QPM period."
---

# Base class

All dispersion and phase-matching methods live on `Medium` and are shared by
every crystal and glass. Individual media only supply the Sellmeier expression
and their fixed orientation.

::: ndispers._baseclass.Medium
    options:
      members: true
      show_root_heading: true
