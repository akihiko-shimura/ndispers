---
title: "Point groups: d tensor, Miller scaling and d_eff"
description: "How ndispers computes the effective nonlinear coefficient d_eff: point-group classes, d-tensor components and Miller-rule scaling to the working wavelengths."
---

# Point groups (second-order nonlinearity)

`ndispers.groups` adds Miller-rule scaling of the d tensor and the effective
coefficient `deff_sfg` to `Medium`. A crystal inherits the class of its point
group and supplies reference measurements; the geometry is generic.

::: ndispers.groups.NonlinearGroup
    options:
      members: ["miller_delta", "d_sfg", "deff_sfg"]
      show_root_heading: true

::: ndispers.groups.Uniax_3m
::: ndispers.groups.Uniax_32
::: ndispers.groups.Uniax_42m
::: ndispers.groups.Uniax_4mm
::: ndispers.groups.Biax_mm2
