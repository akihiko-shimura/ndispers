---
title: 'ndispers: refractive-index dispersion and phase matching of nonlinear-optics media in Python'
tags:
  - Python
  - nonlinear optics
  - ultrafast optics
  - refractive index
  - dispersion
  - phase matching
authors:
  - name: Akihiko Shimura
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
  - name: Independent researcher, Japan
    index: 1
date: 30 August 2026
bibliography: paper.bib
---

# Summary

Frequency conversion of laser light — second-harmonic generation, sum- and
difference-frequency mixing, optical parametric amplification and
oscillation — depends on the refractive-index dispersion of the nonlinear
crystal: the propagation angles, temperatures and periods at which the
interacting waves stay in phase are all set by n(λ, θ, T). `ndispers` is a
Python library that computes these quantities from Sellmeier equations and
thermo-optic coefficients transcribed from the literature. It provides
refractive indices, group-velocity quantities and their temperature
derivatives, walk-off, phase-matching angles, acceptance widths,
quasi-phase-matching periods, parametric tuning curves, and effective
nonlinear coefficients with Miller wavelength scaling, for 60 media spanning
nonlinear crystals, linear birefringent crystals and optical glasses. All
methods accept and return numpy arrays [@harris2020], so the quantities can
be evaluated over wavelength, angle and temperature grids inside a larger
simulation; the evaluated functions are generated ahead of time from sympy
expressions [@meurer2017], leaving numpy as the only runtime dependency.

# Statement of need

The standard tools for these calculations are interactive applications —
most notably SNLO [@smith2018], with web and mobile calculators covering
parts of the same ground [@refractiveindexinfo; @lightcon]. An application
answers the query typed into it; it cannot be called from a program, so
simulation work (pulse propagation, OPO design, thermal analysis) either
re-implements the dispersion layer or copies numbers by hand. `ndispers`
provides that layer as a library.

Two design choices address reproducibility. First, every coefficient is
traceable: each medium class documents the paper, the validity range, and
separate sources for its Sellmeier equation, thermo-optic coefficients and
d tensor, so naming the class and package version specifies a calculation
completely. Where the literature offers several parameterisations of one
crystal, each is a separate class and the disagreement between sources is
visible rather than hidden behind a default. Second, correctness is
maintained by machinery: coefficients are cross-checked by independent
re-extraction from the sources, computed values are pinned by regression
tests, and a published validation record compares the package against
fidelity and independent references.

The library is also designed to be operable by AI assistants: media are
enumerable programmatically, docstrings carry the metadata an agent needs,
out-of-range and misuse conditions raise instructive warnings and errors,
and a single-page API reference is published in the llms.txt convention.

`ndispers` is intended for researchers, engineers and students in nonlinear
and ultrafast optics who need dispersion and phase-matching quantities
inside their own code, with sources they can cite.

# Acknowledgements

The author thanks the authors of the cited dispersion measurements, on whose
published coefficients this package rests.

# References
