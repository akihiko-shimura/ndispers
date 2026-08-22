# Browser apps

**Live: <https://akihiko-shimura.github.io/ndispers/>**

Two [marimo](https://marimo.io) notebooks that run entirely in the browser —
no server, no install. ndispers is a pure-Python wheel and numpy, scipy and
sympy all ship with Pyodide, so the whole calculation happens client-side.

| app | what it does |
|---|---|
| `explorer.py` | plot any of n, n_g, GVD, TOD, dn/dT or walk-off against any of wavelength, temperature or angle, for every medium |
| `phasematching.py` | collinear SFG: phase-matching angles, acceptance bandwidths, walk-off, effective nonlinearity |

## Run locally

```
uv sync --group notebook
uv run marimo edit apps/explorer.py
```

## Build the static site

```
uv run marimo export html-wasm apps/explorer.py -o site/explorer --mode run
uv run marimo export html-wasm apps/phasematching.py -o site/phasematching --mode run
```

Each export is self-contained and can be served by any static host. Measured
on the explorer: about 3.7 s from page load to first result (Pyodide plus the
scientific stack is ~28 MB, cached by the browser afterwards), then ~60 ms to
redraw a 400-point dispersion curve and ~0.6 s when switching crystal, since
that rebuilds the symbolic derivatives.

## Notes on the physics

`phasematching.py` computes acceptance bandwidths as the full width of
sinc²(ΔkL/2) at a chosen threshold, solved with Brent's method rather than by
stepping to the first zero. Cross-checked against the step-and-`fullWidth`
implementation in `opt-workspace/notebooks/pmbandwidth_calc.py`, it agrees to
better than 0.1%.

The spectral acceptance is reported two ways, because they differ by exactly a
factor of two and only one matches published tables. Holding λ₂ at its stated
value and tuning λ₁ alone is what SNLO reports: β-BBO Type-I SHG of 800 nm
gives 1.00 nm·cm against the ~0.9 nm·cm in the literature. Tuning λ₁ and λ₂
together makes Δk move twice as fast — Δk depends on them only through 1/λ₁ and
1/λ₂ with equal weight — so that width is half, 0.50 nm·cm. It is the relevant
one when a single broadband beam supplies both photons. At 1.064 µm the two are
4.23 and 2.12 nm·cm.

Lengths default to 10 mm, so the tabulated ×L products read directly in
per-centimetre units. The units of the acceptance table are selectable —
wavelength in pm, nm or µm and angle in mdeg, mrad, µrad or deg — because the
values span decades across crystals and interactions; the defaults are pm and
mdeg.

The azimuthal cut φ is a control rather than a fixed value. A refractive index
in a uniaxial crystal does not depend on φ — which is why ndispers reports
`phi_rad` as `'arb'` — but d_eff does, through a sin(3φ) or cos(3φ) factor. At
φ = 0 the d22 term of a 3m crystal vanishes entirely: β-BBO Type-I SHG gives
0.018 pm/V there against 1.96 pm/V at the conventional φ = 30° cut.
