# /// script
# requires-python = ">=3.11"
# dependencies = ["ndispers==0.9.1", "numpy", "scipy", "matplotlib"]
# ///
"""Phase-matching calculator for sum-frequency generation.

Phase-matching angles, acceptance bandwidths (wavelength, angle, temperature),
walk-off and effective nonlinearity. Runs in the browser via marimo's WASM
export; see apps/README.md.
"""
import marimo

app = marimo.App(width="medium", app_title="ndispers phase matching")


@app.cell(hide_code=True)
def _():
    import marimo as mo
    return (mo,)


@app.cell(hide_code=True)
def _():
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.optimize import brentq
    import ndispers as nd
    import ndispers.media.crystals as _C

    # only media that can phase-match: one of the two angles must be free
    XTALS = {}
    for _name in sorted(n for n in dir(_C) if not n.startswith("_")):
        _cls = getattr(_C, _name)
        if not isinstance(_cls, type):
            continue
        _x = _cls()
        if "var" in (_x.theta_rad, _x.phi_rad):
            XTALS[_name] = _cls
    return XTALS, brentq, nd, np, plt


@app.cell(hide_code=True)
def _(XTALS, mo, nd):
    mo.md(
        f"""
    # Phase-matching calculator

    Collinear sum-frequency generation in {len(XTALS)} birefringent crystals —
    phase-matching angles, acceptance bandwidths, walk-off and effective
    nonlinearity, computed in your browser. Powered by
    [ndispers {nd.__version__}](https://github.com/akihiko-shimura/ndispers).
    """
    )
    return


@app.cell(hide_code=True)
def _(XTALS, mo):
    xtal = mo.ui.dropdown(options=list(XTALS), value="BetaBBO_Eimerl1987", label="Crystal")
    shg = mo.ui.checkbox(value=True, label="Second harmonic (λ₁ = λ₂)")
    wl1 = mo.ui.number(0.3, 4.0, value=1.064, step=0.001, label="λ₁ (µm)")
    wl2 = mo.ui.number(0.3, 4.0, value=1.064, step=0.001, label="λ₂ (µm)")
    T = mo.ui.number(-50, 300, value=25, step=1, label="Temperature (°C)")
    L = mo.ui.number(0.1, 50.0, value=1.0, step=0.1, label="Crystal length (mm)")
    return L, T, shg, wl1, wl2, xtal


@app.cell(hide_code=True)
def _(L, T, mo, shg, wl1, wl2, xtal):
    mo.hstack(
        [
            mo.vstack([xtal, shg]),
            mo.vstack([wl1, wl2 if not shg.value else mo.md("<sub>λ₂ = λ₁</sub>")]),
            mo.vstack([T, L]),
        ],
        justify="start",
        gap=2,
    )
    return


@app.cell(hide_code=True)
def _(XTALS, T, shg, wl1, wl2, xtal):
    x = XTALS[xtal.value]()
    WL1 = wl1.value
    WL2 = WL1 if shg.value else wl2.value
    WL3 = 1.0 / (1.0 / WL1 + 1.0 / WL2)
    T_C = T.value
    # which angle this crystal varies; the other one is fixed internally
    ANGLE_KEY = "theta" if x.theta_rad == "var" else "phi"
    return ANGLE_KEY, T_C, WL1, WL2, WL3, x


@app.cell(hide_code=True)
def _(ANGLE_KEY, T_C, WL1, WL2, WL3, mo, x):
    # Type I: both inputs share a polarization. Type II: they differ.
    TYPE_OF = {
        "ooe": "I (negative)", "eeo": "I (positive)",
        "oee": "II (negative)", "eoe": "II (negative)",
        "eoo": "II (positive)", "oeo": "II (positive)",
    }
    _pm = x.pmAngles_sfg(WL1, WL2, T_C, deg=True)
    SOLUTIONS = {k: _pm[k][ANGLE_KEY][0] for k in TYPE_OF
                 if isinstance(_pm.get(k), dict) and _pm[k].get(ANGLE_KEY)}

    _rows = "\n".join(
        f"| `{k}` | Type {TYPE_OF[k]} | "
        + (f"**{SOLUTIONS[k]:.3f}°**" if k in SOLUTIONS else "— none")
        + " |"
        for k in TYPE_OF
    )
    mo.md(
        f"""
    ## Phase-matching angles

    λ₁ = {WL1:.4f} µm, λ₂ = {WL2:.4f} µm → **λ₃ = {WL3:.4f} µm**
    &nbsp;·&nbsp; the varying angle of this crystal is **{ANGLE_KEY}**

    | polarizations (1,2,3) | type | angle |
    |---|---|---|
    {_rows}
    """
    )
    return (SOLUTIONS,)


@app.cell(hide_code=True)
def _(SOLUTIONS, mo):
    pmtype = mo.ui.dropdown(
        options=list(SOLUTIONS) or ["ooe"],
        value=(list(SOLUTIONS) or ["ooe"])[0],
        label="Interaction to analyse",
    )
    pmtype if SOLUTIONS else mo.md("*No phase-matching solution — change λ, T or crystal.*")
    return (pmtype,)


@app.cell(hide_code=True)
def _(SOLUTIONS, np, pmtype):
    HAS_PM = bool(SOLUTIONS) and pmtype.value in SOLUTIONS
    ANGLE_PM = np.radians(SOLUTIONS[pmtype.value]) if HAS_PM else np.nan
    POLS = tuple(pmtype.value) if HAS_PM else ("o", "o", "e")
    return ANGLE_PM, HAS_PM, POLS


@app.cell(hide_code=True)
def _(ANGLE_PM, HAS_PM, L, POLS, T_C, WL1, WL2, mo, np, plt, x):
    if HAS_PM:
        _w = 3.0 * np.pi / 180 * 0.05
        _a = np.linspace(ANGLE_PM - _w, ANGLE_PM + _w, 400)
        _dk = x.dk_sfg(WL1, WL2, _a, T_C, *POLS)
        _sinc = x.pmFactor_sfg(WL1, WL2, _a, T_C, *POLS, L.value)

        _fig, (_a1, _a2) = plt.subplots(2, 1, sharex=True, figsize=(7, 4.4))
        _a1.plot(np.degrees(_a), _dk)
        _a1.axhline(0, color="gray", lw=0.8)
        _a1.axvline(np.degrees(ANGLE_PM), color="r", lw=0.8, ls="--")
        _a1.set_ylabel(r"$\Delta k$ (rad/µm)")
        _a2.plot(np.degrees(_a), _sinc)
        _a2.axvline(np.degrees(ANGLE_PM), color="r", lw=0.8, ls="--")
        _a2.set_xlabel("angle (deg)")
        _a2.set_ylabel(r"sinc$^2(\Delta k L/2)$")
        for _ax in (_a1, _a2):
            _ax.grid(alpha=0.3)
        plt.tight_layout()
        _out = mo.vstack([mo.md("## Phase mismatch and conversion"), _fig])
    else:
        _out = mo.md("")
    _out
    return


@app.cell(hide_code=True)
def _(brentq, np):
    def acceptance(factor, x0, span0, threshold=0.5):
        """Full width of factor(x) at `threshold`, centred on the peak at x0.

        factor is sinc^2, so it falls monotonically from 1 either side of x0
        until its first zero: bracket that and solve, rather than stepping.
        span0 is a starting half-width, grown until the value drops below the
        threshold so the root is bracketed.
        """
        def f(u):
            return factor(u) - threshold

        widths = []
        for _sign in (+1, -1):
            span = span0
            for _ in range(60):                      # grow until sign change
                if f(x0 + _sign * span) < 0:
                    break
                span *= 1.6
            else:
                return np.nan
            widths.append(abs(brentq(f, x0, x0 + _sign * span, xtol=span0 * 1e-6) - x0))
        return sum(widths)
    return (acceptance,)


@app.cell(hide_code=True)
def _(mo):
    threshold = mo.ui.slider(
        0.1, 0.9, value=0.5, step=0.05, label="Acceptance threshold", show_value=True
    )
    threshold
    return (threshold,)


@app.cell(hide_code=True)
def _(
    ANGLE_PM, HAS_PM, L, POLS, T_C, WL1, WL2, acceptance, mo, np, shg, threshold, x,
):
    if HAS_PM:
        _th = threshold.value

        _d_ang = acceptance(
            lambda a: float(x.pmFactor_sfg(WL1, WL2, a, T_C, *POLS, L.value)),
            ANGLE_PM, 1e-5, _th)
        _d_T = acceptance(
            lambda t: float(x.pmFactor_sfg(WL1, WL2, ANGLE_PM, t, *POLS, L.value)),
            T_C, 0.05, _th)
        if shg.value:
            _fn = lambda w: float(x.pmFactor_sfg(w, w, ANGLE_PM, T_C, *POLS, L.value))
        else:
            _fn = lambda w: float(x.pmFactor_sfg(w, WL2, ANGLE_PM, T_C, *POLS, L.value))
        _d_wl = acceptance(_fn, WL1, 1e-5, _th)

        _fmt = lambda v, s, u: f"**{v * s:.4g}** {u}" if np.isfinite(v) else "—"
        _out = mo.md(
            f"""
    ## Acceptance bandwidths

    Full width at {_th:.2f} of sinc²(ΔkL/2), for L = {L.value:g} mm.
    Products with the crystal length, which are length-independent, are the
    figures usually quoted.

    | | full width | × L |
    |---|---|---|
    | angle {"(external ≈ n× this)" if False else ""} | {_fmt(_d_ang, 1e6, "µrad")} | {_fmt(_d_ang * L.value, 1e6, "µrad·mm")} |
    | wavelength λ₁ {"(both, λ₁=λ₂)" if shg.value else ""} | {_fmt(_d_wl, 1e6, "pm")} | {_fmt(_d_wl * L.value, 1e6, "pm·mm")} |
    | temperature | {_fmt(_d_T, 1, "°C")} | {_fmt(_d_T * L.value, 1, "°C·mm")} |
    """
        )
    else:
        _out = mo.md("")
    _out
    return


@app.cell(hide_code=True)
def _(ANGLE_PM, HAS_PM, POLS, T_C, WL1, WL2, WL3, mo, np, x):
    if HAS_PM:
        _rows = []
        for _lbl, _wl, _p in (("λ₁", WL1, POLS[0]), ("λ₂", WL2, POLS[1]), ("λ₃", WL3, POLS[2])):
            _rho = float(x.woa_theta(_wl, ANGLE_PM, T_C, pol=_p))
            _n = float(x.n(_wl, ANGLE_PM, T_C, pol=_p))
            _rows.append(f"| {_lbl} = {_wl:.4f} µm | {_p} | {_n:.5f} | {np.degrees(_rho):.4f}° |")
        _out = mo.md(
            "## Refractive indices and walk-off\n\n"
            "at the phase-matching angle. The ordinary ray of a uniaxial crystal "
            "has no walk-off.\n\n"
            "| wave | pol | n | walk-off |\n|---|---|---|---|\n" + "\n".join(_rows)
        )
    else:
        _out = mo.md("")
    _out
    return


@app.cell(hide_code=True)
def _(mo, x):
    # For a 3m crystal d_eff carries a sin(3φ) or cos(3φ) factor, so the
    # azimuthal cut sets it — even though the refractive index does not depend
    # on φ at all, which is why ndispers reports phi_rad as 'arb'.
    phi_cut = mo.ui.slider(
        0, 90, value=30, step=1, label="Azimuthal cut φ (deg)", show_value=True
    )
    phi_cut if hasattr(x, "deff_sfg") and x.phi_rad == "arb" else mo.md("")
    return (phi_cut,)


@app.cell(hide_code=True)
def _(ANGLE_PM, HAS_PM, POLS, T_C, WL1, WL2, mo, np, phi_cut, plt, x):
    _has_deff = hasattr(x, "deff_sfg")
    if HAS_PM and _has_deff:
        _phi = np.radians(phi_cut.value) if x.phi_rad == "arb" else float(x.phi_rad)
        _a = np.linspace(0.02, np.pi / 2 - 0.02, 300)
        _deff = np.array([
            abs(float(x.deff_sfg(WL1, WL2, _t, _phi, T_C, *POLS))) for _t in _a
        ])
        _at_pm = abs(float(x.deff_sfg(WL1, WL2, ANGLE_PM, _phi, T_C, *POLS)))

        _fig, _ax = plt.subplots(figsize=(7, 3.2))
        _ax.plot(np.degrees(_a), _deff)
        _ax.axvline(np.degrees(ANGLE_PM), color="r", lw=0.8, ls="--")
        _ax.set_xlabel("angle (deg)")
        _ax.set_ylabel(r"|$d_\mathrm{eff}$| (pm/V)")
        _ax.grid(alpha=0.3)
        plt.tight_layout()
        _out = mo.vstack([
            mo.md(f"## Effective nonlinearity\n\n"
                  f"|d_eff| at the phase-matching angle: **{_at_pm:.4g} pm/V**, "
                  f"for the φ = {np.degrees(_phi):.0f}° cut. Scaled from the "
                  "1.064 µm value by Miller's rule."),
            _fig,
        ])
    elif HAS_PM:
        _out = mo.md(
            "## Effective nonlinearity\n\n"
            f"*`{type(x).__name__}` carries no second-order coefficients, so "
            "d_eff cannot be computed. Available for the β-BBO parameterisations, "
            "KBBF, SLN and SLT.*"
        )
    else:
        _out = mo.md("")
    _out
    return


@app.cell(hide_code=True)
def _(mo, x):
    mo.accordion(
        {"Sellmeier equation, validity range and references":
            mo.md("```\n" + (x.__doc__ or "").strip() + "\n```")}
    )
    return


if __name__ == "__main__":
    app.run()
