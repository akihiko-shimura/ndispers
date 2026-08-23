# /// script
# requires-python = ">=3.11"
# dependencies = ["ndispers==0.17.0", "numpy", "plotly"]
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
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    from ndispers.helper import brentq
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
    return XTALS, brentq, go, make_subplots, nd, np


@app.cell(hide_code=True)
def _(XTALS, mo, nd):
    mo.md(
        f"""
    # Phase-matching calculator

    Collinear three-wave mixing — sum-frequency generation, or difference-frequency
    generation / optical parametric amplification and oscillation — in {len(XTALS)}
    birefringent crystals: phase-matching angles, tuning curves, acceptance
    bandwidths, walk-off and effective nonlinearity, computed in your browser. Powered by
    [ndispers {nd.__version__}](https://github.com/akihiko-shimura/ndispers).
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    process = mo.ui.radio(
        options=["SFG / SHG", "DFG / OPA / OPO"], value="SFG / SHG", label="Process", inline=True)
    process
    return (process,)


@app.cell(hide_code=True)
def _(XTALS, mo, process):
    IS_DFG = process.value.startswith("DFG")
    xtal = mo.ui.dropdown(options=list(XTALS), value="BetaBBO_Eimerl1987", label="Crystal")
    T = mo.ui.number(-50, 300, value=25, step=1, label="Temperature (°C)")
    L = mo.ui.number(0.1, 100.0, value=10.0, step=0.1, label="Crystal length (mm)")
    # SFG: the two inputs
    shg = mo.ui.checkbox(value=True, label="Second harmonic (λ₁ = λ₂)")
    wl1 = mo.ui.number(200.0, 4000.0, value=1064.0, step=0.1, label="λ₁ (nm)")
    wl2 = mo.ui.number(200.0, 4000.0, value=1064.0, step=0.1, label="λ₂ (nm)")
    # DFG / OPA / OPO: the pump and one of the two generated waves
    wlp = mo.ui.number(200.0, 4000.0, value=800.0, step=0.1, label="λ_pump (nm)")
    given = mo.ui.radio(options=["signal", "idler"], value="signal", label="given wave", inline=True)
    wlg = mo.ui.number(200.0, 20000.0, value=1300.0, step=0.1, label="λ_signal or λ_idler (nm)")
    return IS_DFG, L, T, given, shg, wl1, wl2, wlg, wlp, xtal


@app.cell(hide_code=True)
def _(IS_DFG, L, T, given, mo, shg, wl1, wl2, wlg, wlp, xtal):
    if IS_DFG:
        _mid = mo.vstack([wlp, wlg, given], gap=0.25)
    else:
        _mid = mo.vstack([wl1, wl2 if not shg.value else mo.md("<sub>λ₂ = λ₁</sub>")])
    mo.hstack(
        [
            mo.vstack([xtal] + ([] if IS_DFG else [shg]), gap=0.25),
            _mid,
            mo.vstack([T, L], gap=0.25),
        ],
        justify="start",
        gap=0.75,
    )
    return


@app.cell(hide_code=True)
def _(IS_DFG, XTALS, T, given, mo, shg, wl1, wl2, wlg, wlp, xtal):
    x = XTALS[xtal.value]()
    # ndispers works in µm; the interface is in nm. Internally everything is
    # the SFG of waves 1 and 2 into 3; for DFG / OPA / OPO those are
    # (signal, idler, pump), so every cell below serves both processes.
    if IS_DFG:
        mo.stop(wlg.value <= wlp.value,
                mo.md("*The signal or idler must be longer than the pump — "
                      "1/λ_p = 1/λ_s + 1/λ_i.*"))
        WL3 = wlp.value * 1e-3
        _g = wlg.value * 1e-3
        _other = x.wl_idler(WL3, _g)
        WL1, WL2 = (_g, _other) if given.value == "signal" else (_other, _g)
        WAVES = ("signal", "idler", "pump")
    else:
        WL1 = wl1.value * 1e-3
        WL2 = WL1 if shg.value else wl2.value * 1e-3
        WL3 = 1.0 / (1.0 / WL1 + 1.0 / WL2)
        WAVES = ("λ₁", "λ₂", "λ₃")
    T_C = T.value
    ANGLE_KEY = "theta" if x.theta_rad == "var" else "phi"
    return ANGLE_KEY, T_C, WAVES, WL1, WL2, WL3, x


@app.cell(hide_code=True)
def _(mo, x):
    # the docstring is preformatted, but its reference lines are long; render it
    # in a wrapping block rather than a code fence so it does not scroll sideways
    mo.md(
        '<div style="white-space:pre-wrap; overflow-wrap:anywhere; '
        'font-family:ui-monospace,monospace; font-size:0.82em; line-height:1.5">'
        + (x.__doc__ or "").strip().replace("&", "&amp;").replace("<", "&lt;")
        + "</div>"
    )
    return


@app.cell(hide_code=True)
def _(HAS_PM, POLS, mo, x):
    if HAS_PM:
        _md = "\n\n".join(
            f"$$n_{{{_p}}} = {x.n_latex(_p)}$$"
            for _p in dict.fromkeys(POLS)
        )
        _out = mo.accordion({
            "Sellmeier equation as evaluated": mo.md(
                "What ndispers differentiates and evaluates, with the "
                "coefficients already substituted — not a quotation from the "
                "paper. **λ is in µm and T in °C in this expression**, unlike "
                "the nm inputs above.\n\n" + _md
            )
        })
    else:
        _out = mo.md("")
    _out
    return


@app.cell(hide_code=True)
def _(ANGLE_KEY, IS_DFG, T_C, WAVES, WL1, WL2, WL3, mo, x):
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
        + (f"**{SOLUTIONS[k]:.3f}°** | {SOLUTIONS[k] * 3.141592653589793 / 180:.5f} rad"
           if k in SOLUTIONS else "— none | —")
        + " |"
        for k in TYPE_OF
    )
    mo.md(
        f"""
    ## Phase-matching angles

    {(f"pump {WL3 * 1e3:.2f} nm → signal {WL1 * 1e3:.2f} nm + **idler {WL2 * 1e3:.2f} nm**"
      if IS_DFG else
      f"λ₁ = {WL1 * 1e3:.2f} nm, λ₂ = {WL2 * 1e3:.2f} nm → **λ₃ = {WL3 * 1e3:.2f} nm**")}
    &nbsp;·&nbsp; the varying angle of this crystal is **{ANGLE_KEY}**

    | polarizations ({", ".join(WAVES) if IS_DFG else "1,2,3"}) | type | angle | |
    |---|---|---|---|
    {_rows}
    """
    )
    return SOLUTIONS, TYPE_OF


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
def _(ANGLE_KEY, SOLUTIONS, T_C, np, pmtype, xtal):
    HAS_PM = bool(SOLUTIONS) and pmtype.value in SOLUTIONS
    ANGLE_PM = np.radians(SOLUTIONS[pmtype.value]) if HAS_PM else np.nan
    POLS = tuple(pmtype.value) if HAS_PM else ("o", "o", "e")
    # the state every plot below was computed at
    PLOT_TITLE = (
        f"{xtal.value}   T = {T_C:g} °C, {ANGLE_KEY} = "
        f"{SOLUTIONS[pmtype.value]:.3f}°, pol = {pmtype.value}"
        if HAS_PM else ""
    )
    PLOT_LAYOUT = dict(
        title=dict(font_size=11, x=0, xanchor="left", y=0.97, yanchor="top"),
        width=560, margin=dict(l=60, r=20, t=40, b=45), font_size=11,
        template="plotly_white", showlegend=False, modebar=dict(orientation="v"),
    )
    return ANGLE_PM, HAS_PM, PLOT_LAYOUT, PLOT_TITLE, POLS


@app.cell(hide_code=True)
def _(ANGLE_PM, HAS_PM, L, PLOT_LAYOUT, PLOT_TITLE, POLS, T_C, WL1, WL2, make_subplots, mo, np, x):
    if HAS_PM:
        _w = np.radians(0.15)
        _a = np.linspace(ANGLE_PM - _w, ANGLE_PM + _w, 400)
        _dk = x.dk_sfg(WL1, WL2, _a, T_C, *POLS)
        _sinc = x.pmFactor_sfg(WL1, WL2, _a, T_C, *POLS, L.value)

        _deg = np.degrees(_a)
        _fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.06)
        _fig.add_scatter(x=_deg, y=_dk, mode="lines", row=1, col=1,
                         hovertemplate="%{x:.4f}°<br>Δk %{y:.5g} rad/µm<extra></extra>")
        _fig.add_scatter(x=_deg, y=_sinc, mode="lines", row=2, col=1,
                         hovertemplate="%{x:.4f}°<br>sinc² %{y:.4f}<extra></extra>")
        _fig.add_hline(y=0, line=dict(color="gray", width=0.8), row=1, col=1)
        _fig.add_vline(x=np.degrees(ANGLE_PM), line=dict(color="red", width=0.8, dash="dash"))
        _fig.update_yaxes(title_text="Δk (rad/µm)", row=1, col=1)
        _fig.update_yaxes(title_text="sinc²(ΔkL/2)", row=2, col=1)
        _fig.update_xaxes(title_text="angle (deg)", row=2, col=1)
        _fig.update_layout(**PLOT_LAYOUT, title_text=PLOT_TITLE, height=360)
        _out = mo.vstack([mo.md("## Phase mismatch and conversion"), _fig])
    else:
        _out = mo.md("")
    _out
    return


@app.cell(hide_code=True)
def _(ANGLE_KEY, ANGLE_PM, HAS_PM, IS_DFG, PLOT_LAYOUT, POLS, T_C, WL1, WL2, WL3, go, make_subplots, mo, np, x):
    # OPO / OPA tuning curves: the zero contour of dk over (angle, signal
    # frequency) and over (temperature, signal frequency) - one vectorised
    # evaluation each, every branch drawn, no root finder needed for a plot.
    if IS_DFG and HAS_PM:
        import re as _re
        # longest idler to draw: the crystal's transparency edge (the Sellmeier
        # validity range is usually much narrower than where an OPO actually tunes)
        _doc = type(x).__doc__ or ""
        _sec = (_re.search(r"Transparency range\s*:([^\n]*)", _doc)
                or _re.search(r"Validity range\s*\n\s*-+\s*\n(.{0,400})", _doc, _re.S))
        _m = _sec and _re.search(r"(\d+\.?\d*)\s*(?:µm|um)?\s*(?:to|-|–)\s*(\d+\.?\d*)", _sec.group(1))
        _wl_i_max = min(float(_m.group(2)), 20.0) if _m else 5.0
        _wl_i_max = max(_wl_i_max, 1.2 * WL2)          # the current idler must be inside
        _nu_p = 1.0 / WL3
        _nu_s = np.linspace(_nu_p / 2, _nu_p - 1.0 / _wl_i_max, 400)   # even in frequency
        _wl_s = 1.0 / _nu_s
        _wl_i = 1.0 / (_nu_p - _nu_s)
        _ang = np.linspace(0.0, 90.0, 361)
        _T = np.linspace(-50.0, 300.0, 141)
        with np.errstate(all="ignore"):
            _dk_ang = np.asarray(x.dk_dfg(WL3, _wl_s[None, :], np.radians(_ang)[:, None], T_C, *POLS), float)
            _dk_T = np.asarray(x.dk_dfg(WL3, _wl_s[None, :], ANGLE_PM, _T[:, None], *POLS), float)
        def _zero_curve(xs, z):
            """(x, index) points where z changes sign along its second axis, with
            the crossing interpolated linearly - the dk = 0 locus as a sparse
            point set (a contour trace would ship the whole grid to the browser)."""
            _px, _py = [], []
            for _i, _row in enumerate(z):
                _ok = np.isfinite(_row)
                for _j in np.nonzero(np.diff(np.signbit(_row)) & _ok[:-1] & _ok[1:])[0]:
                    _f = _row[_j] / (_row[_j] - _row[_j + 1])
                    _px.append(xs[_i]); _py.append(_j + _f)
            return np.array(_px), np.array(_py)

        _fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.12,
                             subplot_titles=(f"angle tuning at T = {T_C:g} °C",
                                             f"temperature tuning at {ANGLE_KEY} = {np.degrees(ANGLE_PM):.3f}°"))
        _idx = np.arange(_wl_s.size)
        for _col, _xs, _dk in ((1, _ang, _dk_ang), (2, _T, _dk_T)):
            _px, _pj = _zero_curve(_xs, _dk)
            for _y, _name, _color in ((_wl_s * 1e3, "signal", "#1f77b4"), (_wl_i * 1e3, "idler", "#d62728")):
                _fig.add_scatter(
                    x=_px, y=np.interp(_pj, _idx, _y) if _px.size else [], mode="markers",
                    marker=dict(color=_color, size=3), row=1, col=_col,
                    name=_name, showlegend=(_col == 1),
                    hovertemplate=f"%{{x:.3f}}<br>{_name} %{{y:.1f}} nm<extra></extra>")
        _fig.add_scatter(x=[np.degrees(ANGLE_PM)], y=[WL1 * 1e3], mode="markers",
                         marker=dict(color="#1f77b4", size=8), row=1, col=1, showlegend=False,
                         hovertemplate="current signal<extra></extra>")
        _fig.add_scatter(x=[np.degrees(ANGLE_PM)], y=[WL2 * 1e3], mode="markers",
                         marker=dict(color="#d62728", size=8), row=1, col=1, showlegend=False,
                         hovertemplate="current idler<extra></extra>")
        _fig.add_scatter(x=[T_C], y=[WL1 * 1e3], mode="markers", marker=dict(color="#1f77b4", size=8),
                         row=1, col=2, showlegend=False, hovertemplate="current signal<extra></extra>")
        _fig.add_scatter(x=[T_C], y=[WL2 * 1e3], mode="markers", marker=dict(color="#d62728", size=8),
                         row=1, col=2, showlegend=False, hovertemplate="current idler<extra></extra>")
        _fig.update_xaxes(title_text=f"{ANGLE_KEY} (deg)", range=[0, 90], row=1, col=1)
        _fig.update_xaxes(title_text="T (°C)", range=[-50, 300], row=1, col=2)
        _fig.update_yaxes(title_text="wavelength (nm)", range=[WL3 * 1e3, _wl_i_max * 1e3], row=1, col=1)
        _fig.update_yaxes(range=[WL3 * 1e3, _wl_i_max * 1e3], row=1, col=2)
        _fig.update_layout(**{**PLOT_LAYOUT, "showlegend": True, "width": 760}, height=340,
                           legend=dict(x=0.01, y=0.99, bgcolor="rgba(255,255,255,0.6)"))
        _fig.update_annotations(font_size=11)          # subplot titles at axis-label size
        _out = mo.vstack([mo.md(
            "## Tuning curves\n\n"
            f"Signal and idler that phase-match for pump {WL3 * 1e3:.2f} nm and polarizations "
            f"`{''.join(POLS)}`, as the angle or the temperature is tuned; the dots are the present "
            "operating point. Curves are drawn up to the crystal's transparency edge "
            f"({_wl_i_max:g} µm); outside the Sellmeier validity range they extrapolate the fit. "
            "A medium without a temperature term gives a flat temperature curve."),
            _fig])
    else:
        _out = mo.md("")
    _out
    return


@app.cell(hide_code=True)
def _(IS_DFG, PLOT_LAYOUT, T_C, WL1, WL3, go, mo, np, x):
    # Quasi-phase-matching: the poling period that phase-matches each signal
    # for all waves extraordinary along a principal axis (d33 of LiNbO3, LiTaO3,
    # KTP). Drawn for every crystal with an e-ray; meaningful for the ferroelectrics.
    if IS_DFG:
        try:
            _nu_p = 1.0 / WL3
            _nu_s = np.linspace(_nu_p / 2, _nu_p - 1.0 / min(5.5, 12 * WL3), 300)
            _wl_s = 1.0 / _nu_s
            _wl_i = 1.0 / (_nu_p - _nu_s)
            with np.errstate(all="ignore"):
                _per = np.asarray(x.qpm_period_dfg(WL3, _wl_s, np.pi / 2, T_C, "e", "e", "e"), float)
            _per_now = float(x.qpm_period_dfg(WL3, WL1, np.pi / 2, T_C, "e", "e", "e"))
            _fig = go.Figure(go.Scatter(
                x=_wl_s * 1e3, y=_per, mode="lines", customdata=_wl_i * 1e3,
                hovertemplate="signal %{x:.1f} nm, idler %{customdata:.1f} nm<br>Λ = %{y:.2f} µm<extra></extra>"))
            _fig.add_scatter(x=[WL1 * 1e3], y=[_per_now], mode="markers", marker=dict(size=8, color="#d62728"),
                             hovertemplate="current signal<extra></extra>")
            _fig.update_layout(**PLOT_LAYOUT, height=280, xaxis_title="signal (nm)",
                               yaxis_title="poling period Λ (µm)",
                               title_text=f"QPM, e+e→e along the axis, T = {T_C:g} °C, pump {WL3 * 1e3:.2f} nm")
            _out = mo.vstack([mo.md(
                "## Quasi-phase-matching period\n\n"
                "First-order poling period for all three waves extraordinary, propagating along "
                "a principal axis (the PPLN / PPSLT / PPKTP geometry, using d33): "
                f"**Λ = {_per_now:.2f} µm** for the present signal. Temperature tuning follows the "
                "medium's dn/dT — choose a parameterisation with one (SLN, SLT) for that."),
                _fig])
        except ValueError:
            _out = mo.md("")
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
        """
        def f(u):
            return factor(u) - threshold

        widths = []
        for _sign in (+1, -1):
            span = span0
            for _ in range(80):
                if f(x0 + _sign * span) < 0:
                    break
                span *= 1.6
            else:
                return np.nan
            widths.append(abs(brentq(f, x0, x0 + _sign * span, xtol=span0 * 1e-8) - x0))
        return sum(widths)
    return (acceptance,)


@app.cell(hide_code=True)
def _(mo):
    threshold = mo.ui.slider(
        0.1, 0.9, value=0.5, step=0.05, label="Acceptance threshold", show_value=True
    )
    # scale from ndispers' native unit (µm for wavelength, rad for angle)
    WL_UNITS = {"pm": 1e6, "nm": 1e3, "µm": 1.0}
    ANG_UNITS = {"mdeg": 180e3 / 3.141592653589793, "mrad": 1e3,
                 "µrad": 1e6, "deg": 180 / 3.141592653589793}
    wl_unit = mo.ui.dropdown(options=list(WL_UNITS), value="pm", label="Wavelength unit")
    ang_unit = mo.ui.dropdown(options=list(ANG_UNITS), value="mrad", label="Angle unit")
    mo.hstack([threshold, wl_unit, ang_unit], justify="start", gap=0.75)
    return ANG_UNITS, WL_UNITS, ang_unit, threshold, wl_unit


@app.cell(hide_code=True)
def _(
    ANGLE_PM, ANG_UNITS, HAS_PM, IS_DFG, L, POLS, T_C, WL1, WL2, WL3, WL_UNITS,
    acceptance, ang_unit, mo, np, shg, threshold, wl_unit, x,
):
    if HAS_PM:
        _th = threshold.value
        _L = L.value
        _au, _as = ang_unit.value, ANG_UNITS[ang_unit.value]
        _wu, _ws = wl_unit.value, WL_UNITS[wl_unit.value]

        _d_ang = acceptance(
            lambda a: float(x.pmFactor_sfg(WL1, WL2, a, T_C, *POLS, _L)),
            ANGLE_PM, 1e-6, _th)
        _d_T = acceptance(
            lambda t: float(x.pmFactor_sfg(WL1, WL2, ANGLE_PM, t, *POLS, _L)),
            T_C, 0.01, _th)
        if IS_DFG:
            # pump fixed: tuning the signal moves the idler the other way
            # (energy conservation) - the parametric gain bandwidth, whose
            # leading term is the signal-idler group-velocity mismatch
            _d_sig = acceptance(
                lambda w: float(x.pmFactor_dfg(WL3, w, ANGLE_PM, T_C, *POLS, _L)),
                WL1, 1e-6, _th)
            # signal fixed, pump tuned (idler follows): the pump bandwidth
            _d_pump = acceptance(
                lambda w: float(x.pmFactor_dfg(w, WL1, ANGLE_PM, T_C, *POLS, _L)),
                WL3, 1e-6, _th)
        else:
            # λ₂ held at its stated value, only λ₁ tuned — the SNLO convention
            _d_wl1 = acceptance(
                lambda w: float(x.pmFactor_sfg(w, WL2, ANGLE_PM, T_C, *POLS, _L)),
                WL1, 1e-6, _th)
            # both tuned together; meaningful for SHG of a broadband source
            _d_wboth = acceptance(
                lambda w: float(x.pmFactor_sfg(w, w, ANGLE_PM, T_C, *POLS, _L)),
                WL1, 1e-6, _th) if shg.value else np.nan

        def _row(label, value, scale, unit):
            """value is in ndispers' native unit; L-product is per centimetre."""
            if not np.isfinite(value):
                return f"| {label} | — | — |"
            return (f"| {label} | **{value * scale:.4g}** {unit} | "
                    f"**{value * scale * _L * 0.1:.4g}** {unit}·cm |")

        # (label, width in the chosen unit, unit) - the report cell reads this too
        ACCEPT = [
            ("angle", _d_ang * _as, _au),
            ("temperature", _d_T, "°C"),
        ]
        if IS_DFG:
            ACCEPT.append(("signal (pump fixed, idler follows)", _d_sig * _ws, _wu))
            ACCEPT.append(("pump (signal fixed, idler follows)", _d_pump * _ws, _wu))
        else:
            ACCEPT.append(("λ₁ only, λ₂ fixed", _d_wl1 * _ws, _wu))
            if shg.value:
                ACCEPT.append(("λ₁ and λ₂ tuned together", _d_wboth * _ws, _wu))
        _rows = [_row(_lbl, _v, 1, _u) for _lbl, _v, _u in ACCEPT]

        _out = mo.md(
            f"""
    ## Acceptance bandwidths

    Full width at {_th:.2f} of sinc²(ΔkL/2) for L = {_L:g} mm, and the
    length product, which is what is usually tabulated.

    | | full width at L = {_L:g} mm | × L |
    |---|---|---|
    """ + "\n".join("    " + r for r in _rows) + """

    """ + (
    """The **signal** row is the parametric gain bandwidth: the pump is fixed and
    the idler follows the signal through 1/λ_i = 1/λ_p − 1/λ_s, so to first order
    it is set by the signal–idler group-velocity mismatch — it opens up towards
    degeneracy, and fully in a noncollinear geometry. This is a different quantity
    from SNLO's SFG acceptance, which tunes one input with the other fixed.
    """ if IS_DFG else
    """The wavelength row to compare against tabulated values is **λ₁ only** —
    that is what SNLO reports. Tuning λ₁ and λ₂ together makes Δk move twice
    as fast, so that width is exactly half; it is the relevant one when a
    single broadband beam supplies both photons.
    """)
        )
    else:
        ACCEPT = []
        _out = mo.md("")
    _out
    return (ACCEPT,)


@app.cell(hide_code=True)
def _(ANGLE_PM, HAS_PM, POLS, T_C, WAVES, WL1, WL2, WL3, mo, np, x):
    if HAS_PM:
        _rows = []
        for _lbl, _wl, _p in zip(WAVES, (WL1, WL2, WL3), POLS):
            _rho = float(x.woa_theta(_wl, ANGLE_PM, T_C, pol=_p))
            _n = float(x.n(_wl, ANGLE_PM, T_C, pol=_p))
            _rows.append(
                f"| {_lbl} = {_wl * 1e3:.2f} nm | {_p} | {_n:.5f} | "
                f"{np.degrees(_rho):.4f}° | {_rho * 1e3:.4g} mrad |"
            )
        _out = mo.md(
            "## Refractive indices and walk-off\n\n"
            "at the phase-matching angle. The ordinary ray of a uniaxial crystal "
            "has no walk-off, and the three indices coinciding is the "
            "phase-matching condition itself.\n\n"
            "| wave | pol | n | walk-off | |\n|---|---|---|---|---|\n" + "\n".join(_rows)
        )
    else:
        _out = mo.md("")
    _out
    return


@app.cell(hide_code=True)
def _(ANGLE_PM, HAS_PM, POLS, T_C, WL1, WL2, np, x):
    # For a uniaxial crystal the refractive index does not depend on φ (ndispers
    # reports phi_rad as 'arb'), but |d_eff| does through sin 3φ, cos 2φ, ...,
    # so the azimuthal cut is a free choice. A principal plane of a biaxial
    # crystal fixes one angle; deff_sfg takes None for that one.
    HAS_DEFF = hasattr(x, "deff_sfg")
    PHI_FREE = HAS_DEFF and x.phi_rad == "arb"

    def deff_args(angle, phi):
        """(theta, phi) for deff_sfg: the swept angle plus the fixed/chosen one."""
        if x.theta_rad == "var":
            return (angle, phi if PHI_FREE else None)
        return (None, angle)

    if HAS_PM and PHI_FREE:
        _p = np.linspace(0, 180, 1801)
        _d = np.abs(x.deff_sfg(WL1, WL2, ANGLE_PM, np.radians(_p), T_C, *POLS))
        _peak = _d.max()
        _best = _p[_d > _peak * 0.9999]
        PHI_OPT = [_best[0]]
        for _v in _best[1:]:
            if _v - PHI_OPT[-1] > 1.0:
                PHI_OPT.append(_v)
        DEFF_MAX = _peak
    else:
        PHI_OPT, DEFF_MAX = [], np.nan
    return DEFF_MAX, HAS_DEFF, PHI_FREE, PHI_OPT, deff_args


@app.cell(hide_code=True)
def _(PHI_FREE, PHI_OPT, mo):
    phi_cut = mo.ui.slider(
        0, 180, value=int(round(PHI_OPT[-1])) if PHI_OPT else 0, step=1,
        label="Azimuthal cut φ (deg)", show_value=True,
    )
    phi_cut if PHI_FREE else mo.md("")
    return (phi_cut,)


@app.cell(hide_code=True)
def _(
    ANGLE_PM, DEFF_MAX, HAS_DEFF, HAS_PM, PHI_FREE, PHI_OPT, PLOT_LAYOUT, PLOT_TITLE, POLS,
    T_C, WL1, WL2, deff_args, go, mo, np, phi_cut, x,
):
    # closed forms per point group (walk-off angles ρ_m of the e-waves included);
    # the library evaluates the tensor contraction directly, these are for reading
    _EXPR = {
        ("Uniax_3m", "ooe"): r"d_{31}\sin(\theta+\rho_3) - d_{22}\cos(\theta+\rho_3)\sin 3\phi",
        ("Uniax_3m", "oee"): r"d_{22}\cos(\theta+\rho_2)\cos(\theta+\rho_3)\cos 3\phi",
        ("Uniax_3m", "eoe"): r"d_{22}\cos(\theta+\rho_1)\cos(\theta+\rho_3)\cos 3\phi",
        ("Uniax_32", "ooe"): r"d_{11}\cos(\theta+\rho_3)\cos 3\phi",
        ("Uniax_32", "oee"): r"d_{11}\cos(\theta+\rho_2)\cos(\theta+\rho_3)\sin 3\phi",
        ("Uniax_32", "eoe"): r"d_{11}\cos(\theta+\rho_1)\cos(\theta+\rho_3)\sin 3\phi",
        ("Uniax_42m", "ooe"): r"d_{36}\sin(\theta+\rho_3)\sin 2\phi",
        ("Uniax_42m", "oee"): r"d_{36}\sin(2\theta+\rho_2+\rho_3)\cos 2\phi",
        ("Uniax_42m", "eoe"): r"d_{36}\sin(2\theta+\rho_1+\rho_3)\cos 2\phi",
        ("Uniax_4mm", "ooe"): r"d_{31}\sin(\theta+\rho_3)",
    }
    _group = type(x).__mro__[1].__name__
    _pols = "".join(POLS)
    if HAS_PM and HAS_DEFF:
        _phi = np.radians(phi_cut.value) if PHI_FREE else None
        try:
            _at = abs(float(x.deff_sfg(WL1, WL2, *deff_args(ANGLE_PM, _phi), T_C, *POLS)))
            # point group 2 (BiBO): theta and pi - theta phase-match alike but
            # d_eff differs, so sweep the full half-circle there
            _amax = np.pi if _group == "Biax_2" else np.pi / 2
            _a = np.linspace(0.02, _amax - 0.02, 200 if _amax < 2 else 400)
            _curve = np.abs(x.deff_sfg(WL1, WL2, *deff_args(_a, _phi), T_C, *POLS))
            _fig = go.Figure(go.Scatter(
                x=np.degrees(_a), y=_curve, mode="lines",
                hovertemplate="%{x:.2f}°<br>|d_eff| %{y:.4g} pm/V<extra></extra>"))
            _fig.add_vline(x=np.degrees(ANGLE_PM), line=dict(color="red", width=0.8, dash="dash"))
            if _group == "Biax_2":
                _fig.add_vline(x=180 - np.degrees(ANGLE_PM), line=dict(color="red", width=0.8, dash="dot"))
            _fig.update_layout(
                **PLOT_LAYOUT, height=260,
                title_text=PLOT_TITLE + (f", φ = {phi_cut.value:.0f}°" if PHI_FREE else ""),
                xaxis_title="angle (deg)", yaxis_title="|d_eff| (pm/V)")

            _forbidden = float(np.max(_curve)) < 1e-9 * max(
                abs(x.d_sfg(il, WL1, WL2, T_C)) for il in x._d_ref)
            _opt = (
                " &nbsp;·&nbsp; |d_eff| peaks at **φ = "
                + ", ".join(f"{v:.0f}°" for v in PHI_OPT)
                + f"** with {DEFF_MAX:.4g} pm/V"
                if PHI_OPT else ""
            )
            _expr = _EXPR.get((_group, _pols))
            _formula = (f"$$|d_\\mathrm{{eff}}| = \\left|{_expr}\\right|$$\n\n" if _expr else
                        "d_eff is the contraction of the d tensor with the field "
                        "directions of the three waves, walk-off included.\n\n")
            # the tensor components, at their reference wavelengths and scaled to these
            _rows = "\n".join(
                f"| {il} | {d0:g} | {wl1 * 1e3:g} + {wl2 * 1e3:g} | "
                f"{x.d_sfg(il, WL1, WL2, T_C):.4g} |"
                for il, (d0, wl1, wl2) in x._d_ref.items())
            _out = mo.vstack([
                mo.md(
                    "## Effective nonlinearity\n\n" + _formula
                    + f"At the phase-matching angle{f' and φ = {phi_cut.value:.0f}°' if PHI_FREE else ''}: "
                    f"**{_at:.4g} pm/V**{_opt}.\n\n"
                    + ("This interaction phase-matches but its d_eff vanishes "
                       "identically: the tensor has no component that couples "
                       f"these polarizations in the {x.plane} plane, so no light "
                       "is generated. Choose another interaction or plane.\n\n"
                       if _forbidden else "")
                    + ("In a monoclinic crystal the mirror angle 180° − θ phase-matches "
                       "too but with a different d_eff: "
                       f"**{abs(float(x.deff_sfg(WL1, WL2, *deff_args(np.pi - ANGLE_PM, _phi), T_C, *POLS))):.4g} pm/V** "
                       "(dotted line).\n\n" if _group == "Biax_2" else "")
                    + "| component | reference (pm/V) | measured at (nm) | scaled to these wavelengths |\n"
                    "|---|---|---|---|\n" + _rows + "\n\n"
                    f"<sub>Reference values and how far to trust the Miller-rule scaling: {x._d_note}</sub>"
                ),
                _fig,
            ])
        except ValueError as _e:
            _out = mo.md(f"## Effective nonlinearity\n\n*Not available for `{_pols}`: {_e}*")
    elif HAS_PM:
        _out = mo.md(
            "## Effective nonlinearity\n\n"
            f"*`{type(x).__name__}` is centrosymmetric: no second-order nonlinearity.*"
        )
    else:
        _out = mo.md("")
    _out
    return


@app.cell(hide_code=True)
def _(
    ACCEPT, ANGLE_KEY, ANGLE_PM, DEFF_MAX, HAS_DEFF, HAS_PM, IS_DFG, L, PHI_FREE, PHI_OPT,
    POLS, SOLUTIONS, T_C, TYPE_OF, WAVES, WL1, WL2, WL3, deff_args, mo, nd, np, phi_cut,
    pmtype, shg, threshold, x, xtal,
):
    # A plain-text report of everything above, in the spirit of SNLO's QMIX
    # output but carrying its own provenance: which parameterisation, which
    # papers, which version - so the numbers can be traced later.
    import datetime as _dt

    def _section(doc, head):
        """The lines of a docstring section (`head` underlined with ---)."""
        _lines = [l.strip() for l in (doc or "").splitlines()]
        try:
            _i = _lines.index(head)
        except ValueError:
            return []
        _body = []
        for l in _lines[_i + 2:]:
            if l and set(l) <= {"-"}:
                _body.pop()                      # the previous line was the next heading
                break
            _body.append(l)
        while _body and not _body[-1]:
            _body.pop()
        return _body

    if HAS_PM:
        _deg = np.degrees(ANGLE_PM)
        _cut = f"theta = {_deg:.3f} deg" if ANGLE_KEY == "theta" else f"phi = {_deg:.3f} deg"
        if x.theta_rad != "var":
            _cut = f"theta = {np.degrees(x.theta_rad):g} deg, " + _cut
        elif PHI_FREE:
            _cut += f", phi = {phi_cut.value:g} deg"
        elif not isinstance(x.phi_rad, str):
            _cut += f", phi = {np.degrees(x.phi_rad):g} deg"
        _cut_note = "   (internal angles; phi is a free choice here)" if PHI_FREE else "   (internal angles)"
        _pols = "".join(POLS)
        _type = TYPE_OF[_pols]

        _waves = []
        _wn = ("signal", "idler", "pump") if IS_DFG else ("wl1", "wl2", "wl3")
        for _lbl, _wl, _p in zip(_wn, (WL1, WL2, WL3), POLS):
            _f = lambda q, _wl=_wl, _p=_p: float(getattr(x, q)(_wl, ANGLE_PM, T_C, pol=_p))
            _waves.append((_lbl, _wl, _p, _f("n"), _f("ng"), _f("GV"), _f("GVD"), _f("TOD"),
                           _f("woa_theta") * 1e3, _f("dndT")))
        _gvm = 1 / _waves[2][5] - 1 / _waves[0][5]          # fs/um, 3 vs 1
        _gvm_si = 1 / _waves[0][5] - 1 / _waves[1][5]       # 1 vs 2
        _gvm_pi = 1 / _waves[2][5] - 1 / _waves[1][5]       # 3 vs 2
        _L_um = L.value * 1e3

        # "- Point group : 3m (C3v)" in the docstring; the group mixin's name otherwise
        _pg = next((l.split(":", 1)[1].split(";")[0].strip() for l in (x.__doc__ or "").splitlines()
                    if "Point group" in l), None)
        if _pg is None:
            _pg = next((c.__name__.split("_", 1)[1] for c in type(x).__mro__
                        if c.__name__.startswith(("Uniax_", "Biax_"))), "n/a")
        _lines = [
            "ndispers phase-matching report",
            f"generated {_dt.date.today().isoformat()}  by ndispers {nd.__version__}  "
            "(akihiko-shimura.github.io/ndispers/phasematching)",
            "",
            "CRYSTAL",
            f"  {xtal.value}   point group {_pg}"
            + (f", {x.plane} plane" if x.plane != "arb" else ""),
        ]
        for l in _section(x.__doc__, "Ref"):
            if l:                                    # headings flush, references indented
                _lines.append(("  " if l.endswith(":") else "    ") + l)
        _vr = _section(x.__doc__, "Validity range")
        if _vr:
            _lines.append("  validity:   " + " ".join(_vr))
        _lines += [
            "",
            "CONDITIONS",
            ("  process        DFG / OPA / OPO   1/wl_p = 1/wl_s + 1/wl_i" if IS_DFG else
             "  process        SFG  1/wl3 = 1/wl1 + 1/wl2" + ("   (SHG: wl1 = wl2)" if shg.value else "")),
            (f"  pump, signal, idler  {WL3 * 1e3:.2f}, {WL1 * 1e3:.2f}, {WL2 * 1e3:.2f} nm   (vacuum)" if IS_DFG else
             f"  wl1, wl2, wl3  {WL1 * 1e3:.2f}, {WL2 * 1e3:.2f}, {WL3 * 1e3:.2f} nm   (vacuum)"),
            f"  temperature    {T_C:g} degC",
            f"  length         {L.value:g} mm",
            f"  polarizations  {' '.join(POLS)}   (Type {_type}; order = {', '.join(_wn)})",
            f"  cut            {_cut}{_cut_note}",
            "",
            f"PHASE MATCHING   (all solutions at this T and these wavelengths; {ANGLE_KEY} varies)",
        ]
        for _k in TYPE_OF:
            _lines.append(f"  {_k}  Type {TYPE_OF[_k]:<14s} "
                          + (f"{ANGLE_KEY} = {SOLUTIONS[_k]:.3f} deg" if _k in SOLUTIONS else "none"))
        _lines += [
            "",
            f"AT THE SELECTED SOLUTION  ({_pols}, {ANGLE_KEY} = {_deg:.3f} deg)",
            "  wave     wl(nm)  pol      n       n_g    v_g(um/fs)  GVD(fs^2/mm)  TOD(fs^3/mm)  walk-off(mrad)  dn/dT(1/K)",
        ]
        for _w in _waves:
            _lines.append(f"  {_w[0]:<6s}{_w[1] * 1e3:9.2f}   {_w[2]}   {_w[3]:8.5f} {_w[4]:8.5f}   "
                          f"{_w[5]:9.6f}   {_w[6]:11.2f}   {_w[7]:11.1f}   {_w[8]:12.3f}    {_w[9]:10.3e}")
        _n1, _n2, _n3 = _wn
        _lines.append(f"  group-velocity mismatch  1/v_g({_n3}) - 1/v_g({_n1}) = {_gvm:+.4f} fs/um   "
                      f"(over L: {_gvm * _L_um:+.0f} fs)")
        _lines.append(f"                           1/v_g({_n3}) - 1/v_g({_n2}) = {_gvm_pi:+.4f} fs/um   "
                      f"(over L: {_gvm_pi * _L_um:+.0f} fs)")
        _lines.append(f"                           1/v_g({_n1}) - 1/v_g({_n2}) = {_gvm_si:+.4f} fs/um   "
                      f"(over L: {_gvm_si * _L_um:+.0f} fs)"
                      + ("   <- sets the parametric gain bandwidth" if IS_DFG else ""))
        _lines += ["", f"ACCEPTANCE  (full width at {threshold.value:.2f} of sinc^2(dk L/2);  L = {L.value:g} mm)"]
        _tag = {"λ₁ only, λ₂ fixed": "wl1 only, wl2 fixed  <- what SNLO and tables report",
                "λ₁ and λ₂ tuned together": "wl1 = wl2 tuned together   (broadband SHG)",
                "signal (pump fixed, idler follows)": "signal, pump fixed   (idler follows; gain bandwidth)",
                "pump (signal fixed, idler follows)": "pump, signal fixed   (idler follows)"}
        for _lbl, _v, _u in ACCEPT:
            _name = _tag.get(_lbl, _lbl)
            _u2 = {"°C": "degC", "µrad": "urad", "µm": "um"}.get(_u, _u)
            if np.isfinite(_v):
                _lines.append(f"  {_name.split('  <-')[0].split('   (')[0]:<26s} {_v:10.4g} {_u2:<5s} "
                              f"{_v * L.value * 0.1:10.4g} {_u2}.cm"
                              + ("   <- " + _name.split("<- ")[1] if "<-" in _name else "")
                              + ("   " + _name[_name.index("("):] if "(" in _name else ""))
        _lines += ["", "EFFECTIVE NONLINEARITY"]
        if HAS_DEFF:
            try:
                _phi = np.radians(phi_cut.value) if PHI_FREE else None
                _de = float(x.deff_sfg(WL1, WL2, *deff_args(ANGLE_PM, _phi), T_C, *POLS))
                _lines.append(f"  d_eff = {_de:+.4g} pm/V   at {_cut}"
                              + (f"   (|d_eff| max {DEFF_MAX:.4g} pm/V at phi = "
                                 + ", ".join(f"{v:.0f}" for v in PHI_OPT) + " deg)" if PHI_OPT else ""))
                _lines.append("  " + ", ".join(
                    f"{il} = {x.d_sfg(il, WL1, WL2, T_C):.4g}" for il in x._d_ref)
                    + f" pm/V   (Miller-scaled to {WL1 * 1e3:g} + {WL2 * 1e3:g} nm from "
                    + "; ".join(f"{il} = {d0:g} @ {w1 * 1e3:g}+{w2 * 1e3:g}" for il, (d0, w1, w2) in x._d_ref.items())
                    + ")")
                _lines.append(f"  note: {x._d_note}")
            except ValueError as _e:
                _lines.append(f"  not available for {_pols}: {_e}")
        else:
            _lines.append("  centrosymmetric crystal: no second-order nonlinearity")
        _lines += [
            "",
            "NOTES",
            "  Angles are internal and measured from the optic/dielectric axes; wavelengths are in vacuum.",
            "  The sign of d_eff is a convention; walk-off is included via theta + rho.",
            (f"  Reproduce:  nd.media.crystals.{xtal.value}().pmAngles_dfg({WL3:g}, {WL1:g}, {T_C:g}, deg=True)"
             if IS_DFG else
             f"  Reproduce:  nd.media.crystals.{xtal.value}().pmAngles_sfg({WL1:g}, {WL2:g}, {T_C:g}, deg=True)"),
        ]
        REPORT = "\n".join(_lines)
        _out = mo.vstack([
            mo.md("## Report\n\nEverything above as plain text — copy it into a lab notebook "
                  "or alongside the data it belongs to."),
            mo.ui.code_editor(REPORT, language="text", disabled=True, min_height=200),
        ])
    else:
        REPORT = ""
        _out = mo.md("")
    _out
    return (REPORT,)


if __name__ == "__main__":
    app.run()
