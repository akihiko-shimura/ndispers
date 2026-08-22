# /// script
# requires-python = ">=3.11"
# dependencies = ["ndispers==0.9.1", "numpy", "plotly"]
# ///
"""Refractive index explorer — browse ndispers media interactively.

Runs in the browser via marimo's WASM export; see apps/README.md.
"""
import marimo

app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _():
    import marimo as mo
    return (mo,)


@app.cell(hide_code=True)
def _():
    import time

    _t0 = time.perf_counter()
    import numpy as np
    import plotly.graph_objects as go
    import sympy
    import ndispers as nd
    import ndispers.media.crystals as _C
    import ndispers.media.glasses as _G

    MEDIA = {}
    for _mod in (_C, _G):
        for _name in sorted(n for n in dir(_mod) if not n.startswith("_")):
            _obj = getattr(_mod, _name)
            if isinstance(_obj, type):
                MEDIA[_name] = _obj
    IMPORT_S = time.perf_counter() - _t0
    return IMPORT_S, MEDIA, go, nd, np, sympy, time


@app.cell(hide_code=True)
def _(np):
    # method -> (label, unit, post-processing from ndispers' native unit)
    QUANTITIES = {
        "n": ("refractive index n", "", lambda v: v),
        "ng": ("group index n_g", "", lambda v: v),
        "GVD": ("GVD", "fs²/mm", lambda v: v),
        "TOD": ("TOD", "fs³/mm", lambda v: v),
        "dndT": ("dn/dT", "1/K", lambda v: v),
        "woa_theta": ("walk-off angle", "deg", np.degrees),
    }
    # variable -> (label, unit, default sweep)
    VARIABLES = {
        "wl": ("wavelength", "nm", (250.0, 2000.0)),
        "T": ("temperature", "°C", (-50.0, 300.0)),
        "angle": ("angle", "deg", (0.0, 90.0)),
    }
    return QUANTITIES, VARIABLES


@app.cell(hide_code=True)
def _(IMPORT_S, MEDIA, mo, nd):
    mo.md(
        f"""
    # ndispers explorer

    Refractive index and its derived quantities for {len(MEDIA)} media, computed
    in your browser — no server, no install. Powered by
    [ndispers {nd.__version__}](https://github.com/akihiko-shimura/ndispers).

    <sub>import took {IMPORT_S:.2f} s</sub>
    """
    )
    return


@app.cell(hide_code=True)
def _(MEDIA, mo):
    medium = mo.ui.dropdown(
        options=list(MEDIA), value="BetaBBO_Eimerl1987", label="Medium"
    )
    pol = mo.ui.radio(options=["o", "e"], value="o", inline=True)
    T = mo.ui.number(-50, 300, value=25, step=1)
    angle = mo.ui.number(0, 90, value=0, step=0.1)
    medium
    return T, angle, medium, pol


@app.cell(hide_code=True)
def _(MEDIA, T, angle, medium, np, pol):
    x = MEDIA[medium.value]()
    IS_ISOTROPIC = len(x.symbols) == 2

    def call(meth, wl):
        """Media differ in signature: isotropic ones take no angle or pol."""
        f = getattr(x, meth)
        if IS_ISOTROPIC:
            return f(wl, T.value)
        return f(wl, np.radians(angle.value), T.value, pol=pol.value)

    return IS_ISOTROPIC, call, x


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
def _(IS_ISOTROPIC, mo, pol, sympy, x):
    _pols = ["o"] if IS_ISOTROPIC else [pol.value]
    _md = "\n\n".join(
        f"$$n_{{{_p}}} = {sympy.latex(x.n_expr(_p))}$$" for _p in _pols
    )
    mo.accordion({
        "Sellmeier equation as evaluated": mo.md(
            "What ndispers differentiates and evaluates, with the coefficients "
            "already substituted — not a quotation from the paper. **λ is in µm "
            "and T in °C in this expression.**\n\n" + _md
        )
    })
    return


@app.cell(hide_code=True)
def _(mo):
    wl0 = mo.ui.number(200.0, 4000.0, value=1064.0, step=0.1)
    return (wl0,)


@app.cell(hide_code=True)
def _(IS_ISOTROPIC, QUANTITIES, VARIABLES, mo):
    quantity = mo.ui.dropdown(
        options={v[0]: k for k, v in QUANTITIES.items()},
        value="refractive index n", label="Plot",
    )
    _vars = {v[0]: k for k, v in VARIABLES.items()
             if not (IS_ISOTROPIC and k == "angle")}
    variable = mo.ui.dropdown(options=_vars, value="wavelength", label="against")
    return quantity, variable


@app.cell(hide_code=True)
def _(IS_ISOTROPIC, T, angle, mo, quantity, variable, wl0):
    # held values in the order the methods take them: n(wl, angle, T, pol)
    def _box(label, widget, width="5.5rem"):
        return mo.vstack(
            [mo.md(f"<sub>{label}</sub>"),
             mo.Html(f'<div style="width:{width}">{widget}</div>')],
            gap=0,
        )

    _boxes = [_box("λ (nm)", wl0)]
    if not IS_ISOTROPIC:
        _boxes += [_box("θ or φ (deg)", angle), _box("T (°C)", T),
                   _box("pol", pol, "5rem")]
    else:
        _boxes += [_box("T (°C)", T)]

    mo.vstack([
        mo.hstack(_boxes, justify="start", gap=0.5),
        mo.hstack([quantity, variable], justify="start", gap=0.75),
    ], gap=0.4)
    return


@app.cell(hide_code=True)
def _(
    IS_ISOTROPIC, QUANTITIES, T, VARIABLES, angle, go, mo, np, pol,
    quantity, time, variable, wl0, x,
):
    _q, _v = quantity.value, variable.value
    _qlabel, _qunit, _post = QUANTITIES[_q]
    _vlabel, _vunit, (_lo, _hi) = VARIABLES[_v]

    _sweep = np.linspace(_lo, _hi, 400)
    if _v == "angle":                       # keep off the exact 0 and 90 poles
        _sweep = np.linspace(0.02, 89.98, 400)

    # hold the two variables that are not being swept at their control values
    _held = {"wl": wl0.value, "T": T.value, "angle": angle.value}
    _held[_v] = _sweep
    _held["wl"] = np.asarray(_held["wl"]) * 1e-3          # nm on screen, µm inside

    _t0 = time.perf_counter()
    _f = getattr(x, _q)
    if IS_ISOTROPIC:
        _y = _post(_f(_held["wl"], _held["T"]))
    else:
        _y = _post(_f(_held["wl"], np.radians(_held["angle"]), _held["T"], pol=pol.value))
    _elapsed = time.perf_counter() - _t0

    # a held variable that is not swept is fixed; say so in the subtitle
    _fixed = ", ".join(
        f"{VARIABLES[k][0]} = {v:g} {u}"
        for k, v, u in (("wl", wl0.value, "nm"), ("angle", angle.value, "deg"),
                        ("T", T.value, "°C"))
        if k != _v and not (IS_ISOTROPIC and k == "angle")
    )

    _ylab = _qlabel + (f" ({_qunit})" if _qunit else "")
    _fig = go.Figure(
        go.Scatter(
            x=_sweep, y=np.broadcast_to(_y, _sweep.shape), mode="lines",
            hovertemplate=f"{_vlabel} %{{x:.4g}} {_vunit}<br>"
                          f"{_qlabel} %{{y:.6g}} {_qunit}<extra></extra>",
        )
    )
    _fig.update_layout(
        title=dict(
            text=f"{type(x).__name__}   {_fixed}"
                 + ("" if IS_ISOTROPIC else f", pol = {pol.value}"),
            font_size=11, x=0, xanchor="left", y=0.97, yanchor="top",
        ),
        xaxis_title=f"{_vlabel} ({_vunit})",
        yaxis_title=_ylab,
        width=560,
        margin=dict(l=60, r=20, t=58, b=45),
        height=310,
        font_size=11,
        template="plotly_white",
        showlegend=False,
        modebar=dict(orientation="v"),
    )

    # a flat curve is usually physics, not a bug: say which
    _yy = np.asarray(_y)
    _flat = _yy.ndim == 0 or np.ptp(_yy) <= 1e-12 * max(abs(np.max(_yy)), 1e-30)
    if _flat and _v == "angle":
        _note = ("  ·  flat because the ordinary index of a uniaxial crystal "
                 "has no angle dependence — switch to pol = e"
                 if pol.value == "o" else
                 "  ·  this quantity does not depend on the angle here")
    elif _flat and _v == "T":
        _note = "  ·  flat because this medium's Sellmeier equation has no temperature term"
    elif IS_ISOTROPIC:
        _note = "  ·  isotropic medium: angle and polarization ignored"
    else:
        _note = ""

    mo.vstack([
        _fig,
        mo.md(f"<sub>400 points: {1e3 * _elapsed:.0f} ms{_note}</sub>"),
    ])
    return


@app.cell(hide_code=True)
def _(call, mo, time, wl0):
    _rows = []
    _t0 = time.perf_counter()
    for _label, _meth, _unit in [
        ("n", "n", ""),
        ("dn/dλ", "dn_wl", "1/nm"),
        ("group index n_g", "ng", ""),
        ("group velocity", "GV", "µm/fs"),
        ("GVD", "GVD", "fs²/mm"),
        ("TOD", "TOD", "fs³/mm"),
        ("dn/dT", "dndT", "1/K"),
        ("walk-off θ", "woa_theta", "deg"),
    ]:
        try:
            _v = float(call(_meth, wl0.value * 1e-3))   # nm -> µm
            if _meth == "woa_theta":
                _v = _v * 180 / 3.141592653589793
            elif _meth == "dn_wl":
                _v = _v * 1e-3                       # per µm -> per nm
            _rows.append(f"| {_label} | {_v:.6g} | {_unit} |")
        except Exception as _e:                       # medium may not define it
            _rows.append(f"| {_label} | — | {type(_e).__name__} |")
    _derived_s = time.perf_counter() - _t0

    mo.md(
        "| quantity | value | unit |\n|---|---|---|\n"
        + "\n".join(_rows)
        + f"\n\n<sub>all eight quantities: {1e3 * _derived_s:.0f} ms"
        " (first call per medium builds the symbolic derivatives)</sub>"
    )
    return


if __name__ == "__main__":
    app.run()
