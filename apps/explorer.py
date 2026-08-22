# /// script
# requires-python = ">=3.11"
# dependencies = ["ndispers==0.9.1", "numpy", "matplotlib"]
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
    import matplotlib.pyplot as plt
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
    return IMPORT_S, MEDIA, nd, np, plt, sympy, time


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
    pol = mo.ui.radio(options=["o", "e"], value="o", label="Polarization", inline=True)
    T = mo.ui.slider(-50, 300, value=25, step=5, label="Temperature (°C)", show_value=True)
    angle = mo.ui.slider(
        0, 90, value=0, step=1, label="Angle θ or φ (deg)", show_value=True
    )
    mo.hstack([mo.vstack([medium, pol], gap=0.25),
              mo.vstack([T, angle], gap=0.25)], justify="start", gap=0.75)
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

    def title():
        """The state a plotted curve was computed at, for its title."""
        if IS_ISOTROPIC:
            return f"{medium.value}   T = {T.value} °C"
        _ang = "θ" if x.theta_rad == "var" else "φ"
        return (f"{medium.value}   T = {T.value} °C, "
                f"{_ang} = {angle.value}°, pol = {pol.value}")

    return IS_ISOTROPIC, call, title, x


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
def _(IS_ISOTROPIC, call, mo, np, plt, time, title):
    _t0 = time.perf_counter()
    _wl = np.linspace(0.25, 2.0, 400)
    _n = call("n", _wl)
    _curve_s = time.perf_counter() - _t0

    plt.rcParams.update({"font.size": 8})
    _fig, _ax = plt.subplots(figsize=(5.25, 2.7))
    _ax.plot(_wl, _n)
    _ax.set_xlabel("wavelength (µm)")
    _ax.set_ylabel("refractive index")
    _ax.set_title(title(), fontsize=8)
    _ax.grid(alpha=0.3)
    plt.tight_layout()

    mo.vstack(
        [
            _fig,
            mo.md(
                f"<sub>curve of 400 points: {1e3 * _curve_s:.0f} ms"
                + ("  ·  isotropic medium: angle and polarization ignored" if IS_ISOTROPIC else "")
                + "</sub>"
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo):
    wl0 = mo.ui.number(0.2, 4.0, value=1.064, step=0.001, label="Wavelength (µm)")
    wl0
    return (wl0,)


@app.cell(hide_code=True)
def _(call, mo, time, wl0):
    _rows = []
    _t0 = time.perf_counter()
    for _label, _meth, _unit in [
        ("n", "n", ""),
        ("dn/dλ", "dn_wl", "1/µm"),
        ("group index n_g", "ng", ""),
        ("group velocity", "GV", "µm/fs"),
        ("GVD", "GVD", "fs²/mm"),
        ("TOD", "TOD", "fs³/mm"),
        ("dn/dT", "dndT", "1/K"),
        ("walk-off θ", "woa_theta", "deg"),
    ]:
        try:
            _v = float(call(_meth, wl0.value))
            if _meth == "woa_theta":
                _v = _v * 180 / 3.141592653589793
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
