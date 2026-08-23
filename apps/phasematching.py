# /// script
# requires-python = ">=3.11"
# dependencies = ["ndispers==0.10.0", "numpy", "scipy", "sympy", "matplotlib"]
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
    import sympy
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
    return XTALS, brentq, nd, np, plt, sympy


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
    wl1 = mo.ui.number(200.0, 4000.0, value=1064.0, step=0.1, label="λ₁ (nm)")
    wl2 = mo.ui.number(200.0, 4000.0, value=1064.0, step=0.1, label="λ₂ (nm)")
    T = mo.ui.number(-50, 300, value=25, step=1, label="Temperature (°C)")
    L = mo.ui.number(0.1, 100.0, value=10.0, step=0.1, label="Crystal length (mm)")
    return L, T, shg, wl1, wl2, xtal


@app.cell(hide_code=True)
def _(L, T, mo, shg, wl1, wl2, xtal):
    mo.hstack(
        [
            mo.vstack([xtal, shg], gap=0.25),
            mo.vstack([wl1, wl2 if not shg.value else mo.md("<sub>λ₂ = λ₁</sub>")]),
            mo.vstack([T, L], gap=0.25),
        ],
        justify="start",
        gap=0.75,
    )
    return


@app.cell(hide_code=True)
def _(XTALS, T, shg, wl1, wl2, xtal):
    x = XTALS[xtal.value]()
    # ndispers works in µm; the interface is in nm
    WL1 = wl1.value * 1e-3
    WL2 = WL1 if shg.value else wl2.value * 1e-3
    WL3 = 1.0 / (1.0 / WL1 + 1.0 / WL2)
    T_C = T.value
    ANGLE_KEY = "theta" if x.theta_rad == "var" else "phi"
    return ANGLE_KEY, T_C, WL1, WL2, WL3, x


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
def _(HAS_PM, POLS, mo, sympy, x):
    if HAS_PM:
        _md = "\n\n".join(
            f"$$n_{{{_p}}} = {sympy.latex(x.n_expr(_p))}$$"
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
        + (f"**{SOLUTIONS[k]:.3f}°** | {SOLUTIONS[k] * 3.141592653589793 / 180:.5f} rad"
           if k in SOLUTIONS else "— none | —")
        + " |"
        for k in TYPE_OF
    )
    mo.md(
        f"""
    ## Phase-matching angles

    λ₁ = {WL1 * 1e3:.2f} nm, λ₂ = {WL2 * 1e3:.2f} nm → **λ₃ = {WL3 * 1e3:.2f} nm**
    &nbsp;·&nbsp; the varying angle of this crystal is **{ANGLE_KEY}**

    | polarizations (1,2,3) | type | angle | |
    |---|---|---|---|
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
    return ANGLE_PM, HAS_PM, PLOT_TITLE, POLS


@app.cell(hide_code=True)
def _(ANGLE_PM, HAS_PM, L, PLOT_TITLE, POLS, T_C, WL1, WL2, mo, np, plt, x):
    if HAS_PM:
        _w = np.radians(0.15)
        _a = np.linspace(ANGLE_PM - _w, ANGLE_PM + _w, 400)
        _dk = x.dk_sfg(WL1, WL2, _a, T_C, *POLS)
        _sinc = x.pmFactor_sfg(WL1, WL2, _a, T_C, *POLS, L.value)

        plt.rcParams.update({"font.size": 8})
        _fig, (_a1, _a2) = plt.subplots(2, 1, sharex=True, figsize=(5.25, 3.3))
        _a1.set_title(PLOT_TITLE, fontsize=8)
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
    ang_unit = mo.ui.dropdown(options=list(ANG_UNITS), value="mdeg", label="Angle unit")
    mo.hstack([threshold, wl_unit, ang_unit], justify="start", gap=0.75)
    return ANG_UNITS, WL_UNITS, ang_unit, threshold, wl_unit


@app.cell(hide_code=True)
def _(
    ANGLE_PM, ANG_UNITS, HAS_PM, L, POLS, T_C, WL1, WL2, WL_UNITS,
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

        _rows = [
            _row("angle", _d_ang, _as, _au),
            _row("temperature", _d_T, 1, "°C"),
            _row("λ₁ only, λ₂ fixed", _d_wl1, _ws, _wu),
        ]
        if shg.value:
            _rows.append(_row("λ₁ and λ₂ tuned together", _d_wboth, _ws, _wu))

        _out = mo.md(
            f"""
    ## Acceptance bandwidths

    Full width at {_th:.2f} of sinc²(ΔkL/2) for L = {_L:g} mm, and the
    length product, which is what is usually tabulated.

    | | full width at L = {_L:g} mm | × L |
    |---|---|---|
    """ + "\n".join("    " + r for r in _rows) + """

    The wavelength row to compare against tabulated values is **λ₁ only** —
    that is what SNLO reports. Tuning λ₁ and λ₂ together makes Δk move twice
    as fast, so that width is exactly half; it is the relevant one when a
    single broadband beam supplies both photons.
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
    ANGLE_PM, DEFF_MAX, HAS_DEFF, HAS_PM, PHI_FREE, PHI_OPT, PLOT_TITLE, POLS,
    T_C, WL1, WL2, deff_args, mo, np, phi_cut, plt, x,
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
            plt.rcParams.update({"font.size": 8})
            _fig, _ax = plt.subplots(figsize=(5.25, 2.25))
            _ax.set_title(PLOT_TITLE + (f", φ = {phi_cut.value:.0f}°" if PHI_FREE else ""), fontsize=8)
            _ax.plot(np.degrees(_a), _curve)
            _ax.axvline(np.degrees(ANGLE_PM), color="r", lw=0.8, ls="--")
            if _group == "Biax_2":
                _ax.axvline(180 - np.degrees(ANGLE_PM), color="r", lw=0.8, ls=":")
            _ax.set_xlabel("angle (deg)")
            _ax.set_ylabel(r"|$d_\mathrm{eff}$| (pm/V)")
            _ax.grid(alpha=0.3)
            plt.tight_layout()

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


if __name__ == "__main__":
    app.run()
