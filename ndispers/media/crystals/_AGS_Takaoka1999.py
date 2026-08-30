from ndispers._sym import sympy
from ndispers._baseclass import wl, phi, theta, T
from ndispers.groups import Uniax_42m
from ndispers.helper import vars2

class AGS(Uniax_42m):
    """
    AGS (AgGaS₂, silver thiogallate) crystal, Takaoka & Kato 1999 parameterisation

    DEPRECATED. Use AGS_Kato1996, which reproduces the phase-matching table of
    Kato et al. 2019 to 0.3 degrees where this set deviates by up to 1.35. The
    thermo-optic dispersion formula that is the actual contribution of Takaoka
    & Kato 1999 is not implemented here (the paper was not available for
    transcription), so this class offers nothing the other does not. It is kept
    for one release as an independent check on the transcription of the
    room-temperature set, and will be removed in 1.0.

    - Point group : 4̄2m  (D2d)
    - Crystal system : Tetragonal
    - Dielectric principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Transparency range : about 0.47 to 13 µm

    Sellmeier equation
    ------------------
        n(wl)**2 = A_i + B_i / (wl**2 - C_i) + D_i * wl**2   for i = o, e
        (C in µm**2)

    Validity range
    ---------------
    0.58 to 10.59 µm, at 20 degC

    Note
    ----
    The Sellmeier equation has no temperature term: T_degC is accepted for
    signature uniformity and ignored, and dndT returns 0.
    The source also gives a thermo-optic dispersion formula, which is not
    implemented here (the paper was not available for transcription; only the
    room-temperature Sellmeier set as tabulated by refractiveindex.info,
    main/AgGaS2/nk/Takaoka-o, -e, accessed 2026-08-23). Cross-checked against
    the Kato & Shirahata 1996 set (AGS_Kato1996) to 1e-3 over 1-10 µm.

    Ref
    ---
    Sellmeier equation:
      Takaoka, E., & Kato, K. (1999). Thermo-optic dispersion formula for AgGaS2. Applied Optics, 38(21), 4577-4580. https://doi.org/10.1364/ao.38.004577
    Nonlinear optical coefficient:
      Roberts, D. A. (1992). Simplified characterization of uniaxial and biaxial nonlinear optical crystals: a plea for standardization of nomenclature and conventions. IEEE Journal of Quantum Electronics, 28(10), 2057-2074. https://doi.org/10.1109/3.159516
    """
    _wl_range = (0.58, 10.59)  # um, Sellmeier validity (see docstring)

    __slots__ = ["_A_o", "_B_o", "_C_o", "_D_o",
                 "_A_e", "_B_e", "_C_e", "_D_e"]

    _d_ref = {"d36": (17.5, 1.064, 1.064)}
    _d_note = ("Roberts 1992, Table V: 17.5 pm/V for 1.064 um SHG (relative to quartz), and "
               "11.2 pm/V for 10.6 um SHG, on his scale d36(KDP) = 0.39 pm/V. The 1.064 um "
               "value is held; Miller scaling from it gives 11.1 pm/V for 10.6 um SHG, against "
               "his 11.2 pm/V entry.")

    def __init__(self):
        super().__init__()
        self._plane = 'arb'
        self._theta_rad = 'var'
        self._phi_rad = 'arb'

        """ Constants of dispersion formula """
        # For ordinary ray
        self._A_o = 5.7975
        self._B_o = 0.2311
        self._C_o = 0.0688
        self._D_o = -0.00257
        # For extraordinary ray
        self._A_e = 5.5436
        self._B_e = 0.2230
        self._C_e = 0.0946
        self._D_e = -0.00261

    def n_o_expr(self):
        """ Sympy expression, dispersion formula for o-wave """
        return sympy.sqrt(self._A_o + self._B_o / (wl**2 - self._C_o) + self._D_o * wl**2)

    def n_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-wave """
        return sympy.sqrt(self._A_e + self._B_e / (wl**2 - self._C_e) + self._D_e * wl**2)

    def n_expr(self, pol):
        """ Sympy expression, dispersion formula of a general ray with an angle theta to optic axis. If theta = 0, this expression reduces to 'no_expre'. """
        if pol == 'o':
            return self.n_o_expr()
        elif pol == 'e':
            return self.n_e_expr() / sympy.sqrt( sympy.sin(theta)**2 + (self.n_e_expr()/self.n_o_expr())**2 * sympy.cos(theta)**2 )
        else:
            raise ValueError("pol = '%s' must be 'o' or 'e'" % pol)
