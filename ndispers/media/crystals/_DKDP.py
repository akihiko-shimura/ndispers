from ndispers._sym import sympy

from ndispers._baseclass import T, phi, theta, wl
from ndispers.groups import Uniax_42m
from ndispers.helper import vars2

class DKDP(Uniax_42m):
    """
    DKDP (KD₂PO₄, deuterated potassium dihydrogen phosphate, KD*P) crystal

    - Point group : 4̄2m  (D2d)
    - Crystal system : Tetragonal
    - Dielectric principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Transparency range : 0.2 to 2.1 µm (deuteration moves the infrared edge beyond KDP's 1.57 µm; the source states no figure)

    Sellmeier equation
    ------------------
        n(wl)**2 = A_i + B_i/(1 - C_i/wl**2) + D_i/(1 - E/wl**2)  for i = o, e
        E = 36 µm**2

    Thermo-optic coefficient
    -------------------------
        2 * n_i * dn_i/dT = G_i * Re_i + H_i * Re_i**2 + L_i * Rl_i**2   for i = o, e
        Re_i = wl**2/(wl**2 - wlg_i**2),  Rl_i = wl**2/(wl**2 - wll_i**2)

    The temperature correction dn/dT * (T - 25) vanishes at 25 degC, the
    temperature of the room-temperature KD*P index data the Sellmeier
    equation was fitted to (Ghosh 1992, Sec. 2).

    Validity range
    ---------------
    0.2 to 1.5 µm, the range of the KDP-family index data the fit rests on.
    Ghosh 1992 states no separate range for KD*P; its thermo-optic model is
    plotted over 0.1 to 2.1 µm, but the Sellmeier poles at 0.114 µm (o) and
    0.120 µm (e) make n complex below about 0.13 µm. Sellmeier rms deviation
    7.7e-5 (o), 3.0e-5 (e) in the source's Table I.

    Ref
    ---
    Sellmeier equation and thermo-optic coefficients:
      Ghosh, G. (1992). Dispersion of thermo-optic coefficients and temperature-dependent nonlinear optical devices of some nonlinear crystals. Proc. SPIE, 1622, 49-53. https://doi.org/10.1117/12.637003

    Nonlinear optical coefficient:
      Roberts, D. A. (1992). Simplified characterization of uniaxial and biaxial nonlinear optical crystals: a plea for standardization of nomenclature and conventions. IEEE Journal of Quantum Electronics, 28(10), 2057-2074. https://doi.org/10.1109/3.159516
    """
    __slots__ = ["_A_o", "_B_o", "_C_o", "_D_o",
                 "_A_e", "_B_e", "_C_e", "_D_e", "_E",
                 "_G_o", "_H_o", "_L_o", "_wlg_o", "_wll_o",
                 "_G_e", "_H_e", "_L_e", "_wlg_e", "_wll_e"]

    _d_ref = {"d36": (0.37, 1.064, 1.064)}
    _d_note = ("Roberts 1992, Table V: 0.37 pm/V at 1.064 um SHG, from the absolute and "
               "KDP-relative measurements he compiles (0.36-0.38). Deuteration level affects "
               "the value slightly. Miller scaling untested for KD*P specifically.")

    def __init__(self):
        super().__init__()
        self._plane = 'arb'
        self._theta_rad = 'var'
        self._phi_rad = 'arb'

        """ Constants of dispersion formula """
        # Sellmeier coefficients, Ghosh 1992 Table I
        # For ordinary ray
        self._A_o = 1.4513478
        self._B_o = 0.7905108
        self._C_o = 0.0129717
        self._D_o = 0.5364005
        # For extraordinary ray
        self._A_e = 1.5596463
        self._B_e = 0.5697706
        self._C_e = 0.0144792
        self._D_e = 0.3170391
        # lattice-absorption term, shared by both polarizations
        self._E = 36.0  #µm^2
        # Thermo-optic coefficients, Ghosh 1992 Table II (tabulated in 1e-5/degC)
        self._G_o = -7.357e-5
        self._H_o = 28.990e-5
        self._L_o = -29.426e-5
        self._wlg_o = 0.132  #µm
        self._wll_o = 0.133  #µm
        self._G_e = -13.462e-5
        self._H_e = 28.156e-5
        self._L_e = -19.797e-5
        self._wlg_e = 0.128  #µm
        self._wll_e = 0.134  #µm

    def _n_T25_o_expr(self):
        """ Sympy expression, dispersion formula for o-wave at 25 degC """
        return sympy.sqrt(self._A_o + self._B_o/(1 - self._C_o/wl**2) + self._D_o/(1 - self._E/wl**2))

    def _n_T25_e_expr(self):
        """ Sympy expression, dispersion formula for e-wave at 25 degC """
        return sympy.sqrt(self._A_e + self._B_e/(1 - self._C_e/wl**2) + self._D_e/(1 - self._E/wl**2))

    def dndT_o_expr(self):
        """ Sympy expression for thermo-optic coefficient of o-wave (dn/dT) """
        Re = wl**2 / (wl**2 - self._wlg_o**2)
        Rl = wl**2 / (wl**2 - self._wll_o**2)
        return (self._G_o*Re + self._H_o*Re**2 + self._L_o*Rl**2) / (2 * self._n_T25_o_expr())

    def dndT_e_expr(self):
        """ Sympy expression for thermo-optic coefficient of e-wave (dn/dT) """
        Re = wl**2 / (wl**2 - self._wlg_e**2)
        Rl = wl**2 / (wl**2 - self._wll_e**2)
        return (self._G_e*Re + self._H_e*Re**2 + self._L_e*Rl**2) / (2 * self._n_T25_e_expr())

    def n_o_expr(self):
        """ Sympy expression, dispersion formula for o-wave """
        return self._n_T25_o_expr() + self.dndT_o_expr() * (T - 25)

    def n_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-wave """
        return self._n_T25_e_expr() + self.dndT_e_expr() * (T - 25)

    def n_expr(self, pol):
        """
        Sympy expression,
        dispersion formula of a general ray with an angle theta to optic axis. If theta = 0, this expression reduces to 'no_expre'.

        n(theta) = n_e / sqrt( sin(theta)**2 + (n_e/n_o)**2 * cos(theta)**2 )

        """
        if pol == 'o':
            return self.n_o_expr()
        elif pol == 'e':
            return self.n_e_expr() / sympy.sqrt( sympy.sin(theta)**2 + (self.n_e_expr()/self.n_o_expr())**2 * sympy.cos(theta)**2 )
        else:
            raise ValueError("pol = '%s' must be 'o' or 'e'" % pol)
