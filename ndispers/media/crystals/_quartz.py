from ndispers._sym import sympy
from ndispers._baseclass import wl, phi, theta, T
from ndispers.groups import Uniax_32
from ndispers.helper import vars2

class Quartz(Uniax_32):
    """
    α-quartz (crystalline SiO₂) crystal

    - Point group : 32  (D3)
    - Crystal system : Trigonal
    - Dielectric principal axis, z // c-axis (x, y-axes are arbitrary)
    - Positive uniaxial, with optic axis parallel to z-axis
    - Transparency range : 0.155 to 4.0 µm

    Sellmeier equation
    ------------------
        n(wl)**2 = 1 + sum_k A_k_i * wl**2 / (wl**2 - B_k_i**2),  k = 1..5,  for i = o, e

    Validity range
    ---------------
    0.18 to 0.71 µm (the fit's stated range; the two infrared poles carry it
    smoothly to the infrared edge, but no data constrain it there)

    Note
    ----
    The Sellmeier equation has no temperature term: T_degC is accepted for
    signature uniformity and ignored, and dndT returns 0. Tropf et al. list
    dn/dT = -6.2e-6 /K (o) and -7.0e-6 /K (e) at 0.546 µm (Table 20 of the
    reference). Quartz is optically active; the rotation of the plane of
    polarization along the optic axis is not modelled, as it does not enter
    phase matching away from the axis.

    α-quartz is the absolute-scale reference of second-order nonlinear
    optics: d11 = 0.30 pm/V at 1.064 µm SHG is the value relative to which
    many other crystals' coefficients have been measured.

    Ref
    ---
    Sellmeier equation (as tabulated in):
      Tropf, W. J., Thomas, M. E., & Harris, T. J. (1995). Properties of crystals and glasses. In M. Bass (Ed.), Handbook of Optics, Vol. II, Chapter 33, Table 22. McGraw-Hill.
    Underlying measurements:
      Radhakrishnan, T. (1951). Further studies on the temperature variation of the refractive index of crystals. Proceedings of the Indian Academy of Sciences A, 33, 22-34 (not registered with Crossref).
    Nonlinear optical coefficient:
      Shoji, I., Kondo, T., Kitamoto, A., Shirane, M., & Ito, R. (1997). Absolute scale of second-order nonlinear-optical coefficients. JOSA B, 14(9), 2268-2294. https://doi.org/10.1364/josab.14.002268
    """
    __slots__ = ["_A1_o", "_B1_o", "_A2_o", "_B2_o", "_A3_o", "_B3_o", "_A4_o", "_B4_o", "_A5_o", "_B5_o",
                 "_A1_e", "_B1_e", "_A2_e", "_B2_e", "_A3_e", "_B3_e", "_A4_e", "_B4_e", "_A5_e", "_B5_e"]

    _d_ref = {"d11": (0.30, 1.064, 1.064)}
    _d_note = ("Shoji et al. 1997, absolute, 1.064 um SHG (d11 of alpha-quartz is the "
               "reference of their absolute scale and of Roberts 1992). Alford & Smith "
               "2001 did not test Miller scaling for quartz.")

    def __init__(self):
        super().__init__()
        self._plane = 'arb'
        self._theta_rad = 'var'
        self._phi_rad = 'arb'

        """ Constants of dispersion formula """
        # For ordinary ray
        self._A1_o = 0.663044
        self._B1_o = 0.060
        self._A2_o = 0.517852
        self._B2_o = 0.106
        self._A3_o = 0.175912
        self._B3_o = 0.119
        self._A4_o = 0.565380
        self._B4_o = 8.844
        self._A5_o = 1.675299
        self._B5_o = 20.742
        # For extraordinary ray
        self._A1_e = 0.665721
        self._B1_e = 0.060
        self._A2_e = 0.503511
        self._B2_e = 0.106
        self._A3_e = 0.214792
        self._B3_e = 0.119
        self._A4_e = 0.539173
        self._B4_e = 8.792
        self._A5_e = 1.807613
        self._B5_e = 197.70

    def n_o_expr(self):
        """ Sympy expression, dispersion formula for o-wave """
        return sympy.sqrt(1 + self._A1_o * wl**2 / (wl**2 - self._B1_o**2) + self._A2_o * wl**2 / (wl**2 - self._B2_o**2) + self._A3_o * wl**2 / (wl**2 - self._B3_o**2) + self._A4_o * wl**2 / (wl**2 - self._B4_o**2) + self._A5_o * wl**2 / (wl**2 - self._B5_o**2))

    def n_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-wave """
        return sympy.sqrt(1 + self._A1_e * wl**2 / (wl**2 - self._B1_e**2) + self._A2_e * wl**2 / (wl**2 - self._B2_e**2) + self._A3_e * wl**2 / (wl**2 - self._B3_e**2) + self._A4_e * wl**2 / (wl**2 - self._B4_e**2) + self._A5_e * wl**2 / (wl**2 - self._B5_e**2))

    def n_expr(self, pol):
        """ Sympy expression, dispersion formula of a general ray with an angle theta to optic axis. If theta = 0, this expression reduces to 'no_expre'. """
        if pol == 'o':
            return self.n_o_expr()
        elif pol == 'e':
            return self.n_e_expr() / sympy.sqrt( sympy.sin(theta)**2 + (self.n_e_expr()/self.n_o_expr())**2 * sympy.cos(theta)**2 )
        else:
            raise ValueError("pol = '%s' must be 'o' or 'e'" % pol)
