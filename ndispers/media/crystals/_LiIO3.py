from ndispers._sym import sympy
from ndispers._baseclass import wl, phi, theta, T
from ndispers.groups import Uniax_6
from ndispers.helper import vars2

class LiIO3(Uniax_6):
    """
    LiIO₃ (lithium iodate) crystal

    - Point group : 6  (C6)
    - Crystal system : Hexagonal
    - Dielectric principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Transparency range : 0.38 to 5.5 µm (Tropf et al., Table 20)

    Sellmeier equation
    ------------------
        n(wl)**2 = A_i + B_i * wl**2 / (wl**2 - C_i) + D_i * wl**2 / (wl**2 - E_i)   for i = o, e
        (C, E in µm**2)

    Validity range
    ---------------
    0.5 to 5 µm

    Note
    ----
    The Sellmeier equation has no temperature term: T_degC is accepted for
    signature uniformity and ignored, and dndT returns 0. The set is the one
    Choy & Byer 1976 attribute to R. L. Herbst (private communication), as
    tabulated by Tropf et al.; refractiveindex.info's transcription of the same
    set reads the o-ray UV pole as 0.0350832 µm² where the Handbook text reads
    0.0350823 (a difference of 3e-7 in n). The Umegaki et al. 1971 set agrees
    with it to 1.2e-3 from 0.5 to 2 µm. Large d31 with a large walk-off; a
    classic crystal for SHG of Nd lasers and of ultrashort pulses.

    Ref
    ---
    Sellmeier equation (as tabulated in):
      Tropf, W. J., Thomas, M. E., & Harris, T. J. (1995). Properties of crystals and glasses. In M. Bass (Ed.), Handbook of Optics, Vol. II, Chapter 33, Table 22. McGraw-Hill.
    Underlying measurements and nonlinear coefficients:
      Choy, M. M., & Byer, R. L. (1976). Accurate second-order susceptibility measurements of visible and infrared nonlinear crystals. Physical Review B, 14(4), 1693-1706. https://doi.org/10.1103/physrevb.14.1693
      Roberts, D. A. (1992). Simplified characterization of uniaxial and biaxial nonlinear optical crystals: a plea for standardization of nomenclature and conventions. IEEE Journal of Quantum Electronics, 28(10), 2057-2074. https://doi.org/10.1109/3.159516
    """
    _wl_range = (0.5, 5)  # um, Sellmeier validity (see docstring)

    __slots__ = ["_A_o", "_B_o", "_C_o", "_D_o", "_E_o",
                 "_A_e", "_B_e", "_C_e", "_D_e", "_E_e"]

    _d_ref = {"d31": (4.4, 1.064, 1.064),
              "d33": (4.5, 1.064, 1.064)}
    _d_note = ("Roberts 1992, Table V: d31 = 4.4 pm/V at 1.064 um SHG (his selected value; "
               "Eckardt 1990 absolute 4.2, KDP-relative 4.4), d33 = 4.5 pm/V (d33/d31 = 1.00 at 1.32 um, "
               "Choy & Byer 1976); on his scale d36(KDP) = 0.39 pm/V. Alford & Smith 2001 "
               "find Miller scaling good for LiIO3.")

    def __init__(self):
        super().__init__()
        self._plane = 'arb'
        self._theta_rad = 'var'
        self._phi_rad = 'arb'

        """ Constants of dispersion formula """
        # For ordinary ray
        self._A_o = 2.03132
        self._B_o = 1.37623
        self._C_o = 0.0350823
        self._D_o = 1.06745
        self._E_o = 169.0
        # For extraordinary ray
        self._A_e = 1.83086
        self._B_e = 1.08807
        self._C_e = 0.0313810
        self._D_e = 0.554582
        self._E_e = 158.76

    def n_o_expr(self):
        """ Sympy expression, dispersion formula for o-wave """
        return sympy.sqrt(self._A_o + self._B_o * wl**2 / (wl**2 - self._C_o) + self._D_o * wl**2 / (wl**2 - self._E_o))

    def n_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-wave """
        return sympy.sqrt(self._A_e + self._B_e * wl**2 / (wl**2 - self._C_e) + self._D_e * wl**2 / (wl**2 - self._E_e))

    def n_expr(self, pol):
        """ Sympy expression, dispersion formula of a general ray with an angle theta to optic axis. If theta = 0, this expression reduces to 'no_expre'. """
        if pol == 'o':
            return self.n_o_expr()
        elif pol == 'e':
            return self.n_e_expr() / sympy.sqrt( sympy.sin(theta)**2 + (self.n_e_expr()/self.n_o_expr())**2 * sympy.cos(theta)**2 )
        else:
            raise ValueError("pol = '%s' must be 'o' or 'e'" % pol)
