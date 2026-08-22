import sympy
from ndispers._baseclass import Medium, wl, phi, theta, T
from ndispers.helper import vars2

class Calcite(Medium):
    """
    calcite (CaCO₃) crystal

    - Point group : 3̄m  (D3d)
    - Crystal system : Trigonal
    - Dielectric principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Transparency range : 

    Sellmeier equation
    ------------------
        n(wl) 
        = sqrt(1 + A1_o * wl**2 / (wl**2 - B1_o**2) + A2_o * wl**2 / (wl**2 - B2_o**2) + A3_o * wl**2 / (wl**2 - B3_o**2) + A4_o * wl**2 / (wl**2 - B4_o**2)) for o-wave
        = sqrt(1 + A1_e * wl**2 / (wl**2 - B1_e**2) + A2_e * wl**2 / (wl**2 - B2_e**2) + A3_e * wl**2 / (wl**2 - B3_e**2)) for e-wave
    
    Validity range
    ---------------
    0.2 to 2.2 µm for o-wave
    0.2 to 3.3 µm for e-wave

    Note
    ----
    The Sellmeier equation has no temperature term: T_degC is accepted
    for signature uniformity and ignored, and dndT returns 0.

    Ref
    ---
    Tropf, W. J., Thomas, M. E., & Harris, T. J. (1995). Properties of crystals and glasses. In M. Bass (Ed.), Handbook of Optics, Vol. II, Chapter 33. McGraw-Hill.
    """
    __slots__ = ["_A1_o", "_B1_o", "_A2_o", "_B2_o", "_A3_o", "_B3_o", "_A4_o", "_B4_o",
                 "_A1_e", "_B1_e", "_A2_e", "_B2_e", "_A3_e", "_B3_e"]

    def __init__(self):
        super().__init__()
        self._plane = 'arb'
        self._theta_rad = 'var'
        self._phi_rad = 'arb'

        """ Constants of dispersion formula """
        # For ordinary ray
        self._A1_o = 0.8559
        self._B1_o = 0.0588
        self._A2_o = 0.8391
        self._B2_o = 0.141
        self._A3_o = 0.0009
        self._B3_o = 0.197
        self._A4_o = 0.6845
        self._B4_o = 7.005
        # For extraordinary ray
        self._A1_e = 1.0856
        self._B1_e = 0.07897
        self._A2_e = 0.0988
        self._B2_e = 0.142
        self._A3_e = 0.317
        self._B3_e = 11.468
    

    
    
    def n_o_expr(self):
        """ Sympy expression, dispersion formula for o-wave """
        return sympy.sqrt(1 + self._A1_o * wl**2 / (wl**2 - self._B1_o**2) + self._A2_o * wl**2 / (wl**2 - self._B2_o**2) + self._A3_o * wl**2 / (wl**2 - self._B3_o**2) + self._A4_o * wl**2 / (wl**2 - self._B4_o**2))
    
    def n_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-wave """
        return sympy.sqrt(1 + self._A1_e * wl**2 / (wl**2 - self._B1_e**2) + self._A2_e * wl**2 / (wl**2 - self._B2_e**2) + self._A3_e * wl**2 / (wl**2 - self._B3_e**2))

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
    

    

    
    
    
    
    
