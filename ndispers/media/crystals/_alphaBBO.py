import sympy
from ndispers._baseclass import Medium, wl, phi, theta, T
from ndispers.helper import vars2

class AlphaBBO(Medium):
    """
    α-BBO (α-BaB₂O₄, barium borate) crystal

    - Point group : 3̄m  (D3d); centrosymmetric, so unlike β-BaB₂O₄ this phase
      has no second-order nonlinearity and is used as a birefringent optic
    - Crystal system : Trigonal
    - Dielectric principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Transparency range : 0.19 to 3.5 µm

    Sellmeier equation
    ------------------
        n(wl) = sqrt(A_i + B_i/(wl**2 - C_i) - D_i * wl**2)  for i = o, e
    
    Validity range
    --------------
    Not stated in the source (vendor Sellmeier equation).

    Note
    ----
    The Sellmeier equation has no temperature term: T_degC is accepted
    for signature uniformity and ignored, and dndT returns 0.

    Ref
    ---
    Sellmeier equation: CASTECH vendor data for α-BBO.
    Structure and transparency range:
    Fedorov, P. P., A. E. Kokh, and N. G. Kononova. "Barium borate β-BaB₂O₄ as a material for nonlinear optics." Russian Chemical Reviews 71.8 (2002): 651-671.
    """
    __slots__ = ["_A_o", "_B_o", "_C_o", "_D_o", 
                 "_A_e", "_B_e", "_C_e", "_D_e"]

    def __init__(self):
        super().__init__()
        self._plane = 'arb'
        self._theta_rad = 'var'
        self._phi_rad = 'arb'

        """ Constants of dispersion formula """
        # For ordinary ray
        self._A_o = 2.7471
        self._B_o = 0.01878
        self._C_o = 0.01822
        self._D_o = 0.01354
        # For extraordinary ray
        self._A_e = 2.37153
        self._B_e = 0.01224
        self._C_e = 0.01667
        self._D_e = 0.01516
    

    
    
    def n_o_expr(self):
        """ Sympy expression, dispersion formula for o-wave """
        return sympy.sqrt(self._A_o + self._B_o / (wl**2 - self._C_o) - self._D_o * wl**2)
    
    def n_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-wave """
        return sympy.sqrt(self._A_e + self._B_e / (wl**2 - self._C_e) - self._D_e * wl**2)

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
    

    

    
    
    
    
    
