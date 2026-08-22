import sympy

from ndispers._baseclass import Medium, T, phi, theta, wl
from ndispers.helper import vars2

class KDP(Medium):
    """
    KDP (K H_2 P O_4, Potassium Dihydrogen Phosphate) crystal

    - Point group : -42m  (D_{2d})
    - Crystal system : Tetragonal
    - Dielectic principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Tranparency range : 0.174 to 1.57 µm

    Sellmeier equation
    ------------------
    n(wl_um) = sqrt(A_i + B_i/(wl_um**2 - C_i) - D_i * wl_um**2/(wl_um**2 - 400))  for i = o, e

    Thermo-optic coefficient
    -------------------------
    dn_o/dT = 0 # to be incorporated
    dn_e/dT = 0 # to be incorporated
    
    Validity range
    ---------------
    T = 24.8 degC

    Ref
    ----
    Zernike, Frits. "Refractive indices of ammonium dihydrogen phosphate and potassium dihydrogen phosphate between 2000 Å and 1.5 μ." JOSA 54.10 (1964): 1215-1220.
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
        self._A_o = 2.259276
        self._B_o = 0.01008956
        self._C_o = 0.012942625
        self._D_o = 13.00522
        # For extraordinary ray
        self._A_e = 2.132668
        self._B_e = 0.00863749
        self._C_e = 0.012281043
        self._D_e = 3.22799
    

    
    def dndT_o_expr(self):
        return 0
    
    def dndT_e_expr(self):
        return 0
    
    def n_o_expr(self):
        """ Sympy expression, dispersion formula for o-wave """
        return sympy.sqrt(self._A_o + self._B_o / (wl**2 - self._C_o) - self._D_o / (wl**2 - 400)) + self.dndT_o_expr() * (T - 24.8)
    
    def n_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-wave """
        return sympy.sqrt(self._A_e + self._B_e / (wl**2 - self._C_e) - self._D_e / (wl**2 - 400)) + self.dndT_e_expr() * (T - 24.8)

    def n_expr(self, pol):
        """"
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
    

    

    
    
    
    
    
    
    