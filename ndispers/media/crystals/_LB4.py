import sympy
from ndispers._baseclass import Medium, wl, phi, theta, T
from ndispers.helper import vars2

class LB4(Medium):
    """
    LB4 or LTB (Li_2 B_4 O_7, lithium tetraborate) crystal

    - Point group : 4mm
    - Crystal system : Tetragonal
    - Dielectic principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Tranparency range : 0.16 to 3.6 um

    Dispersion formula for refractive index
    ---------------------------------------
    n(wl) = sqrt(A_i + B_i/(wl**2 - C_i) - D_i * wl**2) for i = o,e

    Thermo-optic coefficient
    -------------------------
    dn/dT = At_i + Bt_i * wl + Ct_i * wl**2 + Dt_i * wl**3 + Et_i * (T - 25)
    for i = o, e
    
    Validity range
    ---------------
    0.18 to 2.3 um
    at 101325 Pa
    dn/dT : 0.436 to 0.644 um, -40 to 60 degC

    Ref
    ----
    T. Sugawara, R. Komatsu, and S. Uda. "Linear and nonlinear optical properties of lithium tetraborate." Solid state communications 107.5 (1998): 233-237.

    
    """
    __slots__ = ["_LB4__plane", "_LB4__theta_rad", "_LB4__phi_rad",
                 "_A_o", "_B_o", "_C_o", "_D_o", 
                 "_A_e", "_B_e", "_C_e", "_D_e",
                 "_At_o", "_Bt_o", "_Ct_o", "_Dt_o", "_Et_o",
                 "_At_e", "_Bt_e", "_Ct_e", "_Dt_e", "_Et_e"]

    def __init__(self):
        super().__init__()
        self._LB4__plane = 'arb'
        self._LB4__theta_rad = 'var'
        self._LB4__phi_rad = 'arb'

        """ Constants of dispersion formula """
        # For ordinary ray
        self._A_o = 2.564310
        self._B_o = 0.012337
        self._C_o = 0.114467**2
        self._D_o = 0.019075
        # For extraordinary ray
        self._A_e = 2.386510
        self._B_e = 0.010664
        self._C_e = 0.113483**2
        self._D_e = 0.012813
        # dn/dT coefficients
        # for o-ray
        self._At_o = 1.893e-5 #1/K
        self._Bt_o = -88.17e-6 #1/(um*K)
        self._Ct_o = 1.497e-4 #1/(um^2*K)
        self._Dt_o = -8.643e-5 #1/(um^3*K)
        self._Et_o = -2.55e-8 #1/(K^2)
        # for e-ray
        self._At_e = 1.297e-5 #1/K
        self._Bt_e = -45.50e-6 #1/(um*K)
        self._Ct_e = 0.714e-4 #1/(um^2*K)
        self._Dt_e = -3.868e-5 #1/(um^3*K)
        self._Et_e = -2.08e-8 #1/(K^2)

    
    @property
    def plane(self):
        return self._LB4__plane

    @property
    def theta_rad(self):
        return self._LB4__theta_rad

    @property
    def phi_rad(self):
        return self._LB4__phi_rad
    
    @property
    def symbols(self):
        return [wl, theta, phi, T]
    
    @property
    def constants(self):
        print(vars2(self))
    
    def dndT_o_expr(self):
        return self._At_o + self._Bt_o * wl + self._Ct_o * wl**2 + self._Dt_o * wl**3 + self._Et_o * (T - 25)
    
    def dndT_e_expr(self):
        return self._At_e + self._Bt_e * wl + self._Ct_e * wl**2 + self._Dt_e * wl**3 + self._Et_e * (T - 25)
    
    def n_o_expr(self):
        """ Sympy expression, dispersion formula for o-ray """
        return sympy.sqrt(self._A_o + self._B_o / (wl**2 - self._C_o) - self._D_o * wl**2) + self.dndT_o_expr() * (T - 25)
    
    def n_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-ray """
        return sympy.sqrt(self._A_e + self._B_e / (wl**2 - self._C_e) - self._D_e * wl**2) + self.dndT_e_expr() * (T - 25)

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
    
    def n(self, wl_um, theta_rad, T_degC, pol='o'):
        """
        Refractive index as a function of wavelength, theta and phi angles for each eigen polarization of light.

        input
        ------
        wl_um     :  float or array_like, wavelength in um
        theta_rad :  float or array_like, 0 to pi radians
        T_degC    :  float or array_like, temperature of crystal in degree C.
        pol       :  {'o', 'e'}, optional, polarization of light

        return
        -------
        Refractive index, float or array_like
        """
        return super().n(wl_um, theta_rad, 0, T_degC, pol=pol)

    def dn_wl(self, wl_um, theta_rad, T_degC, pol='o'):
        return super().dn_wl(wl_um, theta_rad, 0, T_degC, pol=pol)
    
    def d2n_wl(self, wl_um, theta_rad, T_degC, pol='o'):
        return super().d2n_wl(wl_um, theta_rad, 0, T_degC, pol=pol)

    def d3n_wl(self, wl_um, theta_rad, T_degC, pol='o'):
        return super().d3n_wl(wl_um, theta_rad, 0, T_degC, pol=pol)
    
    def GD(self, wl_um, theta_rad, T_degC, pol='o'):
        """Group Delay [fs/mm]"""
        return super().GD(wl_um, theta_rad, 0, T_degC, pol=pol)
    
    def GV(self, wl_um, theta_rad, T_degC, pol='o'):
        """Group Velocity [um/fs]"""
        return super().GV(wl_um, theta_rad, 0, T_degC, pol=pol)
    
    def ng(self, wl_um, theta_rad, T_degC, pol='o'):
        """Group index, c/Group velocity"""
        return super().ng(wl_um, theta_rad, 0, T_degC, pol=pol)
    
    def GVD(self, wl_um, theta_rad, T_degC, pol='o'):
        """Group Delay Dispersion [fs^2/mm]"""
        return super().GVD(wl_um, theta_rad, 0, T_degC, pol=pol)
    
    def TOD(self, wl_um, theta_rad, T_degC, pol='o'):
        """Third Order Dispersion [fs^3/mm]"""
        return super().TOD(wl_um, theta_rad, 0, T_degC, pol=pol)
    
    def woa_theta(self, wl_um, theta_rad, T_degC, pol='e'):
        return super().woa_theta(wl_um, theta_rad, 0, T_degC, pol=pol)
    
    def woa_phi(self, wl_um, theta_rad, T_degC, pol='e'):
        return super().woa_phi(wl_um, theta_rad, 0, T_degC, pol=pol)
    
    def dndT(self, wl_um, theta_rad, T_degC, pol='o'):
        return super().dndT(wl_um, theta_rad, 0, T_degC, pol=pol)