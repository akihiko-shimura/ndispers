import sympy
from sympy.utilities import lambdify
from ndispers._baseclass import Medium, wl, phi, theta, pi

class KDP(Medium):
    """
    KDP (K H_2 P O_4, Potassium Dihydrogen Phosphate) crystal

    - Point group : 42m
    - Crystal system : Tetragonal
    - Dielectic principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Tranparency range : 0.174 - 1.57 um

    Dispersion formula of refractive index
    ---------------------------------------
    n(wl_um) = sqrt(A_i + B_i/(wl_um**2 - C_i) - D_i * wl_um**2/(wl_um**2 - 400))  for i = o, e
    
    Validity range
    ---------------
    

    Ref
    ----
    Zernike, Frits. "Refractive indices of ammonium dihydrogen phosphate and potassium dihydrogen phosphate between 2000 Å and 1.5 μ." JOSA 54.10 (1964): 1215-1220

    Usage
    ------
    >>> kdp = ndispers.media.crystals.KDP()
    >>> kdp.n(0.6, 0, pol='o') # for o-ray, n does not depend on theta.
    >>> kdp.n(0.6, 0.5*pi, pol='e') # along z-axis, it is pure e-ray.
    >>> kdp.n(0.6, 0.23*pi, pol='e')
    >>> kdp.n(0.6, 0*pi, pol='e') # for theta = 0 rad, it corresponds to o-ray.
    >>> kdp.GVD(0.6, 0.23*pi, pol='e')

    @author: Akihiko Shimura
    """
    __slots__ = ["_KDP__plane", "_KDP__theta_rad", "_KDP__phi_rad",
                 "_A_o", "_B_o", "_C_o", "_D_o", 
                 "_A_e", "_B_e", "_C_e", "_D_e"]
                 
    def __init__(self):
        super().__init__()
        self._KDP__plane = 'arb'
        self._KDP__theta_rad = 'var'
        self._KDP__phi_rad = 'arb'

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
    
    @property
    def plane(self):
        return self._KDP__plane
    @property
    def theta_rad(self):
        return self._KDP__theta_rad
    @property
    def phi_rad(self):
        return self._KDP__phi_rad
    @property
    def angles(self):
        msg =  ["plane = %s" % self._KDP__plane]
        msg += ["theta_rad = %s" % self._KDP__theta_rad]
        msg += ["phi_rad = %s" % self._KDP__phi_rad]
        print("\n".join(msg))
    @property
    def symbols(self):
        return [wl, theta, phi]

    @property
    def constants(self):
        msg  = ["A_o = %g" % self._A_o]
        msg += ["B_o = %g" % self._B_o]
        msg += ["C_o = %g" % self._C_o]
        msg += ["D_o = %g" % self._D_o]
        msg += ["A_e = %g" % self._A_e]
        msg += ["B_e = %g" % self._B_e]
        msg += ["C_e = %g" % self._C_e]
        msg += ["D_e = %g" % self._D_e]
        print("\n".join(msg))
    
    def n_o_expr(self):
        """ Sympy expression, dispersion formula for o-ray """
        return sympy.sqrt(self._A_o + self._B_o / (wl**2 - self._C_o) - self._D_o / (wl**2 - 400))
    
    def n_e_expr(self):
        """ Sympy expression, dispersion formula for theta=90 deg e-ray """
        return sympy.sqrt(self._A_e + self._B_e / (wl**2 - self._C_e) - self._D_e / (wl**2 - 400))

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
    
    def n(self, wl_um, theta_rad, pol='o'):
        """
        Refractive index as a function of wavelength, theta and phi angles for each eigen polarization of light.

        input
        ------
        wl_um     :  float, wavelength in um
        theta_rad :  float, 0 to pi radians
        pol       :  str, 'o' or 'e', polarization of light

        return
        -------
        Refractive index, float
        """
        return super().n(wl_um, theta_rad, 0.0, pol=pol)

    def dn_wl(self, wl_um, theta_rad, pol='o'):
        return super().dn_wl(wl_um, theta_rad, 0.0, pol=pol)
    
    def d2n_wl(self, wl_um, theta_rad, pol='o'):
        return super().d2n_wl(wl_um, theta_rad, 0.0, pol=pol)

    def d3n_wl(self, wl_um, theta_rad, pol='o'):
        return super().d3n_wl(wl_um, theta_rad, 0.0, pol=pol)
    
    def GD(self, wl_um, theta_rad, pol='o'):
        """Group Delay [fs/mm]"""
        return super().GD(wl_um, theta_rad, 0.0, pol=pol)
    
    def GV(self, wl_um, theta_rad, pol='o'):
        """Group Velocity [um/fs]"""
        return super().GV(wl_um, theta_rad, 0.0, pol=pol)
    
    def ng(self, wl_um, theta_rad, pol='o'):
        """Group index, c/Group velocity"""
        return super().ng(wl_um, theta_rad, 0.0, pol=pol)
    
    def GVD(self, wl_um, theta_rad, pol='o'):
        """Group Delay Dispersion [fs^2/mm]"""
        return super().GVD(wl_um, theta_rad, 0.0, pol=pol)
    
    def TOD(self, wl_um, theta_rad, pol='o'):
        """Third Order Dispersion [fs^3/mm]"""
        return super().TOD(wl_um, theta_rad, 0.0, pol=pol)