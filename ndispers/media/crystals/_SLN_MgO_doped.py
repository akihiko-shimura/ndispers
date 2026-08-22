import sympy

from ndispers._baseclass import T, phi, theta, wl
from ndispers.groups import Uniax_neg_3m
from ndispers.helper import vars2

class SLN(Uniax_neg_3m):
    """
    1% MgO-doped stoichiometric Lithium niobate (Li Nb O_3) crystal

    - Point group : 3m  (C_{3v})
    - Crystal system : Trigonal
    - Dielectric principal axis, z // c-axis (x, y-axes are arbitrary)
    - Negative uniaxial, with optic axis parallel to z-axis
    - Transparency range : 0.32 µm to 5.2 µm

    Sellmeier equation
    ------------------
        n(wl) = sqrt(a1_i + b1 * f + (a2_i + b2_i * f)/(wl**2 - (a3_i + b3_i * f)**2) + (a4_i + b4_i * f)/(wl**2 - a5_i**2) - a6_i * wl**2) for i=o,e
        f = (T - T0) * (T + T0 + 2 * 273.16) with T0 = 24.5 degC

    Validity range
    --------------
    0.5 to 4 µm (range of the Gayer 2008 coefficients); temperature
    dependence measured up to ~200 degC in the paper.

    Note
    ----
    Sellmeier equation only for e-wave is given.

    Gayer's coefficients were fitted to quasi-phase-matching measurements, i.e. to
    index *differences*. Phase matching is reproduced well (this equation gives
    1544 nm at 40 degC and 1571 nm at 200 degC for SHG in a 19.36 µm period
    crystal, matching Fig. 4 of the paper), but the absolute thermo-optic
    coefficient is not constrained by that data: dn_e/dT comes out around
    +5e-4 /K at 1.064 µm, an order of magnitude above the accepted ~3e-5 /K for
    LiNbO3. The same holds for the other Gayer sets. Use dndT/dndT2 from this
    class for phase-matching work, not as absolute thermo-optic values.

    Ref
    ---
    Gayer, O., et al. "Temperature and wavelength dependent refractive index equations for MgO-doped congruent and stoichiometric LiNbO3." Applied Physics B 91.2 (2008): 343-348.
    Erratum: Applied Physics B 101.2 (2010): 481. (b1 corrected to 4.677e-6)
    https://www.opt-oxide.com/products/sln/

    """
    __slots__ = ["_a1_o", "_a2_o", "_a3_o", "_a4_o",  "_a5_o", "_a6_o",
                 "_a1_e", "_a2_e", "_a3_e", "_a4_e",  "_a5_e", "_a6_e",
                 "_b1_o", "_b2_o", "_b3_o", "_b4_o",
                 "_b1_e", "_b2_e", "_b3_e", "_b4_e",
                 "_d31_1064shg", "_d22_1064shg"]
                 
    _default_pol = 'e'   # no o-ray Sellmeier set exists for SLN

    def __init__(self):
        super().__init__()
        self._plane = 'arb'
        self._theta_rad = 'var'
        self._phi_rad = 'arb'

        """ Constants of dispersion formula """
        # 1% MgO-doped SLN (e-ray only; Gayer 2008 gives no o-ray set for SLN)
        self._a1_e = 5.078
        self._a2_e = 0.0964
        self._a3_e = 0.2065
        self._a4_e = 61.16
        self._a5_e = 10.55
        self._a6_e = 1.59e-2
        # b1 is 4.677e-6, not the 4.677e-7 printed in Table 1 of the 2008 paper:
        # that single value is what the 2010 erratum exists to correct.
        self._b1_e = 4.677e-6
        self._b2_e = 7.822e-8
        self._b3_e = -2.653e-8
        self._b4_e = 1.096e-4

        # Second-order nonlinear optical coefficients
        self._d31_1064shg = 4.4 #pm/V
        self._d22_1064shg = 25 #pm/V

    
    
    def n_e_expr(self):
        """ Sympy expression, dispersion formula for e-wave """
        return sympy.sqrt( self._a1_e + self._b1_e * self.f_expr() + \
            (self._a2_e + self._b2_e * self.f_expr()) / (wl**2 - (self._a3_e + self._b3_e * self.f_expr())**2) + \
                (self._a4_e + self._b4_e * self.f_expr()) / (wl**2 - self._a5_e**2) - self._a6_e * wl**2 )

    def f_expr(self):
        return (T - 24.5) * (T + 24.5 + 2 * 273.16)

    def n_expr(self, pol):
        """
        Sympy expression, 
        dispersion formula,
        only for e-wave

        """
        if pol == 'e':
            return self.n_e_expr()
        else:
            raise ValueError("pol = '%s' must be 'e'. Sellmeier equation for pol='o' is not implemented for this module." % pol)
    

    

    
    
    
    
    
    
    
    
    
    #------------------------------------------------------------------------------------------
    # Wavelength dependence of second-order nonlinear coefficients estimated from Miller's rule
    #------------------------------------------------------------------------------------------
    def d22_sfg(self, wl1o, wl2o, T_degC):
        return super().d22_sfg(wl1o, wl2o, T_degC, delta22=self.delta22(self._d22_1064shg, 1.064, 1.064, T_degC))

    def d31_sfg(self, wl1o, wl2o, T_degC):
        return super().d31_sfg(wl1o, wl2o, T_degC, delta31=self.delta31(self._d31_1064shg, 1.064, 1.064, T_degC))