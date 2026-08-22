from ndispers._baseclass import Medium
import numpy as np

class Uniax_neg_3m(Medium):
    __slots__ = []
    """
    Class of negative uniaxial crystals of point group 3m.

    d = [[   0,   0,   0,   0, d15, -d22 ],
         [ -d22, d22,  0,  d15,  0,   0  ],
         [  d31, d31, d33,  0,   0,   0  ]]
      = [d_il]
        i = 1,2,3 (1=x, 2=y, 3=z)
        l = 1,2,...,6 (1=xx, 2=yy, 3=zz, 4=yz, 5=zx, 6=xy)

    Under the Kleinman symmetry,
    dxzx = dzxx, or d15 = d31.
    Nonzero independent elements are d22(=-d16), d31(=d15) and d33.
    (d33 does not come up in deff since its geometric factor is zero for all the interaction types.)

    """

    def __init__(self):
        super().__init__()

    def delta22(self, d22_known, wl1, wl2, T_degC):
        """
        Miller' delta for d22 coefficient (ijk = yyy).

        Note
        ----
        This value is almost constant and can be used to estimate wavelength dependence of d22
        from a known d22 of SFG with wl1 and wl2 as two pump wavelengths (Miller's rule).
        For yyy, two pump (wl1, wl2) and sum-frequency (wl3) waves are all ordinary waves.

        Parameters
        ----------
        d22_known : float, second-order nonlinear coefficient, d22, known from literature.
        wl1 : float, first pump wavelength in µm for the d22
        wl2 : float, second pump wavelength in µm for the d22
        T_degC : float, crystal temperature in degC.

        return
        ------
        float, Miller's delta for 22 (ijk = yyy) component.

        """
        wl3 = 1./(1./wl1 + 1./wl2)
        theta_rad = 0 # This is an arbitrary value since n of o-wave does not depends on theta.
        # chi : linear susceptibility
        chi_wl1o = self.n(wl1, theta_rad, T_degC, pol='o')**2 - 1
        chi_wl2o = self.n(wl2, theta_rad, T_degC, pol='o')**2 - 1
        chi_wl3o = self.n(wl3, theta_rad, T_degC, pol='o')**2 - 1
        return d22_known / (chi_wl3o * chi_wl2o * chi_wl1o)

    def d22_sfg(self, wl1, wl2, T_degC, delta22=0):
        """
        d22 coefficient as a function of two pump wavelengths for SFG.

        Parameters
        ----------
        wl1 : float, first pump wavelength in µm.
        wl2 : float, second pump wavelength in µm.
        T_degC : crystal temperature in degC.

        return
        ------
        float, d22 coefficient (pm/V).

        """

        wl3 = 1./(1./wl1 + 1./wl2)
        theta_rad = 0 # This is an arbitrary value since n of o-wave does not depends on theta.
        chi_wl1o = self.n(wl1, theta_rad, T_degC, pol='o')**2 - 1
        chi_wl2o = self.n(wl2, theta_rad, T_degC, pol='o')**2 - 1
        chi_wl3o = self.n(wl3, theta_rad, T_degC, pol='o')**2 - 1
        return delta22 * chi_wl3o * chi_wl2o * chi_wl1o
    
    def delta31(self, d31_known, wl1, wl2, T_degC):
        """
        Miller' delta for d31 coefficient (ijk = zxx).

        Note
        ----
        This value is almost constant and can be used to estimate wavelength dependence of d31
        from a known d31 of SFG with wl1 and wl2 as two pump wavelengths (Miller's rule).
        For zxx, two pump waves are ordinary and the SF wave is extra-ordinary.

        Parameters
        ----------
        d31_known : float, second-order nonlinear coefficient, d31, known from literature.
        wl1 : float, first pump wavelength in µm for the d31
        wl2 : float, second pump wavelength in µm for the d31
        T_degC : float, crystal temperature in degC.

        return
        ------
        float, Miller's delta for 31 (ijk = zxx) component.

        """
        wl3 = 1./(1./wl1 + 1./wl2)
        pmAngles = self.pmAngles_sfg(wl1, wl2, T_degC)
        theta_rad = pmAngles['ooe']['theta'][0] # rad
        chi_wl1o = self.n(wl1, theta_rad, T_degC, pol='o')**2 - 1
        chi_wl2o = self.n(wl2, theta_rad, T_degC, pol='o')**2 - 1
        chi_wl3e = self.n(wl3, theta_rad, T_degC, pol='e')**2 - 1
        return d31_known / (chi_wl3e * chi_wl2o * chi_wl1o)
    
    def d31_sfg(self, wl1, wl2, T_degC, delta31=0):
        """
        d31 coefficient as a function of two pump wavelengths for SFG.

        Parameters
        ----------
        wl1 : float, first pump wavelength in µm.
        wl2 : float, second pump wavelength in µm.
        T_degC : float, crystal temperature in degC.

        return
        ------
        float, d31 coefficient (pm/V).

        """

        wl3 = 1./(1./wl1 + 1./wl2)
        pmAngles = self.pmAngles_sfg(wl1, wl2, T_degC)
        theta_rad = pmAngles['ooe']['theta'][0] # rad
        chi_wl1o = self.n(wl1, theta_rad, T_degC, pol='o')**2 - 1
        chi_wl2o = self.n(wl2, theta_rad, T_degC, pol='o')**2 - 1
        chi_wl3e = self.n(wl3, theta_rad, T_degC, pol='e')**2 - 1
        return delta31 * chi_wl3e * chi_wl2o * chi_wl1o
    
    def deff_sfg(self, wl1, wl2, theta_rad, phi_rad, T_degC, pol1, pol2, pol3):
        """
        Effective second-order nonlinear coefficient, d_eff, of SFG.
        Wavelength dependence is estimated from Miller's rule.

        Parameters
        ----------
        wl1  :  float, first pump wavelength in µm.
        wl2  : float, second pump wavelength in µm.
        T_degC : float, crystal temperature in degC.
        pol1 : str, {'o', 'e'}, polarization of the 1st pump wave.
        pol2 : str, {'o', 'e'}, polarization of the 2nd pump wave.
        pol3 : str, {'o', 'e'}, polarization of the sum-frequency wave.

        return
        ------
        float, Effective d coefficient (pm/V).

        """
        wl3 = 1./(1./wl1 + 1./wl2)
        _d31 = self.d31_sfg(wl1, wl2, T_degC)
        _d22 = self.d22_sfg(wl1, wl2, T_degC)
        rho3 = self.woa_theta(wl3, theta_rad, T_degC, pol=pol3)
        # For negative uniaxial type-I (ooe),
        if pol1 == 'o' and pol2 == 'o' and pol3 == 'e':
            deff = _d31 * np.sin(theta_rad + rho3) - _d22 * np.cos(theta_rad + rho3) * np.sin(3*phi_rad)
        
        # For negative uniaxial type-II (oee),
        elif pol1 == 'o' and pol2 == 'e' and pol3 == 'e':
            rho2 = self.woa_theta(wl2, theta_rad, T_degC, pol=pol2)
            deff = _d22 * np.cos(theta_rad + rho2) * np.cos(theta_rad + rho3) * np.cos(3*phi_rad)
        
        # For negative uniaxial type-II (eoe),
        elif pol1 == 'e' and pol2 == 'o' and pol3 == 'e':
            rho1 = self.woa_theta(wl1, theta_rad, T_degC, pol=pol1)
            deff = _d22 * np.cos(theta_rad + rho1) * np.cos(theta_rad + rho3) * np.cos(3*phi_rad)
        
        else:
            raise ValueError("For negative uniaxial crystals,pol1, pol2, pol3 must be among {ooe, oee, eoe}.")

        return deff
