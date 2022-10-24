from ndispers._baseclass import Medium
import numpy as np

class Uniax_neg_32(Medium):
    """
    Class of negative uniaxial crystals of point group 32.

    d = [[ d11, -d11,  0, d14,   0,    0  ],
         [  0,    0,   0,  0,  -d14, -d11 ],
         [  0,    0,   0,  0,    0,    0  ]]
      = [d_il]
        i = 1,2,3 (1=x, 2=y, 3=z)
        l = 1,2,...,6 (1=xx, 2=yy, 3=zz, 4=yz, 5=zx, 6=xy)
    
    Under the Kleinman symmetry,
    dxyz = dyzx, or d14 = d25.
    
    With d25 = - d14, we get d14 = d25 = -d14. Therefore, d14 = 0.
    Nonzero element is only d11.

    """

    def __init__(self):
        super().__init__()
    
    def delta11(self, d11_known, wl1, wl2, T_degC):
        """
        Miller' delta for d11 coefficient (ijk = xxx).

        Parameters
        ----------
        d11_known : float, second-order nonlinear coefficient, d11, known from literature.
        wl1 : float, first pump wavelength in µm for the d11
        wl2 : float, second pump wavelength in µm for the d11
        T_degC : float, crystal temperature in degC.

        return
        ------
        float, Miller's delta for 11 (ijk = xxx) component.

        """
        wl3 = 1./(1./wl1 + 1./wl2)
        theta_rad = 0 # This is an arbitrary value since n of o-wave does not depends on theta.
        # chi : linear susceptibility
        chi_wl1o = self.n(wl1, theta_rad, T_degC, pol='o')**2 - 1
        chi_wl2o = self.n(wl2, theta_rad, T_degC, pol='o')**2 - 1
        chi_wl3o = self.n(wl3, theta_rad, T_degC, pol='o')**2 - 1
        return d11_known / (chi_wl3o * chi_wl2o * chi_wl1o)

    def d11_sfg(self, wl1, wl2, T_degC, delta11=0):
        """
        d11 coefficient as a function of two pump wavelengths for SFG.

        Parameters
        ----------
        wl1 : float, first pump wavelength in µm.
        wl2 : float, second pump wavelength in µm.
        T_degC : crystal temperature in degC.

        return
        ------
        float, d11 coefficient (pm/V).

        """

        wl3 = 1./(1./wl1 + 1./wl2)
        theta_rad = 0 # This is an arbitrary value since n of o-wave does not depends on theta.
        chi_wl1o = self.n(wl1, theta_rad, T_degC, pol='o')**2 - 1
        chi_wl2o = self.n(wl2, theta_rad, T_degC, pol='o')**2 - 1
        chi_wl3o = self.n(wl3, theta_rad, T_degC, pol='o')**2 - 1
        return delta11 * chi_wl3o * chi_wl2o * chi_wl1o
    
    def delta14(self, d14_known, wl1, wl2, T_degC):
        """
        Miller' delta for d14 coefficient (ijk = xyz).

        Parameters
        ----------
        d14_known : float, second-order nonlinear coefficient, d14, known from literature.
        wl1 : float, first pump wavelength in µm for the d14
        wl2 : float, second pump wavelength in µm for the d14
        T_degC : float, crystal temperature in degC.

        return
        ------
        float, Miller's delta for 14 (ijk = xyz) component.

        """
        wl3 = 1./(1./wl1 + 1./wl2)
        theta_rad = 0 # This is an arbitrary value since n of o-wave does not depends on theta.
        # chi : linear susceptibility
        chi_wl1o = self.n(wl1, theta_rad, T_degC, pol='o')**2 - 1
        chi_wl2o = self.n(wl2, theta_rad, T_degC, pol='o')**2 - 1
        chi_wl3o = self.n(wl3, theta_rad, T_degC, pol='o')**2 - 1
        return d14_known / (chi_wl3o * chi_wl2o * chi_wl1o)

    def d14_sfg(self, wl1, wl2, T_degC, delta14=0):
        """
        d11 coefficient as a function of two pump wavelengths for SFG.

        Parameters
        ----------
        wl1 : float, first pump wavelength in µm.
        wl2 : float, second pump wavelength in µm.
        T_degC : crystal temperature in degC.

        return
        ------
        float, d14 coefficient (pm/V).

        """

        wl3 = 1./(1./wl1 + 1./wl2)
        theta_rad = 0 # This is an arbitrary value since n of o-wave does not depends on theta.
        chi_wl1e = self.n(wl1, theta_rad, T_degC, pol='e')**2 - 1
        chi_wl2o = self.n(wl2, theta_rad, T_degC, pol='o')**2 - 1
        chi_wl3e = self.n(wl3, theta_rad, T_degC, pol='e')**2 - 1
        return delta14 * chi_wl3e * chi_wl2o * chi_wl1e
    
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
        _d11 = self.d11_sfg(wl1, wl2, T_degC)
        _d14 = self.d14_sfg(wl1, wl2, T_degC)
        rho3 = self.woa_theta(wl3, theta_rad, T_degC, pol=pol3)
        # For negative uniaxial type-I (ooe),
        if pol1 == 'o' and pol2 == 'o' and pol3 == 'e':
            deff = _d11 * np.cos(3*phi_rad) * np.cos(theta_rad + rho3)
        
        # For negative uniaxial type-II (oee),
        elif pol1 == 'o' and pol2 == 'e' and pol3 == 'e':
            rho2 = self.woa_theta(wl2, theta_rad, T_degC, pol=pol2)
            deff = _d11 * np.sin(3*phi_rad) * np.cos(theta_rad + rho2) * np.cos(theta_rad + rho3)\
                + _d14 * np.sin(theta_rad + rho2) * np.cos(theta_rad + rho3)
        
        # For negative uniaxial type-II (eoe),
        elif pol1 == 'e' and pol2 == 'o' and pol3 == 'e':
            rho1 = self.woa_theta(wl1, theta_rad, T_degC, pol=pol1)
            deff = _d11 * np.sin(3*phi_rad) * np.cos(theta_rad + rho1) * np.cos(theta_rad + rho3)\
                + _d14 * np.sin(theta_rad + rho1) * np.cos(theta_rad + rho3)
        
        else:
            raise ValueError("For negative uniaxial crystals,pol1, pol2, pol3 must be among {ooe, oee, eoe}.")

        return deff
