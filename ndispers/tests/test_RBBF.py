""" 
Compare refractive indices, thermo-optic coefficients, and phase-matching angles to the experimentally measured values.
RBBF (RbBe2BO3F2)
"""
import numpy as np
from numpy import pi

from ndispers.media.crystals import _RBBF_Chen2009

modules = [_RBBF_Chen2009]

def sub(module):
    module_name = module.__name__
    print("[[ Module : %s ]]" % module_name)

    rbbf = module.RBBF()

    # From Table 2 in Chen et al. 2009, JOSAB 26(8), 1519-1525
    # Wavelength (nm)
    wl_list =   [404.7, 435.8, 486.1, 491.6, 546.1, 577.0, 589.3, 656.3, 694.3]
    # Convert to µm
    wl_list = [wl/1000 for wl in wl_list]
    
    # Experimental values for ordinary ray
    n_o_list = [1.49761, 1.49469, 1.49128, 1.49092, 1.48827, 1.48706, 1.48636, 1.48468, 1.48384]
    # Experimental values for extraordinary ray
    n_e_list = [1.41956, 1.41748, 1.41511, 1.41493, 1.41314, 1.41238, 1.41178, 1.41071, 1.41011]

    T_degC = 20.0
    print("T=%.1f degC"  % T_degC)

    print("-"*80)
    print("Wavelength(µm)")
    print("  |  n_o        n_e     |    dn_o/dT       dn_e/dT")
    print("  |  Experimental")
    print("  |  Calculated")
    print("  |  (Cal - Exp)")
    print("-"*80)
    for i, wl in enumerate(wl_list):
        # measured
        n_o_exp = n_o_list[i]
        n_e_exp = n_e_list[i]
        # calculated
        n_o_calc = rbbf.n(wl, 0, T_degC, pol='o')
        n_e_calc = rbbf.n(wl, pi*0.5, T_degC, pol='e')
        dndT_o = rbbf.dndT(wl, 0, T_degC, pol='o')
        dndT_e = rbbf.dndT(wl, 0, T_degC, pol='e')

        print("%.4f" % wl)
        print("  |  {: .5f}   {: .5f}  |".format(n_o_exp, n_e_exp))
        print("  |  {: .5f}   {: .5f}  |  {: .4e}   {: .4e}".format(n_o_calc, n_e_calc, dndT_o, dndT_e))
        print("  | ({: .5f}) ({: .5f}) |".format(n_o_calc - n_o_exp, n_e_calc - n_e_exp))

    print("-"*80)
    
    # Test phase-matching angles for Type I SHG (from Table 3)
    print("\nPhase-matching angles for Type I SHG")
    print("-"*80)
    print("Fundamental Wavelength(µm) | SHG Wavelength(µm) | Exp Angle(deg) | Cal Angle(deg) | Difference")
    print("-"*80)
    
    # Data from Table 3 in the paper
    # Format: [fundamental wavelength (nm), SHG wavelength (nm), experimental angle (deg)]
    pm_data = [
        [354.7, 177.3, 73.38],
        [360.0, 180.0, 70.31],
        [365.0, 182.5, 68.64],
        [370.0, 185.0, 66.92],
        [375.0, 187.5, 65.04],
        [380.0, 190.0, 63.51],
        [385.0, 192.5, 61.85],
        [390.0, 195.0, 60.60],
        [395.0, 197.5, 59.49],
        [400.0, 200.0, 58.04],
        [405.0, 202.5, 56.94],
        [410.0, 205.0, 55.83],
        [415.0, 207.5, 54.88],
        [420.0, 210.0, 54.14],
        [425.0, 212.5, 53.39],
        [430.0, 215.0, 52.37],
        [435.0, 217.5, 51.58],
        [440.0, 220.0, 50.81],
        [515.0, 257.5, 41.17],
        [529.6, 264.8, 39.86],
        [532.0, 266.0, 39.97],
        [549.7, 274.9, 38.24],
        [570.2, 285.1, 36.74],
        [589.7, 294.9, 35.38],
        [610.0, 305.0, 34.00],
        [629.7, 314.9, 32.89],
        [664.5, 332.3, 31.38],
        [730.0, 365.0, 28.55],
        [740.0, 370.0, 28.20],
        [750.0, 375.0, 27.82],
        [750.1, 375.1, 27.83],
        [760.0, 380.0, 27.55],
        [760.8, 380.4, 27.56],
        [770.0, 385.0, 27.16],
        [780.0, 390.0, 26.93],
        [790.0, 395.0, 26.55],
        [799.7, 399.9, 26.30],
        [800.0, 400.0, 26.26],
        [810.0, 405.0, 26.04],
        [812.2, 406.1, 26.12],
        [820.0, 410.0, 25.81],
        [830.0, 415.0, 25.52],
        [840.0, 420.0, 25.28],
        [849.4, 424.7, 25.05],
        [850.0, 425.0, 25.06],
        [860.0, 430.0, 24.80],
        [870.0, 435.0, 24.66],
        [880.0, 440.0, 24.36],
        [897.7, 448.9, 23.93],
        [949.8, 474.9, 23.22],
        [997.7, 498.9, 22.53],
        [1064.0, 532.0, 21.62],
        [1109.0, 554.5, 21.42],
        [1203.1, 601.6, 20.90],
        [1299.2, 649.6, 20.44],
        [1399.5, 699.8, 20.34]
    ]
    
    for data in pm_data:
        fund_wl_nm, shg_wl_nm, exp_angle = data
        fund_wl = fund_wl_nm / 1000  # Convert to µm
        
        # Calculate phase-matching angle using the library
        # For Type I SHG (o + o -> e), we use the 'ooe' key
        pm_angles = rbbf.pmAngles_sfg(fund_wl, fund_wl, T_degC, tol_deg=0.005, deg=True)
        cal_angle = pm_angles['ooe']['theta'][0] if pm_angles['ooe']['theta'] else None
        
        if cal_angle is not None:
            diff = cal_angle - exp_angle
            print("{:.1f} | {:.1f} | {:.2f} | {:.2f} | {:.2f}".format(
                fund_wl_nm, shg_wl_nm, exp_angle, cal_angle, diff))
        else:
            print("{:.1f} | {:.1f} | {:.2f} | No solution found | N/A".format(
                fund_wl_nm, shg_wl_nm, exp_angle))
    
    print("-"*80)

def test_dndT(module):
    """Test thermo-optic coefficients (dn/dT) against experimental values from Zhai et al. 2013"""
    module_name = module.__name__
    print("[[ Module : %s - Thermo-optic Coefficients ]]" % module_name)
    
    rbbf = module.RBBF()
    
    # Data from Table 2 in Zhai et al. 2013, Optical Materials 36(2), 333-336
    # Format: wavelength (µm), experimental dn/dT for o-ray (×10⁻⁶), calculated dn/dT for o-ray (×10⁻⁶),
    #        experimental dn/dT for e-ray (×10⁻⁶), calculated dn/dT for e-ray (×10⁻⁶)
    dndT_data = [
        [0.194, -6.78, -6.78, -7.35, -7.35],
        [0.254, -10.00, -9.99, -9.98, -9.97],
        [0.363, -11.29, -11.36, -9.91, -9.97],
        [0.405, -11.54, -11.54, -9.95, -9.91],
        [0.435, -11.67, -11.63, -9.92, -9.91],
        [0.546, -11.85, -11.84, -10.09, -10.06],
        [0.644, -11.95, -11.96, -10.30, -10.33],
        [0.706, -12.02, -12.03, -10.52, -10.52],
        [1.014, -12.27, -12.27, -11.48, -11.47]
    ]
    
    T_degC = 24.0  # Reference temperature from the paper
    
    print("T=%.1f degC" % T_degC)
    print("-"*100)
    print("Wavelength(µm) | dn_o/dT (×10⁻⁶) |                | dn_e/dT (×10⁻⁶) |               ")
    print("              | Experimental   | Calculated     | Experimental   | Calculated    ")
    print("              |                | (Difference)   |                | (Difference)  ")
    print("-"*100)
    
    for data in dndT_data:
        wl, dndT_o_exp, dndT_o_cal_paper, dndT_e_exp, dndT_e_cal_paper = data
        
        # Calculate dn/dT using our implementation
        dndT_o_calc = rbbf.dndT(wl, 0, T_degC, pol='o') * 1e6  # Convert to ×10⁻⁶
        dndT_e_calc = rbbf.dndT(wl, pi/2, T_degC, pol='e') * 1e6  # Convert to ×10⁻⁶
        
        # Calculate differences
        diff_o = dndT_o_calc - dndT_o_exp
        diff_e = dndT_e_calc - dndT_e_exp
        
        print("{:14.3f} | {:14.2f} | {:14.2f} | {:14.2f} | {:14.2f}".format(
            wl, dndT_o_exp, dndT_o_calc, dndT_e_exp, dndT_e_calc))
        print("              |                | ({:+10.2f})   |                | ({:+10.2f})  ".format(
            diff_o, diff_e))
    
    print("-"*100)
    
    # Calculate average absolute differences
    avg_diff_o = np.mean([abs(data[1] - rbbf.dndT(data[0], 0, T_degC, pol='o') * 1e6) for data in dndT_data])
    avg_diff_e = np.mean([abs(data[3] - rbbf.dndT(data[0], pi/2, T_degC, pol='e') * 1e6) for data in dndT_data])
    
    print("Average absolute difference (o-ray): {:.4f} × 10⁻⁶".format(avg_diff_o))
    print("Average absolute difference (e-ray): {:.4f} × 10⁻⁶".format(avg_diff_e))
    print()

def main():
    print("="*100)
    print("Compare refractive indices, thermo-optic coefficients, and phase-matching angles to the experimentally measured values.")
    print("Refs:\n  - Chen, C., et al. \"Growth, properties, and application to nonlinear optics of a nonlinear optical crystal: RbBe2BO3F2.\" Journal of the Optical Society of America B, 26(8), 1519-1525 (2009)")
    print("  - Zhai, N., et al. \"Measurement of thermal refractive index coefficients of nonlinear optical crystal RbBe2BO3F2.\" Optical Materials, 36(2), 333-336 (2013)")
    print()
    
    for module in modules:
        sub(module)
        test_dndT(module)
    print("="*100)

if __name__ == "__main__":
    main()
