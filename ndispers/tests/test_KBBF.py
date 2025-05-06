""" 
Compare refractive indices to the experimentally measured values.
KBBF (KBe2BO3F2)
"""
from numpy import pi
from ndispers.media.crystals import _KBBF_Li2016

modules = [_KBBF_Li2016]

def sub(module):
    module_name = module.__name__
    print("[[ Module : %s ]]" % module_name)

    kbbf = module.KBBF()

    # From Table 1 in Li et al. 2016, Applied Optics, 55(36), 10423-10426
    # Wavelength (µm)
    wl_list = [
        0.1870, 0.1900, 0.1930, 0.2000, 0.2357, 0.3650, 0.4067, 
        0.4358, 0.4800, 0.5461, 0.5876, 0.6438, 0.7065, 0.8521, 
        1.014, 1.530
    ]
    
    # Experimental values for ordinary ray
    n_o_list = [
        1.59005, 1.58444, 1.57995, 1.56922, 1.52622, 1.49600, 1.49128, 
        1.48844, 1.48538, 1.48208, 1.48048, 1.47875, 1.47727, 1.47456, 
        1.47213, 1.46544
    ]
    
    # Experimental values for extraordinary ray
    n_e_list = [
        1.46299, 1.45991, 1.45718, 1.45132, 1.42962, 1.40641, 1.40336, 
        1.40135, 1.39933, 1.39737, 1.39633, 1.39524, 1.39430, 1.39274, 
        1.39142, 1.38821
    ]

    T_degC = 22.0  # Temperature specified in the paper
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
        n_o_calc = kbbf.n(wl, 0, T_degC, pol='o')
        n_e_calc = kbbf.n(wl, pi*0.5, T_degC, pol='e')
        dndT_o = kbbf.dndT(wl, 0, T_degC, pol='o')
        dndT_e = kbbf.dndT(wl, 0, T_degC, pol='e')

        print("%.4f" % wl)
        print("  |  {: .5f}   {: .5f}  |".format(n_o_exp, n_e_exp))
        print("  |  {: .5f}   {: .5f}  |  {: .4e}   {: .4e}".format(n_o_calc, n_e_calc, dndT_o, dndT_e))
        print("  | ({: .5f}) ({: .5f}) |".format(n_o_calc - n_o_exp, n_e_calc - n_e_exp))

    print("-"*80)
    
    # Test phase-matching angles for Type I SHG (from Table 2)
    print("\nPhase-matching angles for Type I SHG")
    print("-"*80)
    print("Fundamental Wavelength(µm) | SHG Wavelength(µm) | Exp Angle(deg) | Cal Angle(deg) | Difference")
    print("-"*80)
    
    # Data from Table 2 in the paper
    # Format: [fundamental wavelength (nm), experimental angle (deg)]
    pm_data = [
        [340.0, 70.4],
        [345.0, 68.3],
        [352.4, 65.3],
        [354.7, 64.5],
        [361.2, 62.6],
        [363.4, 61.1],
        [367.4, 60.0],
        [369.5, 61.0],
        [374.3, 59.4],
        [378.6, 57.6],
        [384.7, 56.8],
        [388.4, 55.1],
        [410.0, 51.0],
        [440.0, 46.0],
        [460.0, 44.0],
        [480.0, 41.7],
        [500.0, 39.9],
        [520.0, 36.9],
        [550.0, 34.9],
        [589.0, 32.7],
        [600.0, 32.1],
        [680.0, 27.7],
        [770.0, 25.1],
        [850.0, 23.0],
        [900.0, 22.0],
        [950.0, 21.4],
        [1064, 20.2],
        [1200, 19.7],
        [1300, 19.3],
        [1400, 19.3]
    ]
    
    for data in pm_data:
        fund_wl_nm, exp_angle = data
        fund_wl = fund_wl_nm / 1000  # Convert to µm
        shg_wl_nm = fund_wl_nm / 2   # SHG wavelength is half the fundamental
        
        # Calculate phase-matching angle using the library
        # For Type I SHG (o + o -> e), we use the 'ooe' key
        pm_angles = kbbf.pmAngles_sfg(fund_wl, fund_wl, T_degC, tol_deg=0.005, deg=True)
        cal_angle = pm_angles['ooe']['theta'][0] if pm_angles['ooe']['theta'] else None
        
        if cal_angle is not None:
            diff = cal_angle - exp_angle
            print("{:.1f} | {:.1f} | {:.2f} | {:.2f} | {:.2f}".format(
                fund_wl_nm, shg_wl_nm, exp_angle, cal_angle, diff))
        else:
            print("{:.1f} | {:.1f} | {:.2f} | No solution found | N/A".format(
                fund_wl_nm, shg_wl_nm, exp_angle))
    
    print("-"*80)

def main():
    print("="*80)
    print("Compare refractive indices to the experimentally measured values.")
    print("Ref: %s\n"  % ('Li, R., et al. "Dispersion relations of refractive indices suitable for KBe2BO3F2 crystal deep-ultraviolet applications." Applied Optics, 55(36), 10423-10426 (2016)'))
    for module in modules:
        sub(module)
    print("="*80)

if __name__ == "__main__":
    main()
