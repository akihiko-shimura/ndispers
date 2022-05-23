""" 
Compare refractive indices to the experimentaly measured.
LBO
"""
from ndispers.media.crystals import (_LBO_Castech, _LBO_Ghosh1995, _LBO_Newlight, _LBO_KK1994, _LBO_KK2018)
modules = [_LBO_Castech, _LBO_Ghosh1995, _LBO_Newlight, _LBO_KK1994, _LBO_KK2018]

def sub(module):
    module_name = module.__name__
    print("[[ Module : %s ]]" % module_name)

    lbo_yz = module.LBO_yz() # n_x
    lbo_zx = module.LBO_zx() # n_y
    lbo_xy = module.LBO_xy() # n_z

    # from Nikogosyan, Nonlinear Optical Crystals, LBO ref[1]
    wl_list =  [0.2537, 0.3341, 0.5321, 1.0642]
    n_x_list = [1.6335, 1.6043, 1.5785, 1.5656]
    n_y_list = [1.6582, 1.6346, 1.6065, 1.5905]
    n_z_list = [1.6792, 1.6509, 1.6212, 1.6055]

    T_degC = 20.0
    print("T=%.1f degC"  % T_degC)

    print("-"*80)
    print("Wavelength(µm)")
    print("  |  n_x        n_y        n_z    |    dn_x/dT       dn_y/dT       dn_z/dT")
    print("  |  Experimental")
    print("  |  Calculated")
    print("  |  (Cal - Exp)")
    print("-"*80)
    for i, wl in enumerate(wl_list):
        # measured
        n_x_exp = n_x_list[i]
        n_y_exp = n_y_list[i]
        n_z_exp = n_z_list[i]
        # calculated
        n_x_calc = lbo_yz.n(wl, 0, T_degC, pol='o')
        n_y_calc = lbo_zx.n(wl, 0, T_degC, pol='o')
        n_z_calc = lbo_xy.n(wl, 0, T_degC, pol='o')

        dndT_x = lbo_yz.dndT(wl, 0, T_degC, pol='o')
        dndT_y = lbo_zx.dndT(wl, 0, T_degC, pol='o')
        dndT_z = lbo_xy.dndT(wl, 0, T_degC, pol='o')

        print("%.4f" % wl)
        print("  |  {: .4f}   {: .4f}   {: .4f}  |".format(n_x_exp, n_y_exp, n_z_exp))
        print("  |  {: .4f}   {: .4f}   {: .4f}  |  {: .4e}   {: .4e}   {: .4e}".format(n_x_calc, n_y_calc, n_z_calc, dndT_x, dndT_y, dndT_z))
        print("  | ({: .4f}) ({: .4f}) ({: .4f}) |".format(n_x_calc - n_x_exp, n_y_calc - n_y_exp, n_z_calc - n_z_exp))

    print("-"*80)

def main():
    print("="*80)
    print("Compare refractive indices to the experimentaly measured values.")
    print("Ref: %s\n" %('Nikogosyan, D.N., Nonlinear Optical Crystals: A Complete Survey (Springer, 2005)'))
    for module in modules:
        sub(module)

if __name__ == "__main__":
    main()

    