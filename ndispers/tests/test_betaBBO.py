""" 
Compare refractive indices to the experimentaly measured.
beta-BBO
"""
from ndispers.media.crystals import (_betaBBO1987, _betaBBO1995, _betaBBO2010, _betaBBO2018)
modules = [_betaBBO1987, _betaBBO1995, _betaBBO2010, _betaBBO2018]
from numpy import pi

def sub(module):
    module_name = module.__name__
    print("[[ Module : %s ]]" % module_name)

    bbo = module.BetaBBO()

    # from Nikogosyan, D. N. "Beta barium borate (BBO)." Applied Physics A 52.6 (1991): 359-368.
    wl_list =   [0.2128, 0.2660, 0.3541, 0.5321, 1.0642]
    n_o_list = [1.84707, 1.75707, 1.70556, 1.67493, 1.65510]
    n_e_list = [1.67467, 1.61461, 1.57757, 1.55552, 1.54254]

    T_degC = 20.0
    print("T=%.1f degC"  % T_degC)

    print("-"*36)
    print("Wavelength(um)")
    print("  |  n_o        n_e")
    print("  |  Experimental")
    print("  |  Calculated")
    print("  |  (Cal - Exp)")
    print("-"*36)
    for i, wl in enumerate(wl_list):
        # measured
        n_o_exp = n_o_list[i]
        n_e_exp = n_e_list[i]
        # calculated
        n_o_calc = bbo.n(wl, 0, T_degC, pol='o')
        n_e_calc = bbo.n(wl, pi*0.5, T_degC, pol='e')

        print("%.4f" % wl)
        print("  |  {: .4f}   {: .4f}".format(n_o_exp, n_e_exp))
        print("  |  {: .4f}   {: .4f}".format(n_o_calc, n_e_calc))
        print("  | ({: .4f}) ({: .4f})".format(n_o_calc - n_o_exp, n_e_calc - n_e_exp))

    print("-"*36)

def main():
    print("="*42)
    print("Compare refractive indices to the experimentaly measured values.")
    print("Ref: %s\n"  % ('Nikogosyan, D. N. "Beta barium borate (BBO)." Applied Physics A 52.6 (1991): 359-368.'))
    for module in modules:
        sub(module)

if __name__ == "__main__":
    main()

    