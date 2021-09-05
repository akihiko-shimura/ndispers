""" ndispers - test.py """

import numpy as np
from  ndispers import media
# from . import __init__
pi = np.pi

def main():
    # uniaxial
    nlc = media.crystals.BetaBBO1987()
    nlc.help
    nlc.constants
    nlc.angles

    pol = 'o'
    print("\npol=", pol)
    wl_um = 0.6 #um
    print("wl_um =", wl_um)
    ans = nlc.n(wl_um, 0.0, 0.0, pol=pol)
    print("nlc.n(wl_um, 0.0, 0.0, pol='o') =", ans)
    print("type :", type(ans))

    pol = 'e'
    print("\npol=", pol)
    wl_um = 0.6 #um
    print("wl_um =", wl_um)
    ans = nlc.n(wl_um, 0.0, 0.0, pol=pol)
    print("nlc.n(wl_um, 0.0, 0.0, pol='e') =", ans)
    print("type :", type(ans))

    wl_ar = np.arange(0.18, 1.1, 0.1)
    print("\nwl_ar =", wl_ar)
    ans = nlc.n(wl_ar, 0.0, 0.0, pol='o')
    print("nlc.n(wl_ar, 0.0, 0.0, pol='o') =", ans)
    print("type :", type(ans))

    print("\nwl_ar =", wl_ar)
    ans = nlc.n(wl_ar, 0.0, 0.0, pol='e')
    print("nlc.n(wl_ar, 0.0, 0.0, pol='e') =", ans)
    print("type :", type(ans))


    th_ar = np.arange(0, pi, 0.1*pi)
    print("\nth_ar =", th_ar)
    ans = nlc.n(wl_um, th_ar, 0.0, pol='o')
    print("nlc.n(wl_um, th_ar, 0.0, pol='o') =",ans)
    print("type :", type(ans))

    print("\nth_ar =", th_ar)
    ans = nlc.n(wl_um, th_ar, 0.0, pol='e')
    print("nlc.n(wl_um, th_ar, 0.0, pol='e') =", ans)
    print("type :", type(ans))

    ph_ar = np.arange(0, pi, 0.1*pi)
    print("\nph_ar =", ph_ar)
    ans = nlc.n(wl_um, 0.0, ph_ar, pol='o')
    print("nlc.n(wl_um, 0.0, ph_ar, pol='o') =", ans)
    print("type :", type(ans))

    print("\nph_ar =", ph_ar)
    ans = nlc.n(wl_um, 0.0, ph_ar, pol='e')
    print("nlc.n(wl_um, 0.0, ph_ar, pol='e') =", ans)
    print("type :", type(ans))

    nlc._cached_func_dict

    # biaxial
    print("="*40)
    plane = "yz"
    print("\nplane=", plane)
    nlc = media.crystals.KTP_zx()
    nlc.help
    nlc.constants
    nlc.angles

    pol = 'o'
    print("pol=", pol)
    wl_um = 0.6 #um
    print("wl_um =", wl_um)
    ans = nlc.n(wl_um, 0.0, pol=pol)
    print("nlc.n(wl_um, 0.0, pol='o') =", ans)
    print("type :", type(ans))

    pol = 'e'
    print("pol=", pol)
    wl_um = 0.6 #um
    print("wl_um =", wl_um)
    ans = nlc.n(wl_um, 0.0, pol=pol)
    print("nlc.n(wl_um, 0.0, pol='o') =", ans)
    print("type :", type(ans))

    wl_ar = np.arange(0.18, 1.1, 0.1)
    print("\nwl_ar =", wl_ar)
    ans = nlc.n(wl_ar, 0.0, pol='o')
    print("nlc.n(wl_ar, 0.0, pol='o') =", ans)
    print("type :", type(ans))

    print("\nwl_ar =", wl_ar)
    ans = nlc.n(wl_ar, 0.0, pol='e')
    print("nlc.n(wl_ar, 0.0, pol='e') =", ans)
    print("type :", type(ans))

    th_ar = np.arange(0, pi, 0.1*pi)
    print("\nth_ar or phi_ar =", th_ar)
    ans = nlc.n(wl_um, th_ar, pol='o')
    print("nlc.n(wl_um, th_ar, pol='o') =",ans)
    print("type :", type(ans))

    th_ar = np.arange(0, pi, 0.1*pi)
    print("\nth_ar or phi_ar =", th_ar)
    ans = nlc.n(wl_um, th_ar, pol='e')
    print("nlc.n(wl_um, th_ar, pol='e') =",ans)
    print("type :", type(ans))

    nlc._cached_func_dict



main()