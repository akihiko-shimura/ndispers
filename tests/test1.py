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

    wl_um = 0.6 #um
    T_degC = 20
    print("wl_um =", wl_um)
    print("T_degC =", T_degC)
    ans = nlc.n(wl_um, 0, T_degC, pol='o')
    print("nlc.n(wl_um, 0, T_degC, pol='o) =", ans)
    print("type :", type(ans))

    wl_um = 0.6 #um
    print("wl_um =", wl_um)
    ans = nlc.n(wl_um, 0, T_degC, pol='e')
    print("nlc.n(wl_um, 0, T_degC, pol='e) =", ans)
    print("type :", type(ans))

    wl_ar = np.arange(0.18, 1.1, 0.1)
    print("\nwl_ar =", wl_ar)
    ans = nlc.n(wl_um, 0, T_degC, pol='o')
    print("nlc.n(wl_um, 0, T_degC, pol='o') =", ans)
    print("type :", type(ans))

    print("\nwl_ar =", wl_ar)
    ans = nlc.n(wl_ar, 0, T_degC, pol='e')
    print("nlc.n(wl_ar, 0, T_degC, pol='e') =", ans)
    print("type :", type(ans))


    th_ar = np.arange(0, pi, 0.1*pi)
    print("\nth_ar =", th_ar)
    ans = nlc.n(wl_um, th_ar, T_degC, pol='o')
    print("nlc.n(wl_um, th_ar, T_degC, pol='o') =",ans)
    print("type :", type(ans))

    print("\nth_ar =", th_ar)
    ans = nlc.n(wl_um, th_ar, T_degC, pol='e')
    print("nlc.n(wl_um, th_ar, T_degC, pol='e') =", ans)
    print("type :", type(ans))

    ph_ar = np.arange(0, pi, 0.1*pi)
    print("\nph_ar =", ph_ar)
    ans = nlc.n(wl_um, 0, T_degC, pol='o')
    print("nlc.n(wl_um, 0, T_degC, pol='o') =", ans)
    print("type :", type(ans))

    print("\nph_ar =", ph_ar)
    ans = nlc.n(wl_um, 0, T_degC, pol='e')
    print("nlc.n(wl_um, 0, T_degC, pol='e') =", ans)
    print("type :", type(ans))

    nlc._cached_func_dict

    # biaxial
    print("="*40)
    nlc = media.crystals.LBO1990_xy()
    nlc.help
    nlc.constants

    wl_um = 0.6 #um
    T_degC = 20
    print("wl_um =", wl_um)
    print("T_degC =", T_degC)
    ans = nlc.n(wl_um, 0, T_degC, pol='o')
    print("nlc.n(wl_um, 0, T_degC, pol='o') =", ans)
    print("type :", type(ans))

    ans = nlc.n(wl_um, 0, T_degC, pol='e')
    print("nlc.n(wl_um, 0, T_degC, pol='e') =", ans)
    print("type :", type(ans))

    wl_ar = np.arange(0.18, 1.1, 0.1)
    print("\nwl_ar =", wl_ar)
    ans = nlc.n(wl_ar, 0, T_degC, pol='o')
    print("nlc.n(wl_ar, 0, T_degC, pol='o') =", ans)
    print("type :", type(ans))

    print("\nwl_ar =", wl_ar)
    ans = nlc.n(wl_ar, 0, T_degC, pol='e')
    print("nlc.n(wl_ar, 0, T_degC, pol='e') =", ans)
    print("type :", type(ans))

    th_ar = np.arange(0, pi, 0.1*pi)
    print("\nth_ar or phi_ar =", th_ar)
    ans = nlc.n(wl_um, th_ar, T_degC, pol='o')
    print("nlc.n(wl_um, th_ar, T_degC, pol='o') =",ans)
    print("type :", type(ans))

    th_ar = np.arange(0, pi, 0.1*pi)
    print("\nth_ar or phi_ar =", th_ar)
    ans = nlc.n(wl_um, th_ar, T_degC, pol='e')
    print("nlc.n(wl_um, th_ar, T_degC, pol='e') =",ans)
    print("type :", type(ans))

    nlc._cached_func_dict



main()