""" ndispers - test.py """

import numpy as np
from  . import media
# from . import __init__
pi = np.pi

def main():
    bbo = media.crystals.BetaBBO()
    print(bbo.__doc__)
    bbo.property

    pol = 'o'
    print("\npol=", pol)
    wl_um = 0.6 #um
    print("wl_um =", wl_um)
    ans = bbo.n(wl_um, 0.0, 0.0, pol)
    print("bbo.n(wl_um, 0.0, 0.0, pol='o') =", ans)
    print("type :", type(ans))

    pol = 'e'
    print("\npol=", pol)
    wl_um = 0.6 #um
    print("wl_um =", wl_um)
    ans = bbo.n(wl_um, 0.0, 0.0, pol)
    print("bbo.n(wl_um, 0.0, 0.0, pol='e') =", ans)
    print("type :", type(ans))

    wl_ar = np.arange(0.18, 1.1, 0.1)
    print("\nwl_ar =", wl_ar)
    ans = bbo.n(wl_ar, 0.0, 0.0, pol='o')
    print("bbo.n(wl_ar, 0.0, 0.0, pol='o') =", ans)
    print("type :", type(ans))

    print("\nwl_ar =", wl_ar)
    ans = bbo.n(wl_ar, 0.0, 0.0, pol='e')
    print("bbo.n(wl_ar, 0.0, 0.0, pol='e') =", ans)
    print("type :", type(ans))


    th_ar = np.arange(0, pi, 0.1*pi)
    print("\nth_ar =", th_ar)
    ans = bbo.n(wl_um, th_ar, 0.0, pol='o')
    print("bbo.n(wl_um, th_ar, 0.0, pol='o') =",ans)
    print("type :", type(ans))

    print("\nth_ar =", th_ar)
    ans = bbo.n(wl_um, th_ar, 0.0, pol='e')
    print("bbo.n(wl_um, th_ar, 0.0, pol='e') =", ans)
    print("type :", type(ans))

    ph_ar = np.arange(0, pi, 0.1*pi)
    print("\nph_ar =", ph_ar)
    ans = bbo.n(wl_um, 0.0, ph_ar, pol='o')
    print("bbo.n(wl_um, 0.0, ph_ar, pol='o') =", ans)
    print("type :", type(ans))

    print("\nph_ar =", ph_ar)
    ans = bbo.n(wl_um, 0.0, ph_ar, pol='e')
    print("bbo.n(wl_um, 0.0, ph_ar, pol='e') =", ans)
    print("type :", type(ans))

    bbo._cached_func_dict



main()