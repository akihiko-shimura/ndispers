import numpy as np

def vars2(obj):
    slots = [s for klass in type(obj).__mro__ for s in getattr(klass, "__slots__", ())]
    return {k: getattr(obj, k) for k in slots if hasattr(obj, k)}

def returnShape(*args):
    return np.broadcast(*args).shape

def arg_signchange(a):
    a_sign = np.sign(a)
    if_signflip = ((np.roll(a_sign, 1) - a_sign) != 0).astype(int)
    if_signflip[0] = 0
    arg_signflip = np.where(if_signflip == 1)
    return arg_signflip

def resample(x, y, N=1):
    x = np.asarray(x)
    y = np.asarray(y)
    x_new = np.linspace(x[0], x[-1], x.size * N)
    return x_new, np.interp(x_new, x, y)

def fullWidth(x_ar, y_ar, threshold=0.5, N=1):
    x, y = resample(x_ar, y_ar, N=N)
    idx_3dB = np.where(y >= np.max(y) * threshold)
    x_3dB = x[idx_3dB]
    width_3dB = x_3dB[-1] - x_3dB[0]
    return width_3dB

def peak_position(x_ar, y_ar, N=1):
    x, y = resample(x_ar, y_ar, N=N)
    idx_peak = np.argmax(y_ar)
    return x_peak[idx_peak]

def brentq(f, a, b, xtol=2e-12, rtol=8.9e-16, maxiter=100):
    """Root of f in [a, b] by Brent's method (f(a) and f(b) of opposite sign).

    Same algorithm and stopping rule as scipy.optimize.brentq, kept here so the
    package (and its WASM apps) need not pull in scipy for one function.
    """
    fa, fb = f(a), f(b)
    if fa == 0.0:
        return a
    if fb == 0.0:
        return b
    if fa * fb > 0:
        raise ValueError("f(a) and f(b) must have different signs")
    c, fc = a, fa
    d = e = b - a
    for _ in range(maxiter):
        if fb * fc > 0:
            c, fc = a, fa
            d = e = b - a
        if abs(fc) < abs(fb):
            a, b, c = b, c, b
            fa, fb, fc = fb, fc, fb
        tol = 2.0 * rtol * abs(b) + 0.5 * xtol
        m = 0.5 * (c - b)
        if fb == 0.0 or abs(m) <= tol:
            return b
        if abs(e) >= tol and abs(fa) > abs(fb):
            s = fb / fa
            if a == c:                       # secant
                p, q = 2.0 * m * s, 1.0 - s
            else:                            # inverse quadratic interpolation
                q, r = fa / fc, fb / fc
                p = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0))
                q = (q - 1.0) * (r - 1.0) * (s - 1.0)
            if p > 0:
                q = -q
            else:
                p = -p
            if 2.0 * p < min(3.0 * m * q - abs(tol * q), abs(e * q)):
                e, d = d, p / q
            else:
                d = e = m
        else:
            d = e = m
        a, fa = b, fb
        b += d if abs(d) > tol else (tol if m > 0 else -tol)
        fb = f(b)
    raise RuntimeError("brentq: no convergence in %d iterations" % maxiter)
