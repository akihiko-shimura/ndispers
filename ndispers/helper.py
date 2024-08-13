import numpy as np

def vars2(obj):
    return {k: getattr(obj, k) for k in obj.__slots__}

def returnShape(*args):
    return np.broadcast(*args).shape

def arg_signchange(a):
    a_sign = np.sign(a)
    if_signflip = ((np.roll(a_sign, 1) - a_sign) != 0).astype(int)
    if_signflip[0] = 0
    arg_signflip = np.where(if_signflip == 1)
    return arg_signflip

from scipy.interpolate import interp1d
def resample(x, y, N=1, kind='linear'):
    x = np.asarray(x)
    y = np.asarray(y)
    y_func = interp1d(x, y, kind=kind)
    x_new = np.linspace(x[0], x[-1], x.size * N)
    y_new = y_func(x_new)
    return x_new, y_new

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