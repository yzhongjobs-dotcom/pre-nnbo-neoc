"""Noise-free twin: replay logged commands with zero sensor noise (diagnosis only)."""
import numpy as np
import common
from common import P


def twin(p, log, T):
    pl = P.PLATFORMS['PC125']; dt = 5e-4
    p = np.asarray(p, float)
    NC = int(round(pl['Tc']/dt)); nt = len(log)
    out = np.zeros((nt+2, P.LOGC)); s = P.init_state(p, dt)
    zeros = np.zeros((int(T*1000)+10, 2)); extra = np.zeros((int(T*1000)+10, 4))
    args = (NC, 1, .019, 30., pl['VMAX'], dt, zeros)
    n = 0
    while n < nt:
        F_seen = s[6]
        g = np.array([1., 100., 0., 0., 0., log[n, 11]])
        n0 = n
        n = P.advance(s, p, g, 1, *args, out, n, True, extra)
        if n < 0:
            return None
        out[n0, 5] = F_seen
    return out[:n]
