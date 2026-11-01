"""Controller-side finite-footprint contact model (cylinder line contact, Hertz-type elliptic pressure).

Known: probe radius R, length L and reference width LREF (probe geometry). Unknown: the modulus field E(x) along
the known scan path and the indentation depth d. For a fixed half-width a(d), the normal force fe, the pressure
first moment m1 and the static stiffness k = d fe / d d are LINEAR in E; with E on a 1-D grid (linear
interpolation) each is a kernel row times the cell vector. This is the textbook line-contact model (the same form as
the benchmark plant's elastic part); the plant adds standard-linear-solid relaxation, Hunt-Crossley damping, sensor
filtering, noise and bias, and sharp ridge edges that the grid only approximates. Mismatch tests change the plant,
not this model.
"""
import math
import numpy as np

R, L, LREF = 0.03, 0.04, 0.03
NQ = 12
_th = np.arange(1, NQ+1)*math.pi/(NQ+1)


def psi(a):
    a = np.clip(a, 1e-9, LREF*(1-1e-9)); w = np.sqrt(LREF*LREF-a*a)
    return (LREF*LREF-LREF*w+a*a*np.log((LREF+w)/a))/(2*R)


def halfwidth(d):
    top = LREF*.95; dd = np.clip(d, 1e-10, psi(top)); a = np.clip(np.sqrt(R*dd), 1e-8, LREF*.9)
    for _ in range(12):
        w = np.sqrt(LREF*LREF-a*a); der = a/R*np.log((LREF+w)/a)
        a = np.clip(a-(psi(a)-dd)/np.maximum(der, 1e-12), 1e-8, top)
    return a


class Grid:
    def __init__(self, x0, x1, h=1e-3):
        self.x = np.arange(x0, x1+h/2, h); self.h = h; self.n = len(self.x)

    def interp_rows(self, xs):
        """linear-interpolation weights as a dense (len(xs), n) matrix"""
        u = (np.asarray(xs)-self.x[0])/self.h
        i = np.clip(np.floor(u).astype(int), 0, self.n-2); f = np.clip(u-i, 0, 1)
        W = np.zeros((len(xs), self.n)); r = np.arange(len(xs)); W[r, i] = 1-f; W[r, i+1] += f
        return W


def fe_m1_rows(grid, xc, d):
    a = float(halfwidth(d)); xi = a*np.cos(_th); w = math.pi/(NQ+1)*np.sin(_th)**2*a*a*L/(2*R)
    W = grid.interp_rows(xc+xi)
    return w@W, (xi*w)@W


def kernels(grid, xc, d, hd=2e-6):
    """rows (Fe, M1, K): fe = Fe@E, m1 = M1@E, static stiffness k = K@E at centre xc and depth d"""
    F0, M0 = fe_m1_rows(grid, xc, d)
    Fp, _ = fe_m1_rows(grid, xc, d+hd); Fm, _ = fe_m1_rows(grid, xc, max(d-hd, 1e-9))
    return F0, M0, (Fp-Fm)/(2*hd)
