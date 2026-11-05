"""Turns the spectra identified on the Neumann records into the plant's dynamic truth, without using any candidate's form.

Taken from the data: the shape of the response across .5-20 Hz and across the 1-100 s relaxation.
Not taken from the data: the phase at the probe frequency, which is a stratum set by evidence on other tissues, because spinal cord values are not abdominal values.
A shape K(f) goes to a phase through K(f)^p, which multiplies the phase and the slope by p and keeps the shape. Above 20 Hz nothing is measured here, so the slope at 20 Hz simply continues.

neural_white, neural_grey   the two measured shapes
viscous                     flat at low frequency with a rate-like rise, as Boyer reports for skin at 10-60 Hz; a little slow memory
Output: dynamic_truth.json, branch weights per shape and integer phase, weights relative to the equilibrium force.
"""
import json
from pathlib import Path

import numpy as np
from scipy.optimize import nnls

HERE = Path(__file__).resolve().parent
TAU = np.array([.0015, .006, .025, .1, .4, 1.6, 6.4, 25.])       # factor 4, on purpose not the candidates' grid
GRID = np.geomspace(1.e-3, 1.e3, 19)
F = np.geomspace(.002, 200., 120)
F20 = 2.*np.pi*20.8


def branches(f, tau):
    s = 2j*np.pi*np.asarray(f)[..., None]*tau
    return s/(1.+s)


def measured_shape(w):
    """The identified spectrum up to 20 Hz, then a constant slope with its own phase, as causality requires."""
    e = w[0]; a = w[1:20]
    K = e+branches(F, GRID)@a
    k20 = e+branches(20., GRID)@a
    slope = np.angle(k20)/(np.pi/2.)
    high = F > 20.
    K[high] = k20*(F[high]/20.)**slope
    return K


def fit(K, tau=TAU):
    """Nonnegative branches on the plant's grid, relative error."""
    H = branches(F, tau); wt = 1./np.abs(K)
    A = np.c_[np.ones_like(F), H]
    sol, _ = nnls(np.r_[(A.real.T*wt).T, (A.imag.T*wt).T], np.r_[K.real*wt, K.imag*wt], maxiter=20000)
    return sol[0], sol[1:]


def phase(e, a, f=20.8, tau=TAU):
    return float(np.degrees(np.angle(e+branches(f, tau)@a)))


def at_phase(K, target):
    """The power of the shape that lands on the target phase at the probe frequency once on the plant's grid."""
    lo, hi = .02, 2.5
    for _ in range(50):
        p = .5*(lo+hi); e, a = fit(np.abs(K)**p*np.exp(1j*p*np.angle(K)))
        if phase(e, a) < target:
            lo = p
        else:
            hi = p
    return p, e, a


def viscous(K, target):
    """A quarter of the phase from slow memory on the measured shape, the rest from a branch so fast that it acts as a rate term up to the probe frequency."""
    _, e, a = at_phase(K, .25*target)
    a = a.copy(); a[TAU < .1] = 0.
    lo, hi = 0., 200.
    for _ in range(60):
        m = .5*(lo+hi); b = a.copy(); b[0] = m
        if phase(e, b) < target:
            lo = m
        else:
            hi = m
    return e, b


def main():
    fits = json.loads((HERE/'results'/'neumann_fits.json').read_text(encoding='utf-8'))
    pick = {r['tissue']: np.array(r['w']) for r in fits if r['family'] == 'free spectrum' and len(r['train']) == 4}
    out = dict(tau=TAU.tolist(), source='results/neumann_fits.json, free spectrum on all four protocols', shapes={})
    for name, tissue in (('neural_white', 'white'), ('neural_grey', 'grey')):
        K = measured_shape(pick[tissue])
        table = {}
        for target in range(3, 41):
            p, e, a = at_phase(K, float(target))
            table[str(target)] = dict(alpha=(a/e).tolist(), p=p, gap=phase(e, a)-target)
        out['shapes'][name] = table
    K = measured_shape(pick['white'])
    table = {}
    for target in range(3, 41):
        e, b = viscous(K, float(target))
        table[str(target)] = dict(alpha=(b/e).tolist(), p=0., gap=phase(e, b)-target)
    out['shapes']['viscous'] = table
    (HERE/'dynamic_truth.json').write_text(json.dumps(out, indent=1), encoding='utf-8')
    for name, table in out['shapes'].items():
        for target in ('8', '14', '21', '28'):
            x = table[target]; a = np.array(x['alpha'])
            m = lambda f: 1.+branches(f, TAU)@a
            print(f"{name:13s} {target:>2s} deg p {x['p']:.2f} gap {x['gap']:+.2f} sum alpha {a.sum():6.2f} | magnitude at .03, .3, 3, 20.8, 100 Hz relative to .3 Hz: "+' '.join(f"{abs(m(f))/abs(m(.3)):.2f}" for f in (.03, .3, 3., 20.8, 100.))+' | phase '+' '.join(f"{np.degrees(np.angle(m(f))):.1f}" for f in (.3, 3., 20.8, 100.)))


if __name__ == '__main__':
    main()
