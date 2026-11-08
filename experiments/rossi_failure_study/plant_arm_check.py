"""Check the dynamic arm option and measure whether the arm approximation matters for identification. Part 1 of the plant review.

Stage `frf`: stored baseline reproduces with the option off; open-loop response of the deflection against the analytic mass on spring; step-size convergence by comparing dynamic zeta values in the same scene; k_eff column against a noise-free lock-in.
Usage: python -X utf8 plant_arm_check.py frf
"""
import json
import sys
from pathlib import Path

import numpy as np

import overlay_compare as O
import memory_closed_loop as MC
import physics_plant3 as P3

B = O.B
HERE = Path(__file__).resolve().parent
BASE = 874500
CASES = [('A_neural_i6', ('breath_rib', 'rib', 'neural_white', 28., 6), 'nominal'),
         ('slow_neural_i1', ('breath_rib', 'rib', 'neural_white', 28., 1), 'slow_thin'),
         ('slow_viscous_i5', ('breath_rib', 'rib', 'viscous', 14., 5), 'slow_thin')]


def frozen(cell, variant, **extra):
    """Same contact condition held still: no scan motion, no exogenous motion, no sensor noise (replay default)."""
    sc = MC.scene_for(cell, BASE, variant)
    n = len(sc['u'])
    x4 = float(np.interp(4., np.arange(n)*.001, sc['xcmd']))
    sc['xcmd'] = np.full(n, x4); sc['u'] = np.zeros(n); sc['act'] = np.zeros(n); sc['fext'] = np.zeros(n)
    sc['alpha_t0'] = None
    sc.pop('alpha_t0')
    sc.update(extra)
    return sc


def phasor(t, v, f):
    s = t-t.mean()
    A = np.c_[np.cos(2*np.pi*f*t), np.sin(2*np.pi*f*t), np.ones(len(t)), s, s*s]
    b = np.linalg.lstsq(A, v, rcond=None)[0]
    return complex(b[0], -b[1])


def sine_run(sc, f, amp=20e-6, T=6.):
    """Command sine about the settled hold. The command's phase is kept fine with a 1 ms log."""
    pl = dict(P3.PLATFORM, log_dt=.001, tau_a=sc.get('tau_a_true', P3.PLATFORM['tau_a']))
    sc = dict(sc); sc['noise_scale'] = 0.
    s0, buf = P3._init(sc, pl)
    y0 = float(buf[-1])
    n = int(T/.008)
    tc = np.arange(n)*.008
    cmd = y0+amp/abs(1./(1.+2j*np.pi*f*pl['tau_a']))*np.sin(2*np.pi*f*tc)
    L = P3.replay(sc, cmd, seed=0, sensor=False, platform=pl)
    t = L[:, 0]; w = t > T-3.
    return t[w], L[w]


def main():
    out = {}
    for name, cell, variant in CASES:
        O.setup()
        for label, extra in (('old', {}), ('dyn0.1', dict(arm_zeta=.1)), ('dyn0.03', dict(arm_zeta=.03)), ('dyn0.3', dict(arm_zeta=.3))):
            sc = frozen(cell, variant, **extra)
            kr, m = sc['kr'], sc['tool_mass']
            rows = []
            for f in (5., 10., 20.8, 35., 60., 90., 120.):
                t, L = sine_run(sc, f)
                y, Fx, d = L[:, 2], L[:, 18], L[:, 6]
                Y = phasor(t, y, f); F = phasor(t, Fx, f); D = phasor(t, d, f)
                rows.append(dict(f=f, gain=abs(F/Y), phase=float(np.degrees(np.angle(F/Y))), dy=abs(D/Y), keff=float(np.median(L[:, 26])), K=float(np.median(L[:, 21]))))
            out[f'{name} {label}'] = dict(kr=kr, m=m, fn=float(np.sqrt(kr/m)/2/np.pi), rows=rows)
            r = {x['f']: x for x in rows}
            print(f"{name:16s} {label:8s} fn {np.sqrt(kr/m)/2/np.pi:5.1f} | F/y@20.8 {r[20.8]['gain']:8.0f} keff col {r[20.8]['keff']:8.0f} ratio {r[20.8]['gain']/r[20.8]['keff']:.4f} col21 {r[20.8]['K']:8.0f} | F/y " + ' '.join(f"{x['f']:g}:{x['gain']:.0f}" for x in rows), flush=True)
    (HERE/'results'/'plant_arm_check_frf.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
