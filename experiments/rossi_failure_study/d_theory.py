"""What each d treatment does to a window, in closed form. Both are k = Y'Mz / Y'MY, so the error from any d is Y'Md / Y'MY.

Y is the probe's own response in the window (the 20.8 Hz excitation, force units; amplitude of the Method's a_force), sampled at the 8 ms control period.
For each M: information kept (Y'MY / Y'Y, the share of the probe signal left after removing what the prior attributes to d),
and the relative stiffness error that one newton of d of a given shape causes, averaged over phase and position in the window.
Shapes: sinusoids (breathing to above the probe frequency), smooth pulses of the length of flinches and coughs, steps, and components locked to the probe.
"""
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TC, F0, A = .008, 20.8, .1674          # control period, probe frequency, probe force amplitude (Method a_force)


def d4(m, fc):
    rho = (2*math.sin(math.pi*fc*TC))**(-8)
    D = np.diff(np.eye(m), 4, axis=0)
    return np.eye(m)-np.linalg.inv(np.eye(m)+rho*D.T@D)


def poly(m, degree):
    t = np.linspace(-1, 1, m); Q, _ = np.linalg.qr(np.vander(t, degree+1))
    return np.eye(m)-Q@Q.T


ESTIMATORS = {'D4 24 ticks, fc 14 Hz (chosen)': d4(24, 14.), 'D4 24 ticks, fc 10 Hz': d4(24, 10.), 'quintic 24 ticks (chosen)': poly(24, 5),
              'quintic 16 ticks (old)': poly(16, 5), 'cubic 24 ticks': poly(24, 3)}


def probe(m, phase):
    t = np.arange(m)*TC
    return A*np.sin(2*math.pi*F0*t+phase), A*np.cos(2*math.pi*F0*t+phase)


def error(M, d_of_phase, n=24):
    """RMS relative stiffness error over probe phases, per unit of the d shape."""
    out = []
    for ph in np.linspace(0, 2*math.pi, n, endpoint=False):
        Y, _ = probe(M.shape[0], ph)
        out.append(float(Y@M@d_of_phase(ph))/float(Y@M@Y))
    return float(np.sqrt(np.mean(np.square(out))))


def main():
    res = {}
    freqs = [.25, 1., 3., 6., 10., 14., 18., 20.8, 24., 30., 40.]
    print('information kept, Y\'MY / Y\'Y, mean over phase; then relative k error per 1 N of d: sinusoids at f Hz (RMS over both phases)')
    print('  '+' '*34+'info  '+' '.join(f"{f:>6g}" for f in freqs))
    for name, M in ESTIMATORS.items():
        m = M.shape[0]; t = np.arange(m)*TC
        info = np.mean([float(probe(m, ph)[0]@M@probe(m, ph)[0])/float(probe(m, ph)[0]@probe(m, ph)[0]) for ph in np.linspace(0, 2*math.pi, 24, endpoint=False)])
        row = [np.sqrt(np.mean([error(M, lambda ph, f=f, q=q: np.sin(2*math.pi*f*t+q+ph*0)) ** 2 for q in np.linspace(0, 2*math.pi, 8, endpoint=False)])) for f in freqs]
        res[name] = dict(info=float(info), sinusoid=dict(zip(map(str, freqs), map(float, row))))
        print(f"  {name:34s} {info:.3f} "+' '.join(f"{v:6.3f}" for v in row))
    print('\nrelative k error per 1 N of d: smooth sin^2 pulse of the given length placed anywhere in or across the window (RMS), a step anywhere, a ramp of 1 N over the window')
    widths = (.05, .1, .2, .4)
    for name, M in ESTIMATORS.items():
        m = M.shape[0]; t = np.arange(m)*TC; out = {}
        for w in widths:
            vals = []
            for c in np.linspace(-w, m*TC+w, 40):
                d = np.where(np.abs(t-c) < w/2, np.cos(math.pi*(t-c)/w)**2, 0.)
                vals.append(error(M, lambda ph, d=d: d))
            out[f'pulse {w*1e3:.0f} ms'] = float(np.sqrt(np.mean(np.square(vals))))
        out['step'] = float(np.sqrt(np.mean([error(M, lambda ph, j=j: (np.arange(m) >= j).astype(float))**2 for j in range(1, m)])))
        out['ramp'] = error(M, lambda ph: np.arange(m)/(m-1))
        res[name].update(out)
        print(f"  {name:34s} "+' '.join(f"{k} {v:.3f}" for k, v in out.items()))
    print('\ncomponents locked to the probe (unit amplitude relative to the probe response): in phase gives a bias equal to its size under every M; quadrature:')
    for name, M in ESTIMATORS.items():
        m = M.shape[0]
        q = error(M, lambda ph: probe(m, ph)[1]/A)
        res[name]['quadrature'] = q
        print(f"  {name:34s} in phase 1/A each (unidentifiable), quadrature {q:.4f}")
    (HERE/'results'/'d_theory.json').write_text(json.dumps(res, indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
