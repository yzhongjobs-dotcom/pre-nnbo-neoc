"""Where does Huber-D4 beat D4 when the disturbance is known exactly?

One window, 24 ticks, 8 ms: z = k Y + d + e, Y the probe response only, e white at the measured tick noise.
d is a smooth-step (the scene's quintic push) or a sin^2 pulse, with the given rise and amplitude (force units in z),
or an ideal step, slid across the window. Reports, per shape:
  u_pred: largest |D S d| / s from the noise-free first iteration, the prediction of activation (active if > c);
  active: share of windows with any weight below one, noisy;
  D, H: RMS relative stiffness error; leak: the D4 error from d alone, the part a better prior could remove.
Writes results/huber_synthetic_map.json.
"""
import sys
sys.dont_write_bytecode = True
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from huber_prior_check import C_H, NOISE, RHO, Solver

HERE = Path(__file__).resolve().parent
TC, F0 = .008, 20.8
GAIN = .778           # z's gain on a force or motion held over the tick: 1 - a^ns
K = 2000.
AY = .1674*GAIN       # probe force in z
SIGMA = .0231         # measured tick noise of z
M = 24
RISES = (0., .016, .032, .064, .1, .15, .27, .4)
AMPS = (.3, 1., 3., 10., 30.)


def quintic(x):
    x = np.clip(x, 0., 1.)
    return x*x*x*(10.+x*(-15.+6.*x))


def shape(kind, rise, t, c):
    if kind == 'step':
        return (t >= c).astype(float) if rise == 0. else quintic((t-c)/rise)
    if rise == 0.:
        return (np.abs(t-c) < TC/2).astype(float)
    return np.sin(np.pi*np.clip((t-c)/rise, 0., 1.))**2


def one(args):
    kind, rise, amp = args
    S = Solver(); D = S.bands(M)[0]
    rng = np.random.default_rng([int(rise*1e3), int(amp*10), kind == 'step'])
    t = np.arange(M)*TC
    s = NOISE/math.sqrt(RHO)
    span = M*TC+rise
    centres = np.linspace(-rise, M*TC, 16)
    phases = np.linspace(0., 2*np.pi, 6, endpoint=False)
    leak = []; ed = []; eh = []; act = []; umax = 0.; hclean = []; info = []
    for c in centres:
        d = amp*shape(kind, rise, t, c)
        for ph in phases:
            Y = AY/K*np.sin(2*np.pi*F0*t+ph)
            clean = S.fit(Y, K*Y+d, iters=1)
            leak.append(clean['k']/K-1.)
            umax = max(umax, float(np.max(np.abs(D@(clean['sz']-clean['k']*clean['sy'])))/s))
            hclean.append(S.fit(Y, K*Y+d)['k']/K-1.)
            for _ in range(4):
                z = K*Y+d+SIGMA*rng.standard_normal(M)
                a = S.fit(Y, z, iters=1); b = S.fit(Y, z)
                ed.append(a['k']/K-1.); eh.append(b['k']/K-1.); act.append(bool(np.any(b['w'] < 1.)))
                info.append(b['info']/a['info'])
    rms = lambda x: float(np.sqrt(np.mean(np.square(x))))
    return dict(kind=kind, rise_ms=rise*1e3, amp=amp, u_pred=umax, active=float(np.mean(act)), D=rms(ed), H=rms(eh), leak=rms(leak), leak_H=rms(hclean), info=float(np.median(info)))


def influence_bound():
    """Large step in one window. With A the rows at psi = +-c at the solution, sgn their signs, N the rest:
    d-stationarity: (I + rho D_N'D_N) d = z - kY - rho s c D_A' sgn. k-stationarity then gives, for any d in the null space of D_N,
        dk = rho s c sgn' D_A S_N Y / (Y' M_N Y),  S_N = (I + rho D_N'D_N)^-1, M_N = I - S_N,
    which does not contain the step height. rho s = sigma sqrt(rho), so dk/k = (c sigma sqrt(rho) / (k A)) Gamma with Gamma from window geometry only.
    Checked against the solver, and for proportionality in the Huber scale s."""
    S = Solver(); D = S.bands(M)[0]; t = np.arange(M)*TC
    s0 = NOISE/math.sqrt(RHO)
    out = {}
    for factor in (1., SIGMA/NOISE, .25):
        s = factor*s0; rows = []
        for c in range(1, M):
            step = (np.arange(M) >= c).astype(float)
            for ph in np.linspace(0., 2*np.pi, 12, endpoint=False):
                Y = AY/K*np.sin(2*np.pi*F0*t+ph)
                f = S.fit(Y, K*Y+1000.*step, iters=2000, scale=s)
                d = f['sz']-f['k']*f['sy']; u = D@d/s
                act = np.abs(u) >= C_H*(1.-1e-9)
                DN, DA, sg = D[~act], D[act], np.sign(u[act])
                SN = np.linalg.inv(np.eye(M)+RHO*DN.T@DN)
                pred = RHO*s*C_H*float(sg@(DA@(SN@Y)))/float(Y@(Y-SN@Y))
                rows.append((pred/K, f['k']/K-1., int(act.sum())))
        a = np.array(rows)
        out[f'{factor:.3f}'] = dict(formula_vs_solver_max=float(np.max(np.abs(a[:, 0]-a[:, 1]))), median=float(np.median(np.abs(a[:, 1]))), p90=float(np.percentile(np.abs(a[:, 1]), 90)), max=float(np.max(np.abs(a[:, 1]))), rows_at_bound=float(np.median(a[:, 2])),
                                scale=C_H*factor*NOISE*math.sqrt(RHO)/AY)
    # Noise-only row std of u, from the linear D S map, at the Huber scale in use and at the measured tick noise.
    gain = D@np.linalg.inv(np.eye(M)+RHO*D.T@D)
    row = float(np.sqrt(np.max(np.sum(gain*gain, axis=1))))
    out['u_noise_std'] = dict(nominal_scale=row*SIGMA/s0, true_sigma_scale=row*math.sqrt(RHO))
    return out


def main():
    bound = influence_bound()
    print('influence bound', json.dumps(bound))
    tasks = [(k, r, a) for k in ('step', 'pulse') for r in RISES for a in AMPS]
    with ProcessPoolExecutor(8) as ex:
        rows = list(ex.map(one, tasks, chunksize=1))
    floor = SIGMA/math.sqrt(sum((AY*np.sin(2*np.pi*F0*np.arange(M)*TC))**2)*.77)
    print(f'noise-only floor about {floor:.3f}; threshold c = {C_H}')
    print('kind  rise_ms amp  u_pred active  D     H    | noise-free leak D, H | info H/D')
    for r in rows:
        print(f"{r['kind']:5s} {r['rise_ms']:5.0f} {r['amp']:5.1f} {r['u_pred']:8.2f} {r['active']:.2f} {r['D']:.3f} {r['H']:.3f} | {r['leak']:.3f} {r['leak_H']:.3f} | {r['info']:.3f}")
    (HERE/'results/huber_synthetic_map.json').write_text(json.dumps(dict(k=K, sigma=SIGMA, ay=AY, influence_bound=bound, rows=rows), indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
