"""Compliance-domain two-channel filter with a quasi-periodic background model.

C(t) = beta*y(t) + a(t),  a = -beta*b.  The background is breathing (narrow band around f0) plus slow relief and drift:
    a = s + c1,  s: constant-velocity with white acceleration (density qs, relative to C^2),
                 (c1, c2): oscillator at f0 with white drive (density qo, relative to C^2).
Linear Kalman filter, state [s, sd, c1, c2, beta]. With y exact, C's own fast changes pass undistorted (complementary structure):
the error is the probe noise inside the band the background model leaves open, plus background outside the model.
"""
import math

import numpy as np

TC = .008


def fused_osc(data, pars, lag=12, f0=.25):
    t, ys, fs, m, R, valid, finite = data
    qs, qo, qbeta = pars
    dt = TC; tau = lag*TC; w = 2*math.pi*f0
    A = np.eye(5); A[0, 1] = dt
    cw, sw = math.cos(w*dt), math.sin(w*dt)
    A[2, 2], A[2, 3], A[3, 2], A[3, 3] = cw, sw, -sw, cw
    ct, st = math.cos(w*tau), math.sin(w*tau)
    Js = np.array([[dt**3/3, dt**2/2], [dt**2/2, dt]])
    n = len(t)
    out = np.full(n, np.nan); sd = np.full(n, np.nan); gout = np.full(n, np.nan)
    x = np.zeros(5); P = np.zeros((5, 5)); active = False; last_raw = -1000; y0 = 0.
    for j in range(n):
        if finite[j]:
            last_raw = j
        if fs[j] < 1. or j-last_raw > 25 or not np.isfinite(fs[j]+ys[j]):
            active = False
            if j:
                out[j], sd[j], gout[j] = out[j-1], sd[j-1], gout[j-1]
            continue
        if not active:
            if not valid[j]:
                if j:
                    out[j], sd[j], gout[j] = out[j-1], sd[j-1], gout[j-1]
                continue
            C0 = fs[j]/math.exp(m[j]); y0 = ys[j]; beta0 = 1./3.
            x = np.array([C0, 0., 0., 0., beta0])
            P = np.diag([C0*C0*R[j]/24, (C0*.5)**2, (C0*.5)**2, (C0*.5)**2, (beta0*1.1)**2])
            active = True
        else:
            c2 = max(x[0]+x[2]+x[4]*(ys[j]-y0), 1e-6)**2
            Q = np.zeros((5, 5)); Q[:2, :2] = qs*c2*Js; Q[2, 2] = Q[3, 3] = qo*c2*dt; Q[4, 4] = qbeta*x[4]*x[4]
            x = A@x; P = A@P@A.T+Q
            if valid[j] and j >= lag and fs[j-lag] >= 1.:
                Cobs = fs[j-lag]/math.exp(m[j])
                H = np.array([1., -tau, ct, -st, ys[j-lag]-y0])
                Rc = Cobs*Cobs*R[j]
                PH = P@H; S = float(H@PH)+Rc; K = PH/S
                x = x+K*(Cobs-H@x); I = np.eye(5)-np.outer(K, H); P = I@P@I.T+Rc*np.outer(K, K)
        H0 = np.array([1., 0., 1., 0., ys[j]-y0])
        C = max(float(H0@x), 1e-7)
        out[j] = math.log(fs[j])-math.log(C); sd[j] = math.sqrt(max(float(H0@P@H0), 0.))/C; gout[j] = 1./x[4] if x[4] > 0 else np.nan
    return out, sd, gout


def synth(seed=1, g=4., breath=.004, relief=0., dur=24., fb=.25, jitter=0.):
    """Window-like probe output: average of k over 24 ticks, noise sd .17 correlated over the window."""
    t = np.arange(0., dur, TC); ph = np.mod(t, 2.4)
    lam = 7.6+.7*np.clip((ph-.8)/.4, 0, 1)-.7*np.clip((ph-2.)/.4, 0, 1)
    ramp = ((ph >= .8) & (ph < 1.2)) | (ph >= 2.)
    rng = np.random.default_rng(seed)
    # breathing with slowly wandering frequency and amplitude (jitter = relative sd)
    wf = fb*(1+jitter*np.cumsum(rng.normal(size=len(t)))*math.sqrt(TC/4.)); amp = breath*(1+jitter*np.cumsum(rng.normal(size=len(t)))*math.sqrt(TC/4.))
    b = amp*np.sin(2*np.pi*np.cumsum(wf)*TC)+relief*np.sin(2*np.pi*.7*t+1.)
    y = b+g*8/np.exp(lam)
    wn = rng.normal(size=len(t)+24); noise = np.convolve(wn, np.ones(24)/24, mode='valid')[:len(t)]*math.sqrt(24)*.17
    k = np.exp(lam); kwin = np.convolve(np.r_[np.full(23, k[0]), k], np.ones(24)/24, mode='valid')
    ok = np.ones(len(t), bool)
    return (t, y, np.full(len(t), 8.), np.log(kwin)+noise, np.full(len(t), 24*.17**2), ok, ok), lam, ramp


if __name__ == '__main__':
    import itertools
    import dual_channel_replay as D
    for label, kw, f0 in (('exact period', {}, .25), ('filter period 20% off', {}, .30), ('jitter 15%', dict(jitter=.15), .25),
                          ('relief .5 mm at .7 Hz', dict(relief=.0005), .25), ('g 1.5', dict(g=1.5), .25), ('jitter 15% + relief + 20% off', dict(jitter=.15, relief=.0005), .30)):
        res = {}
        for seed in range(4):
            data, lam, ramp = synth(seed=seed, **kw); sc = data[0] > 6.
            r = lambda a, mk: float(np.mean((a[mk & sc]-lam[mk & sc])**2))
            res.setdefault('raw', []).append((r(data[3], ramp), r(data[3], ~ramp), 0, 0))
            a1 = D.probe_filter(data, (1e-6, 3e-3, .99)); res.setdefault('probe-IMM', []).append((r(a1, ramp), r(a1, ~ramp), 0, 0))
            for qs, qo in itertools.product((1e-3, 1e-1), (1e-3, 1e-1, 1e1)):
                a2, sd, gg = fused_osc(data, (qs, qo, 1e-7), f0=f0)
                res.setdefault(f'fused qs {qs:g} qo {qo:g}', []).append((r(a2, ramp), r(a2, ~ramp), gg[-1], np.median(sd[sc])))
        print(label)
        for k, v in res.items():
            v = np.array(v); print(f'   {k:26s} ramp {math.sqrt(v[:, 0].mean()):.3f} flat {math.sqrt(v[:, 1].mean()):.3f} g_end {np.round(v[:, 2], 1)} sd {v[:, 3].mean():.3f}')
