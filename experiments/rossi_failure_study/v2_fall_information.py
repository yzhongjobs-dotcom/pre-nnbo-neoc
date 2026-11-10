"""Why does a falling window carry less information about its endpoint gain? Closed loop, Huber-D4 ungated, four arms:
  base; notch (the force error reaching the admittance feedback passes a causal 20.8 Hz notch, so feedback does not react to the probe);
  exc_true (excitation increment divided by the true k_eff instead of k_good, as v2_excitation_control.py); exc_true_notch (both).
Per tick, from the estimator's own preprocessed log and its final Huber weights, for the 24-row window ending now:
  I_end  = k_T^2 (Y'M_W Y)^2 / ||M_W Y||^2   information about the endpoint gain: 1/I_end is the relative variance of k_hat up to sigma^2
  N      = sum_i (k_T/k_i)^2 m_i^2 / sum m_i^2, m = M_W Y: the within-window distribution, < 1 when earlier rows had a higher gain (falls)
  I_inst = I_end / N                         what the window would carry if every row were scaled to its own gain (excitation, feedback, servo)
Reported by segment relative to the plateau median of the same case. Diagnostic counterfactuals only.
Predictions fixed before running: on falls N is about .86 and I_inst is below the plateau; the notch removes the feedback boost on rises and part
of the loss on falls; exc_true removes the k_good-lag part of I_inst on falls.
Writes results/v2_fall_information.json.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_key] = '1'
import sys
sys.dont_write_bytecode = True
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.signal import iirnotch

import corrected_comparison as CC
import memory_segment_eval as SE
import overlay_compare as O
import rib_transition_eval as R
import v2_verify as V
from d_robust import common

B, M = O.B, SE.M
HERE = Path(__file__).resolve().parent
CASES = (0, 1, 2, 36, 37, 38, 39, 40, 41)
ARMS = tuple(sys.argv[1].split(',')) if len(sys.argv) > 1 else ('base', 'notch', 'exc_true', 'exc_true_notch')
NW = 24


class Notch:
    def __init__(self):
        self.b, self.a = iirnotch(20.8, 2., fs=125.)
        self.x = [0., 0.]; self.y = [0., 0.]

    def __call__(self, e):
        b, a = self.b, self.a
        out = (b[0]*e+b[1]*self.x[0]+b[2]*self.x[1]-a[1]*self.y[0]-a[2]*self.y[1])/a[0]
        self.x = [e, self.x[0]]; self.y = [out, self.y[0]]
        return out


def job(task):
    arm, idx = task
    O.setup()
    case = CC.cases()[idx]
    sc = SE.MC.scene_for(tuple(case['cell']), CC.BASE, case['variant'])
    pol = B.build('Huber-D4', O.chosen()['Huber-D4'], sc)
    CC.own_mechanism(pol)
    patched = []
    est = pol.estimator; inner_est = est.estimate
    truth = {}; info = {}

    def estimate(o, n, p):
        out = inner_est(o, n, p)
        kt = truth.get(n-1)
        if kt is None or not np.isfinite(out[0]):
            return out
        log = est._log(o, n)
        j0 = max(p.ndel+1, n-1-NW); m = n-2-j0+1
        if m < NW:
            return out
        _, Y, z, keep = common._build(log, j0, m, 0, p.ns, p.eb, p.afc, p.ndel, log[n-1, 0], log[n-1, 2], NW*p.ns*1e-3, p.params['enc_y'])
        w = getattr(est, 'weights', None)
        if not keep.all() or w is None or len(w) != NW-4:
            return out
        D = np.diff(np.eye(NW), 4, axis=0)
        Mw = np.eye(NW)-np.linalg.inv(np.eye(NW)+est.rho*(D.T*w)@D)
        my = Mw@Y
        ks = np.array([truth.get(j0+i+1, np.nan) for i in range(NW)])
        if not np.isfinite(ks).all():
            return out
        # Endpoint information = 1 / relative variance of q'z, q = M_W Y / (Y'M_W Y); window distribution factor
        # N = sum (k_T/k_i)^2 m_i^2 / sum m_i^2 with positive weights m_i^2 (< 1 when earlier rows had a higher gain).
        Ie_ = kt*kt*float(Y@my)**2/float(my@my)
        N_ = float(np.sum((kt/ks)**2*my*my)/np.sum(my*my))
        M0 = np.eye(NW)-np.linalg.inv(np.eye(NW)+est.rho*D.T@D); m0 = M0@Y
        I0 = kt*kt*float(Y@m0)**2/float(m0@m0)                       # same with Huber weights off
        tt = np.arange(NW)*p.Tc; car = np.exp(-2j*np.pi*20.8*tt)
        Ip = kt*kt*abs(np.mean(m0*car))**2*NW                          # probe-band part only (20.8 Hz component of M Y)
        info[n-1] = (Ie_, Ie_/N_, I0, Ip, N_)
        return out
    est.estimate = estimate
    sensing = pol.step; inner = sensing.step

    def step(o):
        kt = float(getattr(o, 'truth_keff', np.nan))
        truth[o.n-1] = kt if kt > 0 else np.nan
        if 'notch' in arm and (not patched or patched[0] is not pol.adm):
            # the admittance is created at reset; patch the live object
            adm = pol.adm; inc = adm.increment; nf = Notch()
            adm.increment = lambda e, inc=inc, nf=nf: inc(nf(e))
            patched[:] = [adm]
        x_before = pol.x_prev
        u = inner(o)
        if arm.startswith('exc_true'):
            dx = pol.x_prev-x_before
            k_un = pol.k_good if np.isfinite(pol.k_good) else pol.k_seek
            if kt > 0. and np.isfinite(k_un) and k_un > 0.:
                u = u+dx*(1./kt-1./k_un)
        return u
    sensing.step = step
    log = R.oracle_run(sc, pol, case['seed'])
    m = M.metric(log, pol, sc)
    t = log[:, 0]; k = log[:, 26]
    seg = V.segments(t, k, sc['event_info'], k > 200.)
    Ie, Ii, I0, Ip, NN = (np.full(len(t), np.nan) for _ in range(5))
    for j, (a, b, c, d, e) in info.items():
        if j < len(t):
            Ie[j], Ii[j], I0[j], Ip[j], NN[j] = a, b, c, d, e
    base_e = np.nanmedian(Ie[seg['plateau']]) if np.isfinite(Ie[seg['plateau']]).sum() > 20 else np.nanmedian(Ie)
    out = dict(arm=arm, index=idx, group=case['group'], L=float(np.log1p(m['J'])), seg={})
    for s in ('plateau', 'rising', 'falling', 'reversal'):
        mm = seg[s] & np.isfinite(Ie) & np.isfinite(Ii) & (Ii > 0)
        if mm.sum() < 5:
            continue
        pl = seg['plateau'] & np.isfinite(I0)
        b0 = np.nanmedian(I0[pl]) if pl.sum() > 20 else np.nanmedian(I0); bp = np.nanmedian(Ip[pl]) if pl.sum() > 20 else np.nanmedian(Ip)
        out['seg'][s] = dict(I_end=float(np.median(Ie[mm]/base_e)), I_inst=float(np.median(Ii[mm]/base_e)), N=float(np.median(NN[mm])),
                             I_noW=float(np.median(I0[mm]/b0)), I_probe=float(np.median(Ip[mm]/bp)), I_probe_inst=float(np.median(Ip[mm]/NN[mm]/bp)), n=int(mm.sum()))
    return out


def main():
    with ProcessPoolExecutor(8) as ex:
        rows = list(ex.map(job, [(a, i) for a in ARMS for i in CASES], chunksize=1))
    rep = {}
    print('median over cases of segment medians; I relative to the plateau I_end of the same case')
    print('arm             segment   I_end  I_inst  N (within-window)   [implied sd vs plateau = I_end^-1/2]')
    for a in ARMS:
        for s in ('rising', 'falling', 'reversal'):
            v = [r['seg'][s] for r in rows if r['arm'] == a and s in r['seg']]
            if not v:
                continue
            cell = {q: float(np.nanmedian([x[q] for x in v])) for q in ('I_end', 'I_inst', 'N', 'I_noW', 'I_probe', 'I_probe_inst')}
            rep[f'{a}|{s}'] = cell
            print(f"  {a:15s} {s:9s} {cell['I_end']:.2f}   {cell['I_inst']:.2f}   {cell['N']:.2f}  | W=I {cell['I_noW']:.2f} | probe band {cell['I_probe']:.2f}, per row {cell['I_probe_inst']:.2f}")
    print('closed loop L by arm:', {a: round(float(np.mean([r['L'] for r in rows if r['arm'] == a])), 3) for a in ARMS})
    (HERE/'results/v2_fall_information.json').write_text(json.dumps(dict(summary=rep, rows=rows), indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
