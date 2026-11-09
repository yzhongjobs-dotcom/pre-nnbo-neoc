"""Huber-D4 against D4 on identical regression windows: does each theoretical prediction hold?

Reads results/huber_windows/*.npz (huber_windows.py). No plant, no tuning. Same (Y, z), same window, same rho;
only the prior differs unless a row says otherwise. Truth k_eff is used only for scoring.

Error budget of the D4 estimate on each window, telescoping on ln(k/k_eff):
  noise       noisy window minus the same commands with sensor noise off
  calibration noise off minus noise off with the sensing stage given the true calibration
  k in window k's change inside the window, as the cumulative sum of (k - k_end) dY passed through the same M
  rest        the remainder by subtraction: d through M, model error such as the known +.02 to +.05, anything else. Not an upper bound on what a prior can gain,
              since a different prior also changes how noise and calibration error pass through.
Writes results/huber_prior_check.json.
"""
import sys
sys.dont_write_bytecode = True
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.linalg import solveh_banded

import sweeps  # establishes the existing experiments' shared import paths
from estimators import KMAX, KMIN

HERE = Path(__file__).resolve().parent
WIN = HERE/'results/huber_windows'
TC, FC, F0 = .008, 14., 20.8
OMEGA = 2*math.pi*F0
RHO = (2*math.sin(math.pi*FC*TC))**(-8)
C_H = 1.345
NOISE = .04
SIGMA = .0231         # measured tick noise of z, median of z - z0 in every scene
ORDER = 4
COL = {c: i for i, c in enumerate([0, 1, 2, 3, 5, 11, 14, 18, 21, 26, 27])}
MIN_TICKS = 8
GROUPS = ['A', 'A_slow_thin', 'B', 'A+B', 'cough', 'combo', 'A_large']
SEGS = ['plateau', 'transition', 'event']
SCALES = {'s_sigma': SIGMA/NOISE, 's_.25': .25, 's_.1': .1}


def diff_mat(m):
    return np.diff(np.eye(m), ORDER, axis=0)


class Solver:
    """k = Y'M_w z / Y'M_w Y with M_w = I - (I + rho D'WD)^-1, IRLS as in d_robust.RobustRidge.estimate."""
    def __init__(self):
        self.cache = {}

    def bands(self, m):
        if m not in self.cache:
            D = diff_mat(m)
            basis = np.zeros((ORDER+1, m, len(D)))
            for q in range(ORDER+1):
                basis[q, :m-q] = (D[:, q:]*D[:, :m-q]).T
            self.cache[m] = D, basis.reshape((-1, len(D)))
        return self.cache[m]

    def fit(self, Y, z, noise=NOISE, c=C_H, iters=6, rho=RHO, scale=None):
        m = len(Y)
        D, basis = self.bands(m)
        rhs = np.column_stack((Y, z))
        s = noise/math.sqrt(rho) if scale is None else scale
        w = np.ones(m-ORDER)
        for it in range(iters):
            ab = rho*(basis@w).reshape((ORDER+1, m)); ab[0] += 1.
            smooth = solveh_banded(ab, rhs, lower=True, check_finite=False)
            sy, sz = smooth[:, 0], smooth[:, 1]
            my = Y-sy; info = float(Y@my)
            k = float(my@z)/info if info > 0 else np.nan
            if it == iters-1:
                break
            u = np.abs(D@(sz-k*sy))/s
            w = np.minimum(1., c/np.maximum(u, 1e-12))
        return dict(k=k, info=info, w=w, my=my, sz=sz, sy=sy, D=D)


def clip(k):
    return float(np.clip(k, KMIN, KMAX)) if np.isfinite(k) else np.nan


def load(path):
    with np.load(path) as f:
        data = {k: f[k] for k in f.files}
    data['meta'] = json.loads(str(data['meta']))
    return data


def masks_for(data):
    """Table 2's segments, from its own function, on the source log columns."""
    import memory_segment_eval as SE
    log = np.full((len(data['log']), SE.P3.LOGC), np.nan)
    for c, i in COL.items():
        log[:, c] = data['log'][:, i]
    _, masks = SE.segment_masks(log, dict(event_info=data['meta']['event_info']))
    return log, masks


def tail(a, m):
    x = a[-m:]
    return x[np.isfinite(x)]


def analyse(path):
    data = load(path); meta = data['meta']
    log, masks = masks_for(data)
    n = len(log)
    keff = log[:, 26]
    S = Solver()
    keys = ['D24', 'H24', 'D32', 'H32', 'H24_conv', 'D24_0', 'H24_0', 'D24_i', 'info_D24', 'info_H24', 'frac_act', 'min_w', 'noise_rms', 'e_k', 'u0', 'dk_win', 'a_nyq', 'a_f0', 'a_band']+list(SCALES)+[f'{k}_act' for k in SCALES]
    out = {k: np.full(n, np.nan) for k in keys}
    s0 = NOISE/math.sqrt(RHO)
    for j in range(n):
        Yf, zf = data['Y'][j], data['z'][j]
        if not np.isfinite(Yf).any():
            continue
        for nw, dk, hk in ((24, 'D24', 'H24'), (32, 'D32', 'H32')):
            Y, z = tail(Yf, nw), tail(zf, nw)
            if len(Y) < MIN_TICKS:
                continue
            d = S.fit(Y, z, iters=1); h = S.fit(Y, z)
            out[dk][j] = clip(d['k']); out[hk][j] = clip(h['k'])
            if nw != 24:
                continue
            out['info_D24'][j] = d['info']; out['info_H24'][j] = h['info']
            w = h['w']
            out['frac_act'][j] = float(np.mean(w < 1.)); out['min_w'][j] = float(w.min())
            # Force amplitude of the band-passed motion at the tick Nyquist (62.5 Hz) and at the probe, k_eff times M Y.
            tt = np.arange(len(Y))
            out['a_band'][j] = keff[j]*math.sqrt(2.)*float(np.std(d['my']))
            out['a_nyq'][j] = keff[j]*abs(float(np.mean(d['my']*(-1.)**tt)))
            out['a_f0'][j] = keff[j]*2*abs(np.mean(d['my']*np.exp(-2j*np.pi*F0*TC*tt)))
            out['H24_conv'][j] = clip(S.fit(Y, z, iters=60)['k'])
            for name, f in SCALES.items():
                q = S.fit(Y, z, scale=f*s0)
                out[name][j] = clip(q['k']); out[name+'_act'][j] = float(np.mean(q['w'] < 1.))
            Y0, z0 = tail(data['Y0'][j], nw), tail(data['z0'][j], nw)
            if len(Y0) != len(Y) or keff[j] <= 0:
                continue
            a0 = S.fit(Y0, z0, iters=1)
            out['D24_0'][j] = clip(a0['k']); out['H24_0'][j] = clip(S.fit(Y0, z0)['k'])
            out['noise_rms'][j] = float(np.std(z-z0))
            out['u0'][j] = float(np.max(np.abs(a0['D']@(a0['sz']-keff[j]*a0['sy']))))/s0
            Yi, zi = tail(data['Y0i'][j], nw), tail(data['z0i'][j], nw)
            if len(Yi) == len(Y):
                ai = S.fit(Yi, zi, iters=1)
                out['D24_i'][j] = clip(ai['k'])
                # Window entry i (of 32, right aligned) spans log rows j-32+i and j-32+i+1.
                r = np.arange(j-len(Yi), j)
                if r[0] >= 0:
                    # k_eff is a differential gain, so z - k_end Y from k's change alone is the cumulative sum of (k - k_end) dY,
                    # which, unlike a row form k_i Y_i, does not depend on where Y's origin lies. M removes its constant.
                    kr = .5*(keff[r]+keff[r+1])
                    part = np.r_[0., np.cumsum((kr[1:]-keff[j])*np.diff(Yi))]
                    out['e_k'][j] = math.log(max(1.+float(ai['my']@part)/ai['info']/keff[j], 1e-3))
                    out['dk_win'][j] = math.log(max(keff[j], 1.)/max(keff[r[0]], 1.))
    return dict(meta=meta, out=out, masks=dict(masks), keff=keff, raw=data['raw'], t=log[:, 0])


def robust_sd(x):
    q = np.percentile(x, [25, 75])
    return float((q[1]-q[0])/1.349)


def stats(e):
    e = e[np.isfinite(e)]
    if not len(e):
        return None
    return dict(n=len(e), med_abs=float(np.median(np.abs(e))), p90_abs=float(np.percentile(np.abs(e), 90)), bias=float(np.median(e)), rsd=robust_sd(e), below=float(np.mean(e < math.log(.85))))


def lne(r, key, m):
    return np.log(r['out'][key][m]/r['keff'][m])


def table(res):
    tab = {}
    names = ('D24', 'H24', 'D32', 'H32', 'H24_conv', 'D24_0', 'H24_0')+tuple(SCALES)
    for gr in GROUPS+['all']:
        for sg in SEGS:
            cat = {k: [] for k in names+('frac_act', 'min_w', 'inforatio', 'noise_rms', 'diff', 'act_s_sigma', 'act_s_.25', 'act_s_.1')}
            for r in res:
                if gr != 'all' and r['meta']['group'] != gr:
                    continue
                m = r['masks'][sg] & np.isfinite(r['out']['H24']) & (r['keff'] > 0)
                if not m.any():
                    continue
                for key in names:
                    cat[key].append(lne(r, key, m))
                cat['diff'].append(np.log(r['out']['H24'][m]/r['out']['D24'][m]))
                for key in ('frac_act', 'min_w', 'noise_rms'):
                    cat[key].append(r['out'][key][m])
                for key in SCALES:
                    cat['act_'+key].append(r['out'][key+'_act'][m])
                cat['inforatio'].append(r['out']['info_H24'][m]/r['out']['info_D24'][m])
            if not cat['D24']:
                continue
            c = {k: np.concatenate(v) for k, v in cat.items()}
            act = c['frac_act'] > 0
            row = {k: stats(c[k]) for k in names}
            row.update(ticks=int(len(act)), active_windows=float(np.mean(act)), rows_active=float(np.nanmean(c['frac_act'])),
                       median_min_w_active=float(np.nanmedian(c['min_w'][act])) if act.any() else None,
                       info_ratio_active=float(np.nanmedian(c['inforatio'][act])) if act.any() else None,
                       noise_rms=float(np.nanmedian(c['noise_rms'])), max_abs_lnHD=float(np.nanmax(np.abs(c['diff']))),
                       share_identical=float(np.mean(np.abs(c['diff']) < 1e-12)),
                       active_only=dict(D=stats(c['D24'][act]), H=stats(c['H24'][act]), improved=float(np.mean(np.abs(c['H24'][act]) < np.abs(c['D24'][act]))) if act.any() else None),
                       rows_active_scaled={k: float(np.nanmean(c['act_'+k])) for k in SCALES})
            tab[f'{gr}|{sg}'] = row
    return tab


def budget(res):
    """Telescoping D4 error budget; rms of each signed part."""
    out = {}
    for gr in GROUPS+['all']:
        for sg in SEGS:
            parts = {k: [] for k in ('total', 'noise', 'cal', 'k', 'rest')}
            for r in res:
                if gr != 'all' and r['meta']['group'] != gr:
                    continue
                m = r['masks'][sg] & np.isfinite(r['out']['D24']) & np.isfinite(r['out']['D24_0']) & np.isfinite(r['out']['D24_i']) & np.isfinite(r['out']['e_k']) & (r['keff'] > 0)
                if not m.any():
                    continue
                k0, ki, kd, ek = r['out']['D24_0'][m], r['out']['D24_i'][m], r['out']['D24'][m], r['out']['e_k'][m]
                parts['total'].append(np.log(kd/r['keff'][m])); parts['noise'].append(np.log(kd/k0)); parts['cal'].append(np.log(k0/ki))
                parts['k'].append(ek); parts['rest'].append(np.log(ki/r['keff'][m])-ek)
            if parts['total']:
                c = {k: np.concatenate(v) for k, v in parts.items()}
                out[f'{gr}|{sg}'] = dict(ticks=len(c['total']), **{k: dict(rms=float(np.sqrt(np.mean(v*v))), med_abs=float(np.median(np.abs(v))), p90_abs=float(np.percentile(np.abs(v), 90)), bias=float(np.median(v))) for k, v in c.items()})
    return out


def calibration(res):
    """Predicted mass-calibration error: the sensing stage removes m_hat a, the sensor holds m a, so k_est - k = (m_hat - m) w^2 on the probe."""
    rows = []
    for r in res:
        m = r['masks']['plateau'] & np.isfinite(r['out']['D24_0']) & np.isfinite(r['out']['D24_i']) & (r['keff'] > 0)
        if m.sum() < 5:
            continue
        meas = np.median(np.log(r['out']['D24_0'][m]/r['out']['D24_i'][m]))
        k = float(np.median(r['keff'][m])); dm = r['meta']['tool_mass_hat']-r['meta']['tool_mass']
        rows.append(dict(index=r['meta']['index'], group=r['meta']['group'], k=k, dm=dm, pred=math.log(max(1.+dm*OMEGA**2/k, 1e-3)), measured=float(meas)))
    return rows


def events(res):
    rows = []
    for r in res:
        t = r['t']
        for ev in r['meta']['event_info']:
            m = (t >= ev['start']) & (t <= ev['end']+.2) & np.isfinite(r['out']['H24']) & (r['keff'] > 0)
            if not m.any():
                continue
            k = float(np.median(r['keff'][m]))
            rise = (ev['end']-ev['start']-.15)/2 if ev['kind'] in ('fast', 'slow') else ev['end']-ev['start']
            rows.append(dict(group=r['meta']['group'], index=r['meta']['index'], kind=ev['kind'], amp_mm=ev['amplitude']*1e3, rise_ms=rise*1e3, k=k, d_N=k*ev['amplitude']*.778,
                             u0_max=float(np.nanmax(r['out']['u0'][m])), windows=int(m.sum()), active=float(np.mean(r['out']['frac_act'][m] > 0)), min_w=float(np.nanmin(r['out']['min_w'][m])),
                             D=float(np.median(np.abs(lne(r, 'D24', m)))), H=float(np.median(np.abs(lne(r, 'H24', m)))), D90=float(np.percentile(np.abs(lne(r, 'D24', m)), 90)), H90=float(np.percentile(np.abs(lne(r, 'H24', m)), 90)),
                             D0=float(np.nanmedian(np.abs(lne(r, 'D24_0', m)))), H0=float(np.nanmedian(np.abs(lne(r, 'H24_0', m))))))
    return rows


def absorbed(res):
    """The stated cost: does Huber let d take stiffness signal? Active windows outside events, by the sign of k's change in the window."""
    out = {}
    for name, test in (('rising', lambda d: d > .1), ('falling', lambda d: d < -.1), ('flat', lambda d: np.abs(d) <= .1)):
        dd, hh, ii = [], [], []
        for r in res:
            m = (r['out']['frac_act'] > 0) & np.isfinite(r['out']['dk_win']) & (r['keff'] > 0) & test(np.nan_to_num(r['out']['dk_win'])) & np.isfinite(r['out']['H24'])
            if m.any():
                dd.append(lne(r, 'D24', m)); hh.append(lne(r, 'H24', m)); ii.append(r['out']['info_H24'][m]/r['out']['info_D24'][m])
        if dd:
            d, h, i = map(np.concatenate, (dd, hh, ii))
            out[name] = dict(n=len(d), D_bias=float(np.median(d)), H_bias=float(np.median(h)), D_abs=float(np.median(np.abs(d))), H_abs=float(np.median(np.abs(h))), H_P90=float(np.percentile(np.abs(h), 90)), D_P90=float(np.percentile(np.abs(d), 90)), under_D=float(np.mean(d < math.log(.85))), under_H=float(np.mean(h < math.log(.85))), info_ratio=float(np.median(i)))
    return out


def ringing(res):
    """Do the windows where Huber fails share the mechanism seen one by one: a Nyquist component comparable with the probe?"""
    rows = []
    for r in res:
        m = (r['out']['frac_act'] > 0) & np.isfinite(r['out']['H24']) & (r['keff'] > 0) & np.isfinite(r['out']['a_nyq'])
        for j in np.flatnonzero(m):
            ed = math.log(r['out']['D24'][j]/r['keff'][j]); eh = math.log(r['out']['H24'][j]/r['keff'][j])
            rows.append((abs(eh)-abs(ed), r['out']['a_nyq'][j]/max(r['out']['a_f0'][j], 1e-9), r['out']['a_nyq'][j], r['out']['info_H24'][j]/r['out']['info_D24'][j], r['out']['frac_act'][j], r['keff'][j], abs(ed), abs(eh)))
    a = np.array(rows)
    out = dict(active=len(a))
    bins = ((0., .1), (.1, .3), (.3, 1.), (1., 3.), (3., 1e9))
    out['by_nyq_force_N'] = []
    for lo, hi in bins:
        m = (a[:, 2] >= lo) & (a[:, 2] < hi)
        if m.any():
            out['by_nyq_force_N'].append(dict(lo=lo, hi=hi, n=int(m.sum()), k=float(np.median(a[m, 5])), D=float(np.median(a[m, 6])), H=float(np.median(a[m, 7])), H_P90=float(np.percentile(a[m, 7], 90)), D_P90=float(np.percentile(a[m, 6], 90)), info=float(np.median(a[m, 3])), rows=float(np.median(a[m, 4])), worse_by_02=float(np.mean(a[m, 0] > .2)), better_by_02=float(np.mean(a[m, 0] < -.2))))
    worse = a[:, 0] > .2; better = a[:, 0] < -.2
    out['groups_worse'] = {}
    conv = []
    for r in res:
        m = (r['out']['frac_act'] > 0) & np.isfinite(r['out']['H24']) & (r['keff'] > 0) & np.isfinite(r['out']['a_nyq'])
        idx = np.flatnonzero(m)
        bad = idx[np.abs(np.log(r['out']['H24'][idx]/r['keff'][idx]))-np.abs(np.log(r['out']['D24'][idx]/r['keff'][idx])) > .2]
        if len(bad):
            out['groups_worse'][r['meta']['group']] = out['groups_worse'].get(r['meta']['group'], 0)+len(bad)
            conv.append(r['out']['H24_conv'][bad])
    conv = np.concatenate(conv)
    out['worse_converged_at_floor'] = float(np.mean(conv <= KMIN*(1.+1e-9))); out['worse_converged_n'] = len(conv)
    out['worse'] = dict(n=int(worse.sum()), nyq_force_med=float(np.median(a[worse, 2])), k_med=float(np.median(a[worse, 5])), info_med=float(np.median(a[worse, 3])), rows_med=float(np.median(a[worse, 4])))
    out['better'] = dict(n=int(better.sum()), nyq_force_med=float(np.median(a[better, 2])), k_med=float(np.median(a[better, 5])), info_med=float(np.median(a[better, 3])), rows_med=float(np.median(a[better, 4])))
    return out


def dense_criterion(res):
    """Stiffness misfit dk leaves dk Y_band in z - kY. Through D S (gain G at f0) it reaches the Huber threshold on every row when
    Q = |dk/k| a_band / (c s / G) >= 1, a_band the force amplitude of k M Y. Prediction: below 1 Huber follows D4; above 1 the prior relaxes on all rows,
    the probe is absorbed and the fit runs away. Q is computed from the D4 misfit, on every window, active or not."""
    x = (2*math.sin(math.pi*F0*TC))**4
    gain = x/(1+RHO*x*x)
    limit = C_H*NOISE/math.sqrt(RHO)/gain
    rows = []
    for r in res:
        m = np.isfinite(r['out']['H24']) & (r['keff'] > 0) & np.isfinite(r['out']['a_band'])
        d = r['out']['D24'][m]/r['keff'][m]
        rows.append(np.c_[np.abs(d-1.)*r['out']['a_band'][m]/limit, np.abs(np.log(d)), np.abs(np.log(r['out']['H24'][m]/r['keff'][m])), r['out']['info_H24'][m]/r['out']['info_D24'][m], r['out']['frac_act'][m], r['out']['a_band'][m], r['keff'][m]])
    a = np.concatenate(rows)
    out = dict(gain_f0=gain, limit_N=limit, bins=[])
    for lo, hi in ((0., .25), (.25, .5), (.5, 1.), (1., 2.), (2., 4.), (4., 1e9)):
        m = (a[:, 0] >= lo) & (a[:, 0] < hi)
        if m.any():
            out['bins'].append(dict(lo=lo, hi=hi, n=int(m.sum()), active=float(np.mean(a[m, 4] > 0)), rows=float(np.mean(a[m, 4])), info=float(np.median(a[m, 3])), D=float(np.median(a[m, 1])), H=float(np.median(a[m, 2])), D_P90=float(np.percentile(a[m, 1], 90)), H_P90=float(np.percentile(a[m, 2], 90)), worse_by_02=float(np.mean(a[m, 2]-a[m, 1] > .2)), band_N=float(np.median(a[m, 5])), k=float(np.median(a[m, 6]))))
    return out


def main():
    files = sorted(WIN.glob('case_*.npz'))
    with ProcessPoolExecutor(8) as ex:
        res = list(ex.map(analyse, files))
    # Regression: H24 must equal the estimator output the closed loop produced.
    worst = 0.
    for r in res:
        a, b = r['out']['H24'], r['raw']
        ok = np.isfinite(a) & np.isfinite(b)
        worst = max(worst, float(np.max(np.abs(a[ok]/b[ok]-1.))))
    report = dict(regression_H24_raw=worst, rho=RHO, s=NOISE/math.sqrt(RHO), threshold_N=C_H*NOISE/math.sqrt(RHO), table=table(res), budget=budget(res), calibration=calibration(res), events=events(res), absorbed=absorbed(res), ringing=ringing(res), dense=dense_criterion(res))
    (HERE/'results/huber_prior_check.json').write_text(json.dumps(report, indent=1), encoding='utf-8')
    print('H24 vs closed-loop raw, max rel diff', worst, 'threshold N', report['threshold_N'])
    print('--- same window: active windows, rows, identical share; |e| med/P90 D24 H24 D32; active-only D H improved; info ratio active; min w; rows active at s scales')
    for key, v in report['table'].items():
        a = v['active_only']
        print(f"{key:20s}{v['ticks']:6d} act {v['active_windows']:.3f} rows {v['rows_active']:.3f} same {v['share_identical']:.3f} | {v['D24']['med_abs']:.3f}/{v['D24']['p90_abs']:.3f} {v['H24']['med_abs']:.3f}/{v['H24']['p90_abs']:.3f} {v['D32']['med_abs']:.3f}/{v['D32']['p90_abs']:.3f} | "
              +(f"act {a['D']['med_abs']:.3f}/{a['D']['p90_abs']:.3f} {a['H']['med_abs']:.3f}/{a['H']['p90_abs']:.3f} imp {a['improved']:.2f} info {v['info_ratio_active']:.3f} minw {v['median_min_w_active']:.2f}" if a['D'] else 'no active')
              +' | '+' '.join(f"{k} {v[k]['med_abs']:.3f}/{v[k]['p90_abs']:.3f} rows {v['rows_active_scaled'][k]:.2f}" for k in SCALES))
    print('--- D4 24 budget, rms: total noise cal k rest (rest med_abs)')
    for key, v in report['budget'].items():
        print(f"{key:20s}{v['ticks']:6d} total {v['total']['rms']:.3f} noise {v['noise']['rms']:.3f} cal {v['cal']['rms']:.3f} k {v['k']['rms']:.3f} rest {v['rest']['rms']:.3f} (med_abs {v['total']['med_abs']:.3f} {v['noise']['med_abs']:.3f} {v['cal']['med_abs']:.3f} {v['k']['med_abs']:.3f} {v['rest']['med_abs']:.3f})")
    print('--- calibration, plateau, predicted vs measured')
    for c in report['calibration']:
        if c['group'] in ('B', 'cough', 'A'):
            print(f"{c['group']:6s} {c['index']:2d} k {c['k']:7.0f} dm {c['dm']*1e3:+6.1f} g pred {c['pred']:+.3f} meas {c['measured']:+.3f}")
    print('--- events')
    for e in report['events']:
        print(f"{e['group']:6s} {e['index']:2d} {e['kind']:6s} amp {e['amp_mm']:5.2f} mm rise {e['rise_ms']:4.0f} ms k {e['k']:6.0f} d {e['d_N']:5.1f} N u0 {e['u0_max']:6.2f} act {e['active']:.2f} minw {e['min_w']:.2f} | D {e['D']:.3f}/{e['D90']:.3f} H {e['H']:.3f}/{e['H90']:.3f} | noise-free {e['D0']:.3f} {e['H0']:.3f}")
    print('--- active windows by k change in window'); print(json.dumps(report['absorbed'], indent=0))
    print('--- ringing'); print(json.dumps(report['ringing'], indent=0))
    print('--- dense criterion, limit N', report['dense']['limit_N'], 'gain', report['dense']['gain_f0'])
    for b in report['dense']['bins']:
        print(f"Q {b['lo']:.2f}-{b['hi']:.2f} n {b['n']:6d} active {b['active']:.3f} rows {b['rows']:.3f} info {b['info']:.3f} D {b['D']:.3f}/{b['D_P90']:.3f} H {b['H']:.3f}/{b['H_P90']:.3f} worse {b['worse_by_02']:.3f} band {b['band_N']:.3f} N k {b['k']:.0f}")


if __name__ == '__main__':
    main()
