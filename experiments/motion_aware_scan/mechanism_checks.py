"""Matched mechanism checks after the skill-based self-review (paper draft §五.5, §五.1).

E1 supervisor parity: every comparator wrapped in the same contact-loss supervisor as the proposed
   method (force < 0.3 N: stop calling the controller, approach at 6 mm/s; force >= 0.5 N: resume),
   with its confirmed tuned parameters.
E2 predictor swap with every other parameter fixed at the selected Proposed-IMM values:
   window-polynomial extrapolation (no tracker parameters), single-model Kalman (q_kf, r_scale tuned
   on the development split with the equal 61-evaluation budget), two-model IMM (final).
E3 normalisation removed: feedback uses a fixed stiffness (the two fixed-impedance ends) while the
   estimator, excitation, motion tracker and supervisor stay identical.
Confirmation cohort; paired against the final method's stored confirmation rows.
Run: python -X utf8 mechanism_checks.py  -> results_v2/mechanism_checks.json
"""
import json, math, os
import numpy as np
from multiprocessing import Pool
import common
import comparators as C
import proposed as Pm
import scenes_v2 as SV
import tune_v2 as TV
import variants

H = common.HERE
RES = os.path.join(H, 'results_v2')
BASE = json.load(open(os.path.join(H, 'tuned_v2', 'Proposed-IMM.json')))['best']['params']
TRACK = ('q_smooth', 'q_maneuver', 'p_stay', 'r_scale')
COMMON = {k: v for k, v in BASE.items() if k not in TRACK}
K_ENDS = (1109., 28083.)


def tuned(m):
    return json.load(open(os.path.join(H, 'tuned_v2', m+'.json')))['best']['params']


class Supervised(common.DirectPolicy):
    """the proposed method's contact-loss supervisor around any incremental controller."""
    def __init__(self, inner, v_app=.006):
        self.inner, self.v_app = inner, v_app
        self.name = inner.name+'+sup'

    def reset(self, platform):
        super().reset(platform)
        self.inner.reset(platform)
        self.lost = False

    def step(self, o):
        row = o.log[o.n-1]
        F, sent = float(row[5]), float(row[11])
        if not self.lost and F < .3:
            self.lost = True
        if self.lost:
            if F >= .5:
                self.lost = False
            else:
                return sent+self.v_app*self.Tc
        return self.inner.step(o)


class FixedFeedback(variants.MotionAwareIMM):
    """identical to the final method except that feedback does not scale with the estimate."""
    def __init__(self, k_fb, **kw):
        super().__init__(**kw)
        self.k_fb = k_fb

    def stiffness(self, o, n):
        super().stiffness(o, n)          # keeps k_hat / k_good updates for excitation and motion tracking
        return self.k_fb, True


def build(name, kf=None):
    if name.endswith('+sup'):
        m = name[:-4]
        return Supervised(TV.make(m, tuned(m)))
    if name == 'Pred-poly':
        pol = Pm.MotionAware(**COMMON)
    elif name == 'Pred-KF':
        pol = variants.build('Proposed-KF', dict(COMMON, **kf))
    elif name.startswith('NoNorm:'):
        pol = FixedFeedback(float(name.split(':')[1]), single_model=False, **BASE)
    else:
        raise ValueError(name)
    pol.name = name
    return pol


def _ep(a):
    name, kf, (fam, sid, p, T) = a
    try:
        log, _ = common.run(build(name, kf), p, T, sid)
        m = common.task_metrics(log, 20.8)
    except Exception as ex:
        return dict(method=name, family=fam, sid=sid, ok=False, err=str(ex), J=1e3)
    ok = not m['diverged']
    rip = m['ripple'] if ok else np.nan
    return dict(method=name, family=fam, sid=sid, ok=bool(ok), rmse=m['rmse'], peak=m['peak'], loss_frac=m['loss'],
                loss_max=m['loss_max'], rmse_noexc=float(math.sqrt(max(m['rmse']**2-rip**2, 0.))) if ok else np.nan,
                J=m['J'])


def tune_kf(pool):
    sc = SV.scenes('dev'); r = np.random.default_rng(777)
    lo, hi = variants.BOUNDS['q_kf'], variants.BOUNDS['r_scale']
    draw = lambda: dict(q_kf=float(math.exp(r.uniform(math.log(lo[0]), math.log(lo[1])))),
                        r_scale=float(math.exp(r.uniform(math.log(hi[0]), math.log(hi[1])))))
    def score(cands):
        res = pool.map(_ep, [('Pred-KF', c, s) for c in cands for s in sc], chunksize=1)
        out = []
        for i, c in enumerate(cands):
            rr = res[i*len(sc):(i+1)*len(sc)]
            out.append((float(np.mean([math.log1p(min(x['J'] if np.isfinite(x['J']) else 1e3, 1e3)) for x in rr])), c))
        return out
    rows = score([dict(variants.DEFAULTS['Proposed-KF'])]+[draw() for _ in range(40)])
    best = min(rows, key=lambda x: x[0])
    for _ in range(4):
        cand = [{k: float(np.clip(v*math.exp(r.normal(0, .35)), *variants.BOUNDS[k])) for k, v in best[1].items()}
                for _ in range(5)]
        rows += score(cand); best = min(rows, key=lambda x: x[0])
    return dict(best=best[1], L=best[0], n_eval=len(rows))


def boot(d, n=4000, seed=0):
    r = np.random.default_rng(seed); d = np.asarray(d)
    b = [d[r.integers(0, len(d), len(d))].mean() for _ in range(n)]
    return float(d.mean()), [float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))]


if __name__ == '__main__':
    with Pool(16) as pool:
        kf = tune_kf(pool)
        print('KF matched tuning', kf, flush=True)
        names = ['RossiType+sup', 'WijType+sup', 'DOMPAC+sup', 'RossiDOMPAC+sup', 'Pred-poly', 'Pred-KF'] + \
                ['NoNorm:%g' % k for k in K_ENDS]
        sc = SV.scenes('confirm')
        rows = pool.map(_ep, [(n, kf['best'], s) for n in names for s in sc], chunksize=1)
    ref = [r for r in json.load(open(os.path.join(RES, 'confirm', 'main.json')))]
    prop = {(r['family'], r['sid']): r for r in ref if r['method'] == 'Proposed'}
    summ = {}
    for fam in SV.FAMILIES:
        for n in names:
            R = [r for r in rows if r['method'] == n and r['family'] == fam]
            ok = [r for r in R if r['ok']]
            d = [r['rmse']-prop[(fam, r['sid'])]['rmse'] for r in ok]
            mu, ci = boot(d) if d else (np.nan, [np.nan, np.nan])
            e = dict(n_fail=len(R)-len(ok), rmse=float(np.mean([r['rmse'] for r in ok])) if ok else np.nan,
                     rmse_noexc=float(np.mean([r['rmse_noexc'] for r in ok])) if ok else np.nan,
                     peak=float(max(r['peak'] for r in ok)) if ok else np.nan,
                     loss=float(np.mean([r['loss_frac'] for r in ok])) if ok else np.nan,
                     loss_max=float(max(r['loss_max'] for r in ok)) if ok else np.nan,
                     minus_final=dict(mean=mu, ci=ci, final_better=int(sum(x > 0 for x in d)), n=len(d)))
            summ[fam+'|'+n] = e
            print('%-9s %-16s fail %d rmse %.3f (noexc %.3f) peak %.1f loss %.4f maxloss %.2f | -final %+.3f [%+.3f,%+.3f] final better %d/%d'
                  % (fam, n, e['n_fail'], e['rmse'], e['rmse_noexc'], e['peak'], e['loss'], e['loss_max'], mu, *ci,
                     e['minus_final']['final_better'], len(d)), flush=True)
    json.dump(dict(kf_matched=kf, summary=summ, rows=rows), open(os.path.join(RES, 'mechanism_checks.json'), 'w'), indent=1)
