"""What is happening just before each runaway of the isolated filter? Truth is scoring only."""
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from diagnose import Isolated, TUNED, NOISE, FINAL, common, cohort, MeasuredOnly, RES, ISO
from v1_servo_id import ServoID


class Traced(Isolated):
    def _stiffness_est(self, o, n):
        r = super()._stiffness_est(o, n)
        self.tr.append((float(o.log[n-1, 0]), float(self.kbar), float(self.k_good), float(self.k_hat),
                        float(self.rel) if np.isfinite(self.rel) else np.nan))
        return r

    def reset(self, p):
        super().reset(p); self.tr = []


def analyse(job):
    scene = job; fam, sid, p, T = scene
    q = dict(FINAL)
    q.update(q_lnk0=TUNED['q_lnk0'], q_lnk1=TUNED['q_lnk1'], p_stay_stiff=TUNED['p_stay_stiff'],
             consistency_max=TUNED['consistency_max'], c=TUNED['c'], rho_d=TUNED['rho_d'], gain_rel_max=.1823)
    pol = Traced(NOISE*TUNED['noise_scale'], TUNED['r_scale'], TUNED['p_stay'], single_model=False, **q)
    log, _ = common.run(MeasuredOnly(pol), p, T, sid, platform=dict(common.P.PLATFORMS['PC125']))
    if log is None:
        log_n = 0
        return dict(family=fam, sid=sid, diverged=True)
    tr = np.asarray(pol.tr)
    f = log[:, 3]; t = log[:, 0]; kt = log[:, 8]
    bad = np.flatnonzero((f > 12.) & (t > 1.))
    if not len(bad):
        return dict(family=fam, sid=sid, peak=float(f.max()), onset=None)
    i0 = bad[0]; t0 = t[i0]
    m = (tr[:, 0] > t0-.3) & (tr[:, 0] <= t0)
    idx = np.searchsorted(t, tr[m, 0]); idx = np.clip(idx, 0, len(t)-1)
    kg = kt[idx]/np.maximum(tr[m, 2], 1.); kap = kt[idx]/np.maximum(tr[m, 1], 1.)
    khat = tr[m, 3]/np.maximum(kt[idx], 1.)
    ht = np.array([abs(common.P.surface(float(log[j, 1]), float(log[j, 0]), np.asarray(p, float))[1]) for j in idx])
    return dict(family=fam, sid=sid, peak=float(f.max()), onset=float(t0),
                kg_med=float(np.median(kg)), kg_max=float(np.max(kg)), kappa_med=float(np.median(kap)),
                kappa_max=float(np.max(kap)), khat_over_true_med=float(np.median(khat)), rel_med=float(np.nanmedian(tr[m, 4])),
                surface_speed_max=float(ht.max()), k_true=float(kt[i0]), k_true_ratio_300ms=float(kt[i0]/max(kt[max(i0-37, 0)], 1.)))


def main():
    rows = json.loads((RES/'diagnose_isolated.json').read_text())['rows']
    bad = {(r['family'], r['sid']) for r in rows if r['config'] == 'iso_band_gate' and (r['rmse'] is None or r['peak'] > 50)}
    scenes = [s for s in cohort('diag', 50) if (s[0], s[1]) in bad]
    with ProcessPoolExecutor(max_workers=15) as pool:
        out = list(pool.map(analyse, scenes))
    for r in out:
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
    (RES/'diagnose_failures.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
