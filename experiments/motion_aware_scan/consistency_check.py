"""Theory/fairness consistency checks for the final method (paper draft §三.3, §五.1).

1) rho_d capped so that the commutation allowance used when screening A1 (delta = 0.07,
   delta ~ (tau_d+tau_s) rho_d / Tc) covers the final parameters: rho_d = 0.07*Tc/(Tc+tau_s).
2) recovery supervisor removed (the comparators have no such layer), and 3) reduced to the gate that the
   feed-forward loop analysis requires (FF/excitation off and k_bar frozen while lost; no approach mode).
All other parameters are the selected Proposed-IMM values; confirmation cohort; paired against
the stored confirmation rows of the comparators (results_v2/confirm/main.json).
Run: python -X utf8 consistency_check.py  -> results_v2/consistency.json
"""
import json, math, os
import numpy as np
from multiprocessing import Pool
import common
import scenes_v2 as SV
import variants

H = common.HERE
BASE = json.load(open(os.path.join(H, 'tuned_v2', 'Proposed-IMM.json')))['best']['params']
RHO_CAP = 0.07*0.008/(0.008+0.019)
VARIANTS = {'final': {}, 'rho_capped': dict(rho_d=min(BASE['rho_d'], RHO_CAP)), 'no_recovery': dict(recovery=False),
            'gate_only': dict(recovery='gate')}


def _ep(a):
    v, (fam, sid, p, T) = a
    pol = variants.build('Proposed-IMM', dict(BASE, **VARIANTS[v]))
    log, _ = common.run(pol, p, T, sid)
    m = common.task_metrics(log, 20.8)
    return dict(variant=v, family=fam, sid=sid, ok=not m['diverged'], rmse=m['rmse'], peak=m['peak'],
                loss_frac=m['loss'], J=m['J'])


def boot(d, n=4000, seed=0):
    r = np.random.default_rng(seed); d = np.asarray(d)
    b = [d[r.integers(0, len(d), len(d))].mean() for _ in range(n)]
    return float(d.mean()), [float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))]


if __name__ == '__main__':
    sc = SV.scenes('confirm')
    with Pool(16) as po:
        rows = po.map(_ep, [(v, s) for v in VARIANTS for s in sc])
    ref = json.load(open(os.path.join(H, 'results_v2', 'confirm', 'main.json')))
    out = dict(rho_cap=RHO_CAP, rows=rows, summary={})
    for fam in SV.FAMILIES:
        for v in VARIANTS:
            R = {r['sid']: r for r in rows if r['variant'] == v and r['family'] == fam}
            ok = [r for r in R.values() if r['ok']]
            e = dict(n_fail=len(R)-len(ok), rmse=float(np.mean([r['rmse'] for r in ok])),
                     peak=float(max(r['peak'] for r in ok)), loss=float(np.mean([r['loss_frac'] for r in ok])))
            for m in ('RossiType', 'WijType', 'DOMPAC', 'Proposed'):
                C = {r['sid']: r for r in ref if r['method'] == m and r['family'] == fam and r['ok']}
                d = [C[s]['rmse']-R[s]['rmse'] for s in C if s in R and R[s]['ok']]
                mu, ci = boot(d)
                e['vs_'+m] = dict(mean=mu, ci=ci, wins=int(sum(x > 0 for x in d)), n=len(d))
            out['summary'][fam+'|'+v] = e
            print(fam, v, 'fail', e['n_fail'], 'rmse %.3f peak %.1f loss %.4f' % (e['rmse'], e['peak'], e['loss']),
                  ' '.join('%s %+.3f[%+.3f,%+.3f] %d/%d' % (m[:5], e['vs_'+m]['mean'], *e['vs_'+m]['ci'],
                                                          e['vs_'+m]['wins'], e['vs_'+m]['n'])
                           for m in ('RossiType', 'WijType', 'DOMPAC', 'Proposed')))
    json.dump(out, open(os.path.join(H, 'results_v2', 'consistency.json'), 'w'), indent=1)
