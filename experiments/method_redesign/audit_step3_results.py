"""Audit the delivered result: scope, objective, gates, timing and optional rerun.

python -X utf8 experiments/method_redesign/audit_step3_results.py --rerun --workers 4
The rerun checks one pre-specified scene per family/method, without retuning.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
import json
from pathlib import Path
import numpy as np
from run_step3 import OUT, METHODS, CE, FilePool, episode, sweep
from mpc_robustness import screen


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--rerun', action='store_true')
    ap.add_argument('--workers', type=int, default=4); args = ap.parse_args()
    rec = json.loads((OUT/'step3_mpc.json').read_text(encoding='utf8'))
    rows = rec['rows']; tune = rec['tuning']
    assert len(rows) == 2940 and len(tune['all']) == 61
    assert len({(r['method'], r['family'], r['sid']) for r in rows}) == len(rows)
    assert all(not 130000 <= r['sid'] < 140000 for r in rows)
    for method in METHODS:
        for fam in CE.FAMS:
            assert sum(r['method']==method and r['family']==fam for r in rows) == 50
        for kind, doses in (('ramp', sweep.RAMPS), ('ratio', sweep.RATIOS)):
            for d in doses:
                assert sum(r['method']==method and r['family']==f'{kind}:{d:g}' for r in rows) == 20
    eval_ids = {r['sid'] for r in rows}
    for candidate in tune['all']:
        assert len(candidate['scenes']) == 50
        assert all(190000 <= r['sid'] < 195000 and r['sid'] not in eval_ids for r in candidate['scenes'])
        ll = []
        for fam in CE.FAMS:
            rr = [r for r in candidate['scenes'] if r['family']==fam]
            assert len(rr) == 10
            ll.append(np.mean([np.log1p(min(r['J'], 1e3)) for r in rr]))
        np.testing.assert_allclose(np.mean(ll), candidate['L'], atol=1e-12, rtol=0)
    assert rec['selected'] == min(tune['all'], key=lambda r:r['L'])['params']
    # The delivered implementation retains all dynamic eigenvalues. Recompute
    # every attempted gate, including the rejected candidates.
    max_gate_diff = 0.
    for old in tune['screen_attempts']:
        new = screen(old['params'])
        assert new['accepted'] == old['accepted']
        max_gate_diff = max(max_gate_diff, abs(new['rho_max']-old['rho_max']))
    assert max_gate_diff < 1e-12
    print('PASS counts, no final cohort, no train/evaluation overlap, all objectives and all stability gates')
    times = np.load(OUT/'step3_timings.npz', allow_pickle=False)
    for ss in rec['summary']:
        rr = [r for r in rows if r['family']==ss['family'] and r['method']==ss['method']]
        rm = np.array([r['rmse'] for r in rr if r.get('rmse') is not None])
        np.testing.assert_allclose([np.median(rm), rm.mean()], [ss['rmse_median'], ss['rmse_mean']], atol=1e-12)
        tt = np.concatenate([times[f"{r['method']}__{r['family']}__{r['sid']}__call_ns"] for r in rr])
        np.testing.assert_allclose(np.quantile(tt*1e-6, [.5, .99]), ss['timing_ms'], atol=1e-6, rtol=0)
        if ss['method'].startswith('MPC'):
            qt = np.concatenate([times[f"{r['method']}__{r['family']}__{r['sid']}__solve_ns"] for r in rr])
            np.testing.assert_allclose(np.quantile(qt*1e-6, [.5, .99]), ss['qp_timing_ms'], atol=1e-6, rtol=0)
        base = {r['sid']: r for r in rows if r['method']=='Proposed' and r['family']==ss['family']}
        if ss['paired_bootstrap']:
            delta = [r['rmse']-base[r['sid']]['rmse'] for r in rr if r.get('rmse') is not None and base[r['sid']].get('rmse') is not None]
            np.testing.assert_allclose(np.mean(delta), ss['paired_bootstrap']['mean_delta'], atol=1e-12)
    print('PASS independent summary and individual-duration artifact checks')
    result = dict(scope='2940 evaluation + 3050 tuning episodes; 83 attempted gates',
                  all_counts_objectives_gates_metrics_timings_passed=True,
                  max_gate_recompute_error=max_gate_diff)
    if args.rerun:
        sc = CE.cohort('diag', 50)
        chosen = [sc[50*fi] for fi in range(5)]
        tasks = [(m, rec['selected'], s) for m in METHODS for s in chosen]
        with FilePool(args.workers) as pool: repeated = pool.map(episode, tasks)
        differences = []
        keys = ('J', 'rmse', 'peak', 'loss', 'loss_max', 'diverged')
        for new in repeated:
            old = next(r for r in rows if r['method']==new['method'] and r['sid']==new['sid'])
            assert 'error' not in new
            diff = {key: abs(float(new[key])-float(old[key])) for key in keys}
            differences.append(dict(method=new['method'], family=new['family'], sid=new['sid'], absolute_difference=diff))
        result['rerun_30_scenes'] = differences
        result['rerun_max_absolute_difference'] = {key:max(d['absolute_difference'][key] for d in differences) for key in keys}
        # Numerical solver tolerances do not imply bitwise trajectory equality.
        result['rerun_close'] = all(d['absolute_difference']['rmse'] < 1e-4 and d['absolute_difference']['peak'] < 1e-3 for d in differences)
        print('RERUN', json.dumps(result['rerun_max_absolute_difference']), 'close=', result['rerun_close'])
    rec['audit'] = result
    (OUT/'step3_mpc.json').write_text(json.dumps(rec, ensure_ascii=False, indent=1), encoding='utf8')
    print('AUDIT SAVED')


if __name__ == '__main__': main()
