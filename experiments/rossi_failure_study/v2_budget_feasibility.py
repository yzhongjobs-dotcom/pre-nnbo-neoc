"""Check V2's deterministic Huber-budget interpretation on the original plant.

This is a diagnostic, not an online estimator or an experiment demonstrating
unidentifiability of the complete physical model. Truth is used offline only.
It tests the deliberately broad residual class D={d: sum Huber(D4 d/s)<=Gamma}.
This class includes phase/model error and window variation when calibrated from
z0-k_end*Y0; it must NOT be described as isolated anatomical motion.

Default: reuse/generate four existing development fixtures (A, B, cough,
A_large). Use the first half to choose each group's residual budget; test two
alternative constant gain explanations on the second half, without sensor
noise. Each alternative has its own residual but exactly the same observation.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[_key] = '1'
import argparse
import json
from pathlib import Path
import numpy as np
import huber_windows as W
import huber_prior_check as H

HERE = Path(__file__).resolve().parent
CASES = (0, 12, 24, 36)


def penalty(x, c=1.345):
    x = np.asarray(x); a = np.abs(x)
    return float(np.sum(np.where(a <= c, x*x/2, c*a-c*c/2)))


def check_case(index):
    data = H.load(H.WIN/f'case_{index:02d}.npz')
    log = data['log']; t = log[:, H.COL[0]]; k = log[:, H.COL[26]]
    D = H.diff_mat(24); scale = H.SIGMA/np.sqrt(H.RHO)
    records = []
    for j in range(96, len(t), 6):
        y = H.tail(data['Y0'][j], 24); z = H.tail(data['z0'][j], 24)
        if len(y) != 24 or np.min(k[j-24:j+1]) <= 200:
            continue
        center = float(k[j])
        variants = [(fraction, z-fraction*center*y) for fraction in (.75, 1., 1.25)]
        budgets = [penalty(D@d/scale) for _, d in variants]
        # The alternative residuals reconstruct the identical noiseless z.
        equality_error = max(float(np.max(np.abs(f*center*y+d-z))) for f,d in variants)
        event = any(e['start'] <= t[j] <= e['end']+.2 for e in data['meta']['event_info'])
        records.append(dict(t=float(t[j]), k_offline=center, budget_center=budgets[1],
                            budget_low=budgets[0], budget_high=budgets[2],
                            pair_required=max(budgets[0],budgets[2]),
                            alias_only_budget=penalty(.25*center*(D@y)/scale),
                            equality_error=equality_error, event=event,
                            window_log_change=float(np.log(k[j]/k[j-24]))))
    split = float(t[-1]/2)
    dev = [r for r in records if r['t'] < split]
    test = [r for r in records if r['t'] >= split]
    gamma99 = float(np.percentile([r['budget_center'] for r in dev],99))
    summary = dict(index=index,group=data['meta']['group'],n_dev=len(dev),n_test=len(test),
                   Gamma_dev99=gamma99,
                   test_center_coverage=float(np.mean([r['budget_center'] <= gamma99 for r in test])),
                   test_both_alternatives_feasible_share=float(np.mean([r['pair_required'] <= gamma99 for r in test])),
                   equality_error_max=max(r['equality_error'] for r in records),
                   test_alias_budget_median=float(np.median([r['alias_only_budget'] for r in test])),
                   dev_budget_median=float(np.median([r['budget_center'] for r in dev])))
    return dict(summary=summary,records=records,meta=data['meta'])


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--generate',action='store_true')
    args = ap.parse_args()
    W.OUT.mkdir(exist_ok=True)
    for i in CASES:
        if not (W.OUT/f'case_{i:02d}.npz').exists():
            if not args.generate:
                raise FileNotFoundError(f'case {i}: rerun with --generate')
            W.job(W.CC.cases()[i])
    result = dict(scope='four existing development cases; diagnostic of a broad deterministic residual class, not an online algorithm',
                  repo_commit='1afb6fc920e0c8e298d1ba5c2d749540dfe30890',
                  observation='Y0,z0: recorded commands replayed without sensor noise, original calibrated sensing',
                  sigma_z=H.SIGMA, rho=H.RHO, window=24, hop=6,
                  alternatives='constant gains 0.75*k_end and 1.25*k_end; k_end used offline for diagnosis only',
                  cases=[check_case(i) for i in CASES])
    out = HERE/'results/v2_budget_feasibility.json'
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps([r['summary'] for r in result['cases']],ensure_ascii=False,indent=2))


if __name__ == '__main__':
    main()
