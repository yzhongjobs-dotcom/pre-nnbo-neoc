"""Stage 1: continuous parameter paths in the existing Appendix M.3 root.

No residual/rollout acceptance check in predictor execution. Exact roots are
used ONLY for the initial/scheduled anchors and retrospective evaluation.
This is not a physical R15 validation or a continuous-HJB optimality proof.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import platform

OUT = Path(__file__).resolve().parent
LOW, HIGH = .65, 1.05
TOL = .001  # W error, normalized; diagnostic choice, not robot requirement
M = 2 / (1 + 4 * LOW**2)**1.5
A = 4 * HIGH  # |d[2 eta/(1+2 W)]/dW| <= 4 HIGH for W >= 0


def root(eta):
    return (math.sqrt(1 + 4 * eta * eta) - 1) / 2


def sensitivity(w, eta):
    return 2 * eta / (1 + 2 * w)


def path(kind, n):
    if kind == 'ramp':
        return [LOW + (HIGH - LOW) * i / n for i in range(n + 1)]
    # Four cycles: total variation is 3.2 despite returning to the anchor.
    return [.85 + .2 * math.sin(8 * math.pi * i / n) for i in range(n + 1)]


def execute(method, seq):
    w = root(seq[0]); anchor_eta = seq[0]; anchor_w = w
    anchor_c = sensitivity(w, anchor_eta)
    bound = 0.; refresh = 0; sensitivity_evals = 1
    residual_evals = 0; nonlinear_solves = 1
    trace = [(0, seq[0], w, bound, False)]
    for i in range(1, len(seq)):
        eta = seq[i]; prev = seq[i-1]; delta = eta-prev; d = abs(delta)
        reset = False
        if method in ('frozen_sensitivity', 'budgeted_frozen_NEOC'):
            if method == 'budgeted_frozen_NEOC' and M/2*(eta-anchor_eta)**2 > TOL:
                anchor_eta=eta; anchor_w=root(eta)
                anchor_c=sensitivity(anchor_w,eta)
                nonlinear_solves += 1; sensitivity_evals += 1
                refresh += 1; reset=True
            w = anchor_w + anchor_c * (eta-anchor_eta)
            bound = M / 2 * (eta-anchor_eta)**2
        elif method in ('continuous_NEOC', 'budgeted_anchor_NEOC'):
            proposed = (1 + A*d)*bound + M/2*d*d
            if method == 'budgeted_anchor_NEOC' and proposed > TOL:
                # Predetermined error budget, no candidate residual or true root
                # error queried to decide this. Reset has a real solve cost.
                w = root(eta); bound = 0.; refresh += 1
                nonlinear_solves += 1; reset = True
            else:
                w += sensitivity(w, prev)*delta
                sensitivity_evals += 1; bound = proposed
        elif method == 'one_Newton':
            residual = eta*eta-w-w*w
            w += residual/(1+2*w); residual_evals += 1
        elif method == 'fixed_policy':
            pass
        else:
            raise ValueError(method)
        if not (0 <= w < .7):
            raise ArithmeticError('Inactive-branch domain violated')
        trace.append((i, eta, w, bound, reset))
    counts = dict(anchor_solves=nonlinear_solves, additional_anchors=refresh,
                  sensitivity_evaluations=sensitivity_evals,
                  residual_evaluations=residual_evals,
                  candidate_acceptance_checks=0)
    return trace, counts


def main():
    summaries = []; all_rows = []
    methods = ['frozen_sensitivity', 'budgeted_frozen_NEOC', 'continuous_NEOC',
               'budgeted_anchor_NEOC', 'one_Newton', 'fixed_policy']
    for kind in ['ramp', 'four_cycles']:
        for n in [100, 200, 400, 800]:
            seq = path(kind, n)
            exact = [root(eta) for eta in seq]  # evaluator only
            for method in methods:
                trace, counts = execute(method, seq)
                errs = [abs(row[2]-exact[i]) for i, row in enumerate(trace)]
                certified = method in methods[:4]
                violations = sum(err > row[3] + 1e-12
                                 for err, row in zip(errs, trace)) if certified else None
                if certified:
                    assert violations == 0
                if method in ('budgeted_anchor_NEOC','budgeted_frozen_NEOC'):
                    assert max(row[3] for row in trace) <= TOL + 1e-12
                summaries.append(dict(path=kind, updates=n, method=method,
                    max_step=max(abs(seq[i]-seq[i-1]) for i in range(1, len(seq))),
                    total_variation=sum(abs(seq[i]-seq[i-1]) for i in range(1,len(seq))),
                    max_W_error=max(errs), final_W_error=errs[-1],
                    max_error_bound=max(row[3] for row in trace) if certified else None,
                    bound_violations=violations, **counts))
                for row, err in zip(trace, errs):
                    i, eta, w, bound, reset = row
                    all_rows.append(dict(path=kind, updates=n, method=method,
                        step=i, eta=eta, W=w, reference_W=exact[i], error=err,
                        error_bound=bound if certified else '', anchor_reset=reset))
    # Scalar formula and implicit sensitivity agree; finite difference is audit.
    for eta in [.65,.8,1.,1.05]:
        numerical = (root(eta+1e-6)-root(eta-1e-6))/2e-6
        assert abs(numerical-sensitivity(root(eta),eta)) < 1e-9
    for kind in ['ramp','four_cycles']:
        items = [r for r in summaries if r['path']==kind and r['method']=='continuous_NEOC']
        assert all(items[i+1]['max_W_error'] < items[i]['max_W_error']
                   for i in range(len(items)-1))
    for name, rows in [('summary.csv',summaries), ('trace.csv',all_rows)]:
        with (OUT/name).open('w',newline='',encoding='utf-8') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    result=dict(stage='1_CONTINUOUS_ROOT_MECHANISM_ONLY',
                model='Appendix M.3 inactive finite root W+W^2=eta^2',
                eta_interval=[LOW,HIGH], W_tolerance=TOL,
                M_second_derivative_bound=M, A_sensitivity_error_bound=A,
                online_candidate_acceptance_checks=0,
                physical_closed_loop_tested=False, wall_time_advantage_tested=False,
                python=platform.python_version(),
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                summaries=summaries)
    (OUT/'results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps([r for r in summaries if r['updates'] in [100,800]],indent=2))


if __name__ == '__main__':
    main()
