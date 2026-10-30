"""Run eight deterministic smoke episodes; no tuning or result files.

python -X utf8 experiments/motion_aware_scan/smoke_comparators.py
Timing covers policy calls, including estimation, excluding reset/plant/JIT.
Finite trajectories do not establish satisfactory force tracking or stability.
"""
import sys
sys.dont_write_bytecode = True
import numpy as np
try:
    from . import common
    from .comparators import RossiType, WijType, DOMPAC, RossiDOMPAC
except ImportError:
    import common
    from comparators import RossiType, WijType, DOMPAC, RossiDOMPAC


def main():
    print('method         task rmse[N] peak[N] n_loss osc_ex[N] ripple[N] '
          'k_hat[N/m] k_eff[N/m] mean_us finite')
    failed = False
    for cls in (RossiType, WijType, DOMPAC, RossiDOMPAC):
        for task, scenario_seed in (('T2', 11), ('T3', 12)):
            p, T = common.S.make(task, np.random.default_rng(scenario_seed))
            policy = cls()
            log, policy = common.run(policy, p, T, seed=0)
            metrics = common.task_metrics(log, getattr(policy, 'f_exc', None))
            finite = (log is not None and len(log) == round(T/policy.Tc)
                      and np.isfinite(log).all()
                      and np.isfinite(np.asarray(policy.trace)).all()
                      and len(policy.trace) == len(log)-1)
            failed |= not finite
            k_last = policy.trace[-1][1] if policy.trace else np.nan
            k_true = log[-1, 8] if log is not None else np.nan
            mean_us = 1e6*np.mean(policy.t_call[1:])
            print(f'{policy.name:14s} {task} {metrics["rmse"]:7.4f} '
                  f'{metrics["peak"]:7.3f} {metrics["n_loss"]:6} '
                  f'{metrics["osc_ex"]:9.4f} {metrics["ripple"]:9.4f} '
                  f'{k_last:10.1f} {k_true:10.1f} {mean_us:7.2f} {finite}')
    if failed:
        raise RuntimeError('Incomplete or nonfinite episode/trace')


if __name__ == '__main__':
    main()
