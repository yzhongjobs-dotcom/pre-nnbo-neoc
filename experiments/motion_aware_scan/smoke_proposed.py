import time, numpy as np, common
from proposed import MotionAware, Oracle, FixedEnd
from common import S
for task, sd in (('T1', 13), ('T2', 11), ('T3', 12)):
    p, T = S.make(task, np.random.default_rng(sd))
    for pol in (MotionAware(), MotionAware(feedforward=False, tag='NoFF'), Oracle(p), FixedEnd(5e2), FixedEnd(3e4)):
        t0 = time.time(); log, _ = common.run(pol, p, T, 0)
        m = common.task_metrics(log, 20.8)
        kt = log[-1, 8] if log is not None else np.nan
        print(f"{task} {pol.name:12s} rmse={m['rmse']:.3f} peak={m['peak']:.2f} nloss={m['n_loss']} osc={m['osc_ex']:.3f} "
              f"ripple={m['ripple']:.3f} kbar_end={pol.kbar:.0f} k_true_end={kt:.0f} wall={time.time()-t0:.1f}s")
