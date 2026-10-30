"""Scenario families after the setup audit (2026-10-30).

S2      : static heterogeneous scan (T1 ridges 10-30x, breathing removed) - stiffness span only.
T3m     : spine scan, breathing + 1-2 push events with moderate ramps (0.3-0.5 s).
T2m-soft, T2m-stiff : push/retract on homogeneous soft / stiff tissue, moderate ramps.
T2f-stiff (stress test): stiff tissue with fast ramps (0.15-0.3 s; peak surface speed up to
          ~80 mm/s), reported separately so that a few extreme events do not dominate averages.
Ranges are the benchmark's synthetic engineering values (dose variables), not patient data.
"""
import numpy as np
import common
from common import S

FAMILIES = ('S2', 'T3m', 'T2m-soft', 'T2m-stiff', 'T2f-stiff')
SPLITS = dict(dev=(6, 20000), test=(12, 30000), replication=(12, 40000), select=(10, 80000), confirm=(12, 90000))


def make(family, rng):
    if family == 'S2':
        return S.make('T1', rng, dict(A=0.0))
    if family == 'T3m':
        return S.make('T3', rng, dict(ramp=float(rng.uniform(.3, .5))))
    if family in ('T2m-soft', 'T2m-stiff'):
        return S.make('T2', rng, dict(stiff=family.endswith('stiff'), ramp=float(rng.uniform(.3, .5))))
    if family == 'T2f-stiff':
        return S.make('T2', rng, dict(stiff=True, ramp=float(rng.uniform(.15, .3))))
    raise ValueError(family)


def scenes(split, families=FAMILIES):
    n, base = SPLITS[split]
    out = []
    for fi, fam in enumerate(families):
        rng = np.random.default_rng(base+100*fi)
        for i in range(n):
            p, T = make(fam, rng)
            out.append((fam, base+100*fi+i, p, T))
    return out
