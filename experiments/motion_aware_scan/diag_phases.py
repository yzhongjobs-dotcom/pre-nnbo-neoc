"""Development-set diagnosis: where does the remaining force error come from?
Splits squared force error into push-ramp, push-hold, breathing-only and stiffness-transition
phases (truth used for labelling only). Dev seeds only; test cohorts untouched."""
import numpy as np, json, common, tune_equal as tu
import run_experiments as RX
from common import P

def phases(log, p):
    t = log[:, 0]; lab = np.full(len(t), 'breath', object)
    for k in range(int(p[15])):
        t0, amp, tr, hold = p[16+4*k:20+4*k]
        lab[(t >= t0) & (t < t0+tr)] = 'ramp'
        lab[(t >= t0+tr) & (t < t0+tr+hold)] = 'hold'
        lab[(t >= t0+tr+hold) & (t < t0+2*tr+hold)] = 'ramp'
    kt = log[:, 8]; dk = np.abs(np.gradient(np.log(np.maximum(kt, 1.))))/0.008
    lab[(dk > 2.0) & (lab == 'breath')] = 'stiff-edge'
    return lab

if __name__ == '__main__':
    sc = tu.scenes([('T2', 6, 9200), ('T3', 6, 9300)])
    for name in ('Proposed', 'RossiType', 'WijType', 'Oracle'):
        acc = {}
        for task, sid, p, T in sc:
            tag = task if task != 'T2' else ('T2s' if p[0] > 2e5 else 'T2f')
            log, _ = common.run(RX.build(name, p), p, T, sid)
            if log is None: continue
            sel = log[:, 0] >= 1.0; e2 = (log[sel, 3]-P.FD)**2; lab = phases(log, p)[sel]
            for ph in ('breath', 'ramp', 'hold', 'stiff-edge'):
                m = lab == ph
                a = acc.setdefault((tag, ph), [0., 0]); a[0] += e2[m].sum(); a[1] += m.sum()
        out = {}
        for (tag, ph), (s2, n) in sorted(acc.items()):
            if n: out.setdefault(tag, []).append(f'{ph}: rms={np.sqrt(s2/n):.2f} N, share={100*s2/sum(v[0] for (tg,_),v in acc.items() if tg==tag):.0f}%, time={n*0.008:.1f}s')
        print(name); [print('  ', k, '|', '; '.join(v)) for k, v in out.items()]
