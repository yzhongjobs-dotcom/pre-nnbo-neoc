"""Where does lead-filter selection lose on fast pushes? Phase-resolved RMSE (true force) and selected candidate per
phase, using the push events in the scene vector (p[16+4k]: onset, amplitude, ramp, hold). Diagnosis cohort,
fast-push family (20 scenes) and stiff push. Run: python -X utf8 p2_sel_phases.py
"""
import os, sys
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import p2_hf_attribution as H
from p2_hf_attribution import common, CE, P
CFG = {'Proposed': dict(), 'SelServo0.3': dict(sched='select', servo_crit=True)}
PH = ('quiet', 'ramp-in', 'hold', 'ramp-out')
def phases(t, p):
    lab = np.zeros(len(t), int)
    for k in range(int(p[15])):
        on, amp, ramp, hold = p[16+4*k:20+4*k]
        lab[(t >= on) & (t < on+ramp)] = 1; lab[(t >= on+ramp) & (t < on+ramp+hold)] = 2
        lab[(t >= on+ramp+hold) & (t < on+2*ramp+hold)] = 3
    return lab
class Rec(H.Probe):
    def drift_pred(self, o, n):
        r = super().drift_pred(o, n)
        st = getattr(self, '_sel', None); self.sel_t = getattr(self, 'sel_t', [])
        self.sel_t.append((o.log[n-1, 0], st['i'] if st and st.get('i') is not None else 0))
        return r
def _ep(a):
    cfg, (fam, sid, p, T) = a
    pol = H.mk(**CFG[cfg]); pol.__class__ = Rec; pol.p_scene = p
    log, _ = common.run(pol, p, T, sid)
    t = log[:, 0]; e = log[:, 3]-P.FD; sel = t >= 1.; lab = phases(t, p)
    out = [float(np.sqrt(np.mean(e[sel & (lab == i)]**2))) if np.any(sel & (lab == i)) else np.nan for i in range(4)]
    st = np.array(pol.sel_t) if getattr(pol, 'sel_t', None) else None
    fr = []
    if st is not None and cfg != 'Proposed':
        ls = phases(st[:, 0], p)
        fr = [float(np.mean(st[ls == i, 1] == 0)) if np.any(ls == i) else np.nan for i in range(4)]
    return cfg, fam, out, fr
if __name__ == '__main__':
    sc = [s for s in CE.cohort('diag', 20) if s[0] in ('T2f-stiff', 'T2m-stiff')]
    with Pool(16) as po: rows = po.map(_ep, [(c, s) for c in CFG for s in sc], chunksize=1)
    for fam in ('T2m-stiff', 'T2f-stiff'):
        for c in CFG:
            v = np.array([r[2] for r in rows if r[0] == c and r[1] == fam])
            line = '%-9s %-12s ' % (fam, c)+' '.join('%s %.3f/%.3f' % (PH[i], np.nanmedian(v[:, i]), np.nanmean(v[:, i])) for i in range(4))
            fr = [r[3] for r in rows if r[0] == c and r[1] == fam and r[3]]
            if fr: line += ' | raw share ' + ' '.join('%.2f' % x for x in np.nanmean(fr, 0))
            print(line, flush=True)
