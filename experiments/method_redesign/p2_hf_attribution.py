"""Which path injects the 23-62 Hz force error of the final method? Diagnosis cohort, 20 per family.
NoNoise: plant force-sensor noise scale 0 (diagnosis only). KPF0.5: certified A1 with the same D1 and K_PF1=0.5.
FF-noVel: motion feed-forward uses the predicted position s only (drops tau_s * s_dot). Also logs the IMM
maneuver-mode probability (median over normal-contact ticks). Run: python -X utf8 p2_hf_attribution.py
"""
import os, sys
import numpy as np
from scipy.signal import welch
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'motion_aware_scan')); sys.path.insert(0, os.path.join(HERE, '..', 'cobot_scan_benchmark'))
import common, variants, cohort_eval as CE
from common import P
F = dict(CE.FINAL)
class Probe(variants.MotionAwareIMM):
    novel = False; notch = 0.; vel_lp = 0.; sched = None
    def reset(self, platform):
        super().reset(platform); self.mu_log = []
        if self.notch:
            # known-input notch at the excitation frequency on the motion measurement before the tracker
            w0 = 2*np.pi*self.params['f_exc']*self.Tc; r = self.notch
            b = np.array([1., -2*np.cos(w0), 1.]); a = np.array([1., -2*r*np.cos(w0), r*r]); b *= a.sum()/b.sum()
            tr = self.tracker; up, rs = tr.update, tr.reset; st = {'x': None}
            def reset_(z0):
                st['x'] = [z0, z0]; st['y'] = [z0, z0]; return rs(z0)
            def update_(z):
                if st['x'] is None: st['x'] = [z, z]; st['y'] = [z, z]
                yv = b[0]*z+b[1]*st['x'][0]+b[2]*st['x'][1]-a[1]*st['y'][0]-a[2]*st['y'][1]
                st['x'] = [z, st['x'][0]]; st['y'] = [yv, st['y'][0]]; return up(yv)
            tr.reset, tr.update = reset_, update_
    def drift_pred(self, o, n):
        r = super().drift_pred(o, n)
        if self.tracker.initialized: self.mu_log.append(float(self.tracker._mu[-1]))
        if self.sched == 'wiener' and isinstance(r, tuple):
            # two-band lead: low band (v_LP, 3 Hz) always; high band h = v - v_LP passed with the Wiener gain
            # g = max(0, 1 - N/P_h); P_h short-time power of h, N its noise floor by minimum statistics
            L = (self.ndel+.5)*self.Tc; s, sd = self.tracker.predict(L)
            st = getattr(self, '_wn', None)
            if st is None or self._fin_prev is None: st = self._wn = dict(vl=sd, P=0., N=np.inf)
            st['vl'] += (1-np.exp(-2*np.pi*self.w_fc*self.Tc))*(sd-st['vl']); h = sd-st['vl']
            st['P'] += (1-np.exp(-self.Tc/self.w_T))*(h*h-st['P'])
            st['k'] = st.get('k', 0)+1
            if st['k'] < 10: st['N'] = np.inf                                      # warm-up of P (10 ticks)
            else: st['N'] = min(st['P'], st['N']*np.exp(self.Tc/self.w_floor*np.log(2.)))   # floor may double per w_floor s
            g = (max(0., 1.-self.w_beta*st['N']/st['P']) if st['P'] > 0 else 0.) if np.isfinite(st['N']) else 1.
            self.g_log = getattr(self, 'g_log', [])+[g]
            r = (None, s+self.tau*(st['vl']+g*h)); self._motion_prediction = r
            return r
        if self.sched == 'shrink' and isinstance(r, tuple):
            # posterior shrinkage of the servo-lag lead: g = v^2/(v^2 + c P_vv) (tracker's own mixture covariance)
            L = (self.ndel+.5)*self.Tc; s, sd = self.tracker.predict(L); Pvv = float(self.tracker.P[1, 1])
            g = sd*sd/(sd*sd+self.shrink_c*Pvv+1e-30); self.g_log = getattr(self, 'g_log', [])+[g]
            r = (None, s+self.tau*g*sd); self._motion_prediction = r
            return r
        if self.sched == 'select' and isinstance(r, tuple):
            # causal selection of the lead filter by its own task: predict the motion measurement at the actuation
            # time (L + tau ahead); score = exponentially weighted squared error against the z measured later
            L = (self.ndel+.5)*self.Tc; s, sd = self.tracker.predict(L); h = max(int(round((L+self.tau)/self.Tc)), 1)
            st = getattr(self, '_sel', None)
            if st is None or self._fin_prev is None or not self.tracker.initialized:
                st = self._sel = dict(vf=np.full(len(self.cands), sd), J=np.zeros(len(self.cands)), buf=[]); 
            al = np.array([1. if f == np.inf else 1-np.exp(-2*np.pi*f*self.Tc) for f in self.cands])
            st['vf'] = st['vf']+al*(sd-st['vf'])
            pred = s+self.tau*st['vf']
            z = float(self.tracker.x[0]) if False else self.Z
            lam = np.exp(-self.Tc/self.sel_T)
            if self.servo_crit:
                # exact criterion: candidate feed-forward through the identified servo (delay ndel+1, first-order lag
                # tau) versus the measured motion measurement now
                st['buf'].append(pred.copy()); D = self.ndel+1
                if len(st['buf']) > D:
                    pd = st['buf'].pop(0); asv = 1-np.exp(-self.Tc/max(self.tau, 1e-4))
                    st['yc'] = pd.copy() if st.get('yc') is None else st['yc']+asv*(pd-st['yc'])
                    if self.force_crit:
                        # counterfactual force error of candidate c: measured error minus k (servo FF_c - servo FF_sel)
                        isel = st.get('i', 0) if st.get('i') is not None else 0
                        ku = self.k_good if np.isfinite(self.k_good) else self.kbar
                        ec = (o.log[n-1, 5]-4.0)+ku*(st['yc']-st['yc'][isel])
                        st['J'] = lam*st['J']+(1-lam)*ec**2
                    else:
                        st['J'] = lam*st['J']+(1-lam)*(st['yc']-z)**2
                    if self.fast_T:
                        lf = np.exp(-self.Tc/self.fast_T); st['Jf'] = lf*st.get('Jf', np.zeros_like(st['J']))+(1-lf)*(st['yc']-z)**2
            else:
                st['buf'].append((n+h, pred.copy()))
                while st['buf'] and st['buf'][0][0] <= n:
                    _, pp = st['buf'].pop(0); st['J'] = lam*st['J']+(1-lam)*(pp-z)**2
            i = int(np.argmin(np.maximum(st['J'], st['Jf']) if self.fast_T and 'Jf' in st else st['J'])); self.sel_log = getattr(self, 'sel_log', np.zeros(len(self.cands)))
            self.sel_log[i] += 1; off = st.get('off', 0.)
            if self.bumpless and st.get('i') is not None and st['i'] != i: off += float(pred[st['i']]-pred[i])
            st['i'] = i; st['off'] = off*np.exp(-self.Tc/.05) if self.bumpless else 0.
            r = (None, float(pred[i])+off); self._motion_prediction = r
            return r
        if self.vel_lp and isinstance(r, tuple):
            # first-order low-pass on the servo-lag lead only (filtered differentiator); position prediction unchanged
            s, sd = self.tracker.predict((self.ndel+.5)*self.Tc); a = 1-np.exp(-2*np.pi*self.vel_lp*self.Tc)
            if self.sched == 'oracle':
                kg = max(float(o.log[n-1, 8]), 1.)/self.k_good if np.isfinite(self.k_good) else np.inf
                a = a if abs(np.log(kg)) > np.log(1.15) else 1.
            elif self.sched == 'spatial':
                # diagnosis only: footprint currently spans a modulus boundary (true field, true depth)
                import plant as PL
                xc, d = float(o.log[n-1, 1]), max(float(o.log[n-1, 6]), 1e-6); hw = PL.halfwidth(d)
                Ev = [PL.modulus(xc+u*hw, self.p_scene) for u in np.linspace(-1.2, 1.2, 13)]
                a = a if max(Ev)/min(Ev) > 1.15 else 1.
            elif self.sched == 'est':
                bad = self.fit is not None and np.isfinite(self.fit.get('chi2', np.nan)) and self.fit['chi2'] > self.params['chi2_max']
                if np.isfinite(self.k_check) and np.isfinite(self.k_good): bad = bad or abs(np.log(self.k_check/self.k_good)) > np.log(1.15)
                a = a if bad else 1.
            self.sched_log = getattr(self, 'sched_log', 0)+(a < 1)
            vf = getattr(self, '_vf', None); vf = sd if (vf is None or self._fin_prev is None) else vf+a*(sd-vf)
            self._vf = vf; r = (None, s+self.tau*vf); self._motion_prediction = r
        if self.novel and isinstance(r, tuple):
            s, _ = self.tracker.predict((self.ndel+.5)*self.Tc); r = (None, s); self._motion_prediction = r
        return r
def mk(**kw):
    pol = variants.build('Proposed-IMM', dict(F, **{k: v for k, v in kw.items() if k not in ('novel', 'notch', 'vel_lp', 'sched', 'cands', 'sel_T', 'bumpless', 'servo_crit', 'fast_T', 'force_crit', 'shrink_c', 'w_fc', 'w_T', 'w_floor', 'w_beta')})); pol.__class__ = Probe
    pol.novel = kw.get('novel', False); pol.notch = kw.get('notch', 0.); pol.vel_lp = kw.get('vel_lp', 0.); pol.sched = kw.get('sched'); pol.cands = kw.get('cands', (np.inf, 10., 6., 3.)); pol.sel_T = kw.get('sel_T', .3); pol.bumpless = kw.get('bumpless', False); pol.servo_crit = kw.get('servo_crit', False); pol.fast_T = kw.get('fast_T', 0.); pol.force_crit = kw.get('force_crit', False); pol.shrink_c = kw.get('shrink_c', 1.); pol.w_fc = kw.get('w_fc', 3.); pol.w_T = kw.get('w_T', .03); pol.w_floor = kw.get('w_floor', 1.); pol.w_beta = kw.get('w_beta', 1.); pol.reset_hooks = True; return pol
CFG = {'Proposed': dict(), 'NoNoise': dict(), 'KPF0.5': dict(KPF1=.5), 'FF-noVel': dict(novel=True)}
if len(sys.argv) > 1 and sys.argv[1] == 'force': CFG = {'Proposed': dict(), 'SelF0.3': dict(sched='select', servo_crit=True, force_crit=True), 'SelF0.1': dict(sched='select', servo_crit=True, force_crit=True, sel_T=.1)}
elif len(sys.argv) > 1 and sys.argv[1] == 'calib':
    q0, m0 = F['q_smooth'], F['q_maneuver']
    CFG = {'Proposed': dict(), 'R2': dict(r_scale=2.), 'R2-qm15': dict(r_scale=2., q_maneuver=15*m0), 'R2-qs15-qm15': dict(r_scale=2., q_smooth=15*q0, q_maneuver=15*m0), 'R2-qm100': dict(r_scale=2., q_maneuver=100*m0)}
elif len(sys.argv) > 1 and sys.argv[1] == 'wiener': CFG = {'Wien-b5': dict(sched='wiener', w_beta=5.), 'Wien-b10': dict(sched='wiener', w_beta=10.), 'Wien-b5-T0.06': dict(sched='wiener', w_T=.06, w_beta=5.)}
elif len(sys.argv) > 1 and sys.argv[1] == 'shrink': CFG = {'Proposed': dict(), 'Shrink1': dict(sched='shrink'), 'Shrink4': dict(sched='shrink', shrink_c=4.)}
elif len(sys.argv) > 1 and sys.argv[1] == 'two': CFG = {'Proposed': dict(), 'Sel2-0.04': dict(sched='select', servo_crit=True, fast_T=.04), 'Sel2-0.02': dict(sched='select', servo_crit=True, fast_T=.02)}
elif len(sys.argv) > 1 and sys.argv[1] == 'servo': CFG = {'Proposed': dict(), 'SelServo0.3': dict(sched='select', servo_crit=True), 'SelServo1.0': dict(sched='select', servo_crit=True, sel_T=1.)}
elif len(sys.argv) > 1 and sys.argv[1] == 'bump': CFG = {'Sel0.3': dict(sched='select'), 'Sel0.3-bumpless': dict(sched='select', bumpless=True)}
elif len(sys.argv) > 1 and sys.argv[1] == 'select': CFG = {'Proposed': dict(), 'Sel0.3': dict(sched='select'), 'Sel1.0': dict(sched='select', sel_T=1.), 'Sel0.1': dict(sched='select', sel_T=.1)}
elif len(sys.argv) > 1 and sys.argv[1] == 'spatial': CFG = {'Proposed': dict(), 'Lead3-spatial': dict(vel_lp=3., sched='spatial'), 'Lead6-spatial': dict(vel_lp=6., sched='spatial')}
elif len(sys.argv) > 1 and sys.argv[1] == 'sched': CFG = {'Proposed': dict(), 'Lead3-oracle': dict(vel_lp=3., sched='oracle'), 'Lead3-est': dict(vel_lp=3., sched='est'), 'Lead6-est': dict(vel_lp=6., sched='est')}
elif len(sys.argv) > 1 and sys.argv[1] == 'lead': CFG = {'Proposed': dict(), 'Lead3Hz': dict(vel_lp=3.), 'Lead6Hz': dict(vel_lp=6.), 'Lead10Hz': dict(vel_lp=10.)}
elif len(sys.argv) > 1: CFG = {'Proposed': dict(), 'Notch0.8': dict(notch=.8), 'Notch0.9': dict(notch=.9)}
B = ((0, 5), (5, 10), (10, 23), (23, 62.5))
def _ep(a):
    cfg, (fam, sid, p, T) = a
    p = p.copy()
    if cfg == 'NoNoise': p[14] = 0.
    pol = mk(**CFG[cfg]); pol.p_scene = p; log, _ = common.run(pol, p, T, sid)
    m = common.task_metrics(log, 20.8)
    if m['diverged']: return cfg, fam, None
    e = log[log[:, 0] >= 1., 3]-P.FD; f, S = welch(e-e.mean(), fs=125., nperseg=256); df = f[1]-f[0]
    return cfg, fam, ([float(np.sqrt(np.sum(S[(f >= lo) & (f < hi)])*df)) for lo, hi in B], m['rmse'], m['peak'],
                      m['loss_max'], float(np.median(pol.mu_log)) if pol.mu_log else np.nan, getattr(pol, 'sel_log', None), float(np.median(pol.g_log)) if hasattr(pol, 'g_log') else np.nan)
if __name__ == '__main__':
    fams = ('S2', 'T3m', 'T2m-stiff', 'T2f-stiff')
    sc = [s for s in CE.cohort('diag', 20) if s[0] in fams]
    with Pool(16) as po: rows = po.map(_ep, [(c, s) for c in CFG for s in sc], chunksize=1)
    print('bands Hz', B)
    for fam in fams:
        for c in CFG:
            R = [r[2] for r in rows if r[0] == c and r[1] == fam]; ok = [x for x in R if x is not None]
            v = np.array([x[0] for x in ok]); rm = [x[1] for x in ok]
            print('%-9s %-9s rmse %.3f/%.3f r%d L%d mu_man %.2f | %s' % (fam, c, np.median(rm), np.mean(rm),
                  len(R)-len(ok)+sum(x[2] > 50 for x in ok), sum(x[3] > .5 for x in ok), np.nanmedian([x[4] for x in ok]),
                  ' '.join('%.3f' % b for b in np.median(v, 0))), flush=True)
            sl = [x[5]/x[5].sum() for x in ok if x[5] is not None and np.ndim(x[5]) == 1 and len(x[5]) == 4]
            if sl: print('    selection fractions (inf,10,6,3 Hz):', np.round(np.mean(sl, 0), 2))
            gg = [x[6] for x in ok if np.isfinite(x[6])]
            if gg: print('    median lead gain', np.round(np.median(gg), 2))
