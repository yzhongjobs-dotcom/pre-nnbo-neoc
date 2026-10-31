"""Per-tick trace of one scene for failure analysis: force, true k, k_bar, k_good, latest window estimate and its
precision/consistency flags, IMM maneuver probability, supervisor mode. Prints a decimated table around the worst peak.
Run: python -X utf8 trace_scene.py <family> <sid> [cohort] [config]
"""
import json, os, sys
import numpy as np
import common
import scenes_v2 as SV
import cohort_eval as CE
import variants

H = common.HERE
BASE = json.load(open(os.path.join(H, 'tuned_v2', 'Proposed-IMM.json')))['best']['params']


class Traced(variants.MotionAwareIMM):
    def reset(self, platform):
        super().reset(platform)
        self.tr = {}

    def step(self, o):
        y = super().step(o)
        f = self.fit or {}
        mu = float(self.tracker.mu[1]) if self.tracker.initialized else np.nan
        self.tr[o.n] = (self.kbar, self.k_good, self.k_hat, self.rel, f.get('chi2', np.nan), mu, self.mode)
        return y


if __name__ == '__main__':
    fam, sid = sys.argv[1], int(sys.argv[2])
    coh = sys.argv[3] if len(sys.argv) > 3 else 'diag'
    cfg = sys.argv[4] if len(sys.argv) > 4 else 'final'
    sc = [s for s in CE.cohort(coh, 60) if s[0] == fam and s[1] == sid][0]
    extra = dict(CE.extra_configs()).get(cfg, {})
    pol = Traced(single_model=False, **dict(BASE, **extra))
    log, pol = common.run(pol, sc[2], sc[3], sid)
    t, F, k = log[:, 0], log[:, 3], log[:, 8]
    i0 = int(np.argmax(np.where(t > 1, F, 0)))
    first = int(np.argmax((t > 1) & (F > 12)))
    if os.environ.get('FOCUS') == 'loss':
        lost = (t > 1) & (F < 0.5); best = (0, 0); i = 0
        while i < len(t):
            if lost[i]:
                j = i
                while j < len(t) and lost[j]: j += 1
                best = max(best, (j-i, i)); i = j
            else: i += 1
        first = best[1]; print('longest loss %.2f s from t=%.2f' % (best[0]*0.008, t[first]))
    print('worst peak %.1f N at t=%.2f; first >12 N at t=%.2f; events:' % (F[i0], t[i0], t[first]),
          [(round(sc[2][16+4*j], 2), round(sc[2][18+4*j], 2), round(sc[2][19+4*j], 2)) for j in range(int(sc[2][15]))])
    print('   t     F      y        surf     k_true   k_bar    k_good   k_hat   rel   chi2  mu   mode')
    for i in range(max(0, first-100), min(len(t), first+int(os.environ.get('AFTER', '60'))), 8):
        kb, kg, kh, rel, chi, mu, md = pol.tr.get(i, (np.nan,)*7)
        import plant as Pl
        sf = Pl.surface(log[i, 1], t[i], sc[2])[0]
        print('%6.2f %6.1f %8.2f %8.2f %8.0f %8.0f %8.0f %8.0f %5.2f %5.1f %4.2f %d' % (t[i], F[i], log[i, 2]*1e3, sf*1e3, k[i], kb, kg, kh, rel, chi, mu, md))
