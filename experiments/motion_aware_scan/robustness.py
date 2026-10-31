"""Robustness to model mismatch and out-of-distribution conditions (paper draft §五, item: assumptions violated).

Mismatch: every controller is reset with the NOMINAL platform (8 ms, one-period delay, 19 ms servo lag, 30 Hz force
filter) while the simulated arm runs a perturbed one -- servo lag x0.6 / x1.5, one extra period of command delay,
20 Hz force filter -- or the controller reads force with extra white noise (0.08 N, doubling the variance).
Out of distribution (scene parameters outside the tuning ranges): fast/large breathing (0.4-0.6 Hz, 5-8 mm) on the
spine scan; ridge ratio 40-60x on the static scan; large (8-12 mm) and fast (0.10-0.15 s) pushes on stiff tissue.
Scenes: 30 per family from the untouched final cohort seeds (130000+) for mismatch; OOD scenes use seeds 150000+.
Run: python -X utf8 robustness.py -> results_v2/robustness.json
"""
import json, math, os
import numpy as np
from multiprocessing import Pool
import common
import cohort_eval as CE
import scenes_v2 as SV
from common import P

H = common.HERE
RES = os.path.join(H, 'results_v2')
NOMINAL = dict(P.PLATFORMS['PC125'])
METHODS = tuple(os.environ.get('METHODS', 'Proposed,RossiType,WijType,DOMPAC,YuenType').split(','))
MISMATCH = {'nominal': {}, 'lag x0.6': dict(tau_a=0.019*0.6), 'lag x1.5': dict(tau_a=0.019*1.5),
            'delay +1': dict(ndel=2), 'filter 20 Hz': dict(fc=20.), 'force noise x2': 'noise'}


class NominalModel(common.DirectPolicy):
    """controller sees the nominal platform (and optionally extra force noise); the arm runs the perturbed one."""
    def __init__(self, inner, extra_noise=0., seed=0):
        self.inner, self.extra_noise, self.seed = inner, extra_noise, seed
        self.name = inner.name

    def reset(self, platform):
        super().reset(NOMINAL)
        self.inner.reset(NOMINAL)
        self.rng = np.random.default_rng(self.seed); self.nz = None

    def step(self, o):
        if self.extra_noise <= 0:
            return self.inner.step(o)
        n = o.n
        if self.nz is None or len(self.nz) < len(o.log):
            # same spectrum as the native sensor noise: 1 kHz white noise through the 30 Hz first-order filter,
            # sampled at the control tick (adding it doubles the noise variance the controller sees)
            ns = int(round(NOMINAL['Tc']/1e-3)); afc = 1.-math.exp(-2*math.pi*NOMINAL['fc']*1e-3)
            w = self.rng.standard_normal((len(o.log)+10)*ns)*self.extra_noise; f = 0.; out = np.empty(len(o.log)+10)
            for i in range(len(out)):
                for _ in range(ns): f += afc*(w[i*ns+_]-f)
                out[i] = f
            self.nz = out
        view = o.log.copy(); view[:n, 5] += self.nz[:n]
        return self.inner.step(type(o)(view, n))


def ood_scenes(n=30):
    out = []
    rng = np.random.default_rng(150000)
    for i in range(n):
        p, T = SV.S.make('T3', rng, dict(ramp=float(rng.uniform(.3, .5)), A=float(rng.uniform(.005, .008))))
        p[3] = rng.uniform(.4, .6); out.append(('OOD-breathing', 150000+i, p, T))
    rng = np.random.default_rng(151000)
    for i in range(n):
        p, T = SV.S.make('T1', rng, dict(A=0.0, ratio=float(rng.uniform(40, 60)))); out.append(('OOD-ridge60', 151000+i, p, T))
    rng = np.random.default_rng(152000)
    for i in range(n):
        p, T = SV.S.make('T2', rng, dict(stiff=True, ramp=float(rng.uniform(.10, .15)), amp=float(rng.uniform(.008, .012))))
        out.append(('OOD-bigpush', 152000+i, p, T))
    return out


def _ep(a):
    cond, cfg, (fam, sid, p, T) = a
    spec = MISMATCH.get(cond, {})
    plat = dict(NOMINAL)
    if isinstance(spec, dict): plat.update(spec)
    try:
        pol = NominalModel(CE.build(cfg, p), extra_noise=0.08 if spec == 'noise' else 0., seed=sid)
        log = P.episode(p, pol, T=T, platform=plat, seed=sid)
        m = common.task_metrics(log, 20.8)
    except Exception as ex:
        return dict(cond=cond, cfg=cfg, family=fam, sid=sid, ok=False, err=str(ex))
    return dict(cond=cond, cfg=cfg, family=fam, sid=sid, ok=not m['diverged'], rmse=m['rmse'], peak=m['peak'],
                loss=m['loss'], loss_max=m['loss_max'])


if __name__ == '__main__':
    n = int(os.environ.get('N', '30'))
    fams = ('S2', 'T3m', 'T2m-stiff')
    coh = os.environ.get('COHORT', 'final')
    base = [s for s in CE.cohort(coh, n) if s[0] in fams]
    jobs = [(c, m, s) for c in MISMATCH for m in METHODS for s in base]
    if os.environ.get('OOD', '1') == '1':
        jobs += [('ood', m, s) for m in METHODS for s in ood_scenes(n)]
    with Pool(16) as po:
        rows = po.map(_ep, jobs, chunksize=1)
    json.dump(rows, open(os.path.join(RES, 'robustness_%s.json' % os.environ.get('TAG', coh)), 'w'), indent=1)
    for cond in list(MISMATCH)+['ood']:
        famset = sorted({r['family'] for r in rows if r['cond'] == cond})
        for fam in famset:
            line = '%-15s %-14s' % (cond, fam)
            for m in METHODS:
                R = [r for r in rows if r['cond'] == cond and r['cfg'] == m and r['family'] == fam]
                ok = [r for r in R if r['ok']]
                med = np.median([r['rmse'] for r in ok]) if ok else np.nan
                bad = sum((not r['ok']) or r['peak'] > 50 for r in R)
                line += ' | %s %.3f (%d bad)' % (m[:6], med, bad)
            print(line, flush=True)
