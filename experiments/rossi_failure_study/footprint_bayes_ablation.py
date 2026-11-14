"""Which prior knowledge carries the gain of the exact-generator filter (footprint_bayes_synthetic.py), and how fragile is it?

Truth: rib/gap density (contrast 1.6, level jitter .15, rib 40-90 ticks, gap 20-70 ticks), footprint M = 60 ticks, per-tick noise as before.
The particle filter's model is varied one item at a time; then the truth is varied with a matched and with the nominal model.
Run: python -X utf8 footprint_bayes_ablation.py
"""
from multiprocessing import Pool
import numpy as np
from footprint_prior_synthetic import TC, make, window

NOM = dict(M=60, contrast=1.6, jit=.15, rib=(40, 90), gap=(20, 70), geom=False)
R = 24*.17**2; N = 24000; NP = 4000


def new_segment(typ, m, rng):
    if m['geom']:                                        # memoryless lengths with the same means: no knowledge of rib width
        mean = np.where(typ == 1, np.mean(m['rib']), np.mean(m['gap'])); ln = rng.geometric(1./mean)
    else:
        ln = np.where(typ == 1, rng.integers(*m['rib'], size=typ.size), rng.integers(*m['gap'], size=typ.size))
    return ln, typ*m['contrast']+rng.normal(0., m['jit'], size=typ.size)


def pf(obs, m, rng, npart=NP):
    M = m['M']; typ = rng.integers(0, 2, size=npart); ln, lev = new_segment(typ, m, rng); age = (rng.random(npart)*ln).astype(int)
    buf = np.zeros((npart, M)); ptr = 0; logw = np.zeros(npart); out = np.empty(len(obs)); idx = np.arange(npart)
    for j in range(-3*M, len(obs)):
        done = age >= ln
        if done.any():
            typ[done] = 1-typ[done]; ln[done], lev[done] = new_segment(typ[done], m, rng); age[done] = 0
        buf[:, ptr] = lev; ptr = (ptr+1) % M; age += 1
        if j < 0:
            continue
        k = buf.mean(1); logw += -.5*(obs[j]-k)**2/R; logw -= logw.max(); w = np.exp(logw); w /= w.sum(); out[j] = w@k
        if 1./np.sum(w*w) < .5*npart:
            pick = np.searchsorted(np.cumsum(w), (rng.random()+idx)/npart).clip(0, npart-1)
            typ, ln, lev, age, buf = typ[pick], ln[pick], lev[pick], age[pick], buf[pick]; logw[:] = 0.
    return out


def run(job):
    name, truth, model = job
    t = {**NOM, **truth}; m = {**NOM, **model}; rng = np.random.default_rng(3)
    e, k = make(N, t['M'], rng, rib=t['rib'], gap=t['gap'], contrast=t['contrast'], jitter=t['jit'])
    obs = k+rng.normal(size=N)*np.sqrt(R); s = slice(500, None)
    dk = np.gradient(k, TC); thr = .4*t['contrast']/(t['M']*TC); rise = (dk > thr)[s]; fall = (dk < -thr)[s]
    ew = (window(obs)-k)[s]; eb = (pf(obs, m, np.random.default_rng(11))-k)[s]; rms = lambda x: float(np.sqrt(np.nanmean(x*x)))
    return (f'{name:44s} | {rms(ew):.3f} {rms(eb):.3f} | window {np.nanmean(ew[rise]):+.3f} {np.nanmean(ew[fall]):+.3f}  filter {eb[rise].mean():+.3f} {eb[fall].mean():+.3f}'
            f' | worst 1% under-estimate: window {np.nanpercentile(ew, 1):+.2f} filter {np.percentile(eb, 1):+.2f}')


JOBS = [
    ('exact model', {}, {}),
    ('model: level jitter .4 (levels loosely known)', {}, dict(jit=.4)),
    ('model: contrast 1.2 (truth 1.6)', {}, dict(contrast=1.2)),
    ('model: contrast 2.0 (truth 1.6)', {}, dict(contrast=2.0)),
    ('model: footprint 45 ticks (truth 60)', {}, dict(M=45)),
    ('model: footprint 75 ticks (truth 60)', {}, dict(M=75)),
    ('model: memoryless segment lengths', {}, dict(geom=True)),
    ('model: jitter .4 + memoryless lengths', {}, dict(jit=.4, geom=True)),
    ('truth jitter .4, matched model', dict(jit=.4), dict(jit=.4)),
    ('truth jitter .4, nominal model', dict(jit=.4), {}),
    ('truth contrast .8, matched model', dict(contrast=.8), dict(contrast=.8)),
    ('truth contrast .8, nominal model (1.6)', dict(contrast=.8), {}),
    ('truth contrast .8, model jitter .4 + memoryless', dict(contrast=.8), dict(contrast=.8, jit=.4, geom=True)),
    ('truth footprint 30 ticks (fast scan), matched', dict(M=30), dict(M=30)),
]

if __name__ == '__main__':
    print('case | rms of ln k error: window, filter | mean error on rises, falls | 1st percentile of error')
    with Pool(7) as p:
        for line in p.imap(run, JOBS):
            print(line, flush=True)
