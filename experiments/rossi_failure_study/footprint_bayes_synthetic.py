"""Near Bayes-optimal causal estimate for the synthetic rib/gap footprint data of footprint_prior_synthetic.py.

Question: the minimax bound in research/v2_representation_adjudication_20261114.md is for a Lipschitz class, where linear estimators are
near optimal. The synthetic profile is a structured class (two levels, sharp edges, known footprint), where a nonlinear estimator that
knows the generator exactly could do much better. A particle filter with the exact generator as its model gives the achievable error of
ANY causal estimator on this data (up to particle noise). Compare with the 24-tick window and the tuned random-walk memory.
Run: python -X utf8 footprint_bayes_synthetic.py [particles]
"""
import sys
import numpy as np
from footprint_prior_synthetic import TC, make, window, rw

RIB, GAP, CONTRAST, JIT = (40, 90), (20, 70), 1.6, .15


def new_segment(typ, rng):
    """Lengths and levels for particles that start a segment of type `typ` (0 gap, 1 rib)."""
    ln = np.where(typ == 1, rng.integers(*RIB, size=typ.size), rng.integers(*GAP, size=typ.size))
    return ln, typ*CONTRAST+rng.normal(0., JIT, size=typ.size)


def particle_filter(obs, R, M, npart, rng):
    typ = rng.integers(0, 2, size=npart); ln, lev = new_segment(typ, rng); age = (rng.random(npart)*ln).astype(int)
    buf = np.zeros((npart, M)); ptr = 0; logw = np.zeros(npart); out = np.empty(len(obs)); idx = np.arange(npart)
    for j in range(-3*M, len(obs)):                      # the first 3M ticks only fill the buffers from the prior
        done = age >= ln
        if done.any():
            typ[done] = 1-typ[done]; ln[done], lev[done] = new_segment(typ[done], rng); age[done] = 0
        buf[:, ptr] = lev; ptr = (ptr+1) % M; age += 1
        if j < 0:
            continue
        k = buf.mean(1); logw += -.5*(obs[j]-k)**2/R; logw -= logw.max(); w = np.exp(logw); w /= w.sum()
        out[j] = w@k
        if 1./np.sum(w*w) < .5*npart:                    # systematic resampling
            pick = np.searchsorted(np.cumsum(w), (rng.random()+idx)/npart).clip(0, npart-1)
            typ, ln, lev, age, buf = typ[pick], ln[pick], lev[pick], age[pick], buf[pick]; logw[:] = 0.
    return out


if __name__ == '__main__':
    npart = int(sys.argv[1]) if len(sys.argv) > 1 else 4000
    rng = np.random.default_rng(3); prng = np.random.default_rng(11); R = 24*.17**2; n = 40000; s = slice(500, None)
    print(f'particles {npart}. M | rms of ln k error: window / best random-walk / Bayes (exact generator) | mean error on rises, falls: window, Bayes | rms on rises+falls: window, Bayes | plateau rms: window, Bayes')
    for M in (60, 75, 120):
        e, k = make(n, M, rng); obs = k+rng.normal(size=n)*np.sqrt(R)
        dk = np.gradient(k, TC); rise = (dk > 1.)[s]; fall = (dk < -1.)[s]; flat = ~(rise | fall)
        w = window(obs)
        r = min((rw(obs, R, q) for q in (1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3)), key=lambda est: np.mean((est-k)[s]**2))
        b = particle_filter(obs, R, M, npart, prng)
        ew, er, eb = (w-k)[s], (r-k)[s], (b-k)[s]; rms = lambda x: float(np.sqrt(np.nanmean(x*x)))
        print(f'{M:4d} ({M*TC:.2f} s) | {rms(ew):.3f} / {rms(er):.3f} / {rms(eb):.3f} | window {np.nanmean(ew[rise]):+.3f} {np.nanmean(ew[fall]):+.3f}, Bayes {eb[rise].mean():+.3f} {eb[fall].mean():+.3f}'
              f' | {rms(ew[rise | fall]):.3f} {rms(eb[rise | fall]):.3f} | {rms(ew[flat]):.3f} {rms(eb[flat]):.3f}')
