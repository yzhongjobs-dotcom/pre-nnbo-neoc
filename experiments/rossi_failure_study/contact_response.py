"""What gain and phase the probe actually sees: frequency response from probe motion to contact force at a frozen state.

Noise-free replay with no tissue motion: commanded sine on the encoder, F_true regressed on it. Repeated with a rigid
support (very large body mass) to separate tissue memory from the body yielding, and at three working forces for the
static nonlinearity. Checks whether the phase at the probe frequency predicts the gain at lower frequencies, as a weak
power law with causality would: |G(w)| = |G(wc)| (w/wc)^beta, phase = beta pi/2.
"""
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import physics_plant as PP
import scenes as SC

HERE = Path(__file__).resolve().parent
FREQS = (.5, 1., 2., 4., 8., 12., 16., 20.8, 26., 35., 45.)
TC = .008


def scene(tissue, seed, rigid):
    rng = np.random.default_rng(seed)
    over = dict(field_sigma=0.)
    if rigid:
        over['m'] = 1e4
    sc = SC.make('rest', tissue, rng, T=12., v=1e-6, overrides=over)
    sc['h_s'] = np.zeros_like(sc['h_s'])
    return sc


def response(sc, y0, f, amp):
    n = int(round(sc['T']/TC)); t = np.arange(n)*TC
    log = PP.replay(sc, y0+amp*np.sin(2*np.pi*f*t))
    if log is None:
        return None
    m = log[:, 0] >= sc['T']-max(6., 6./f)
    if m.sum() < 8:
        m = log[:, 0] >= 2.
    tt = log[m, 0]; a = np.c_[np.sin(2*np.pi*f*tt), np.cos(2*np.pi*f*tt), np.ones(m.sum()), tt-tt.mean()]
    cy = np.linalg.lstsq(a, log[m, 2], rcond=None)[0]; cf = np.linalg.lstsq(a, log[m, 3], rcond=None)[0]
    g = complex(cf[0], cf[1])/complex(cy[0], cy[1])
    return abs(g), float(np.angle(g)), float(np.mean(log[m, 3])), float(np.mean(log[m, 8])), float(np.mean(log[m, 9])), float(np.mean(log[m, 10]))


def job(a):
    tissue, seed, rigid, dy = a
    sc = scene(tissue, seed, rigid)
    s, _ = PP._init(sc, PP.PLATFORM)
    y0 = s[2]+dy
    settle = response(sc, y0, 20.8, 1e-7)
    if settle is None:
        return None
    amp = .15/max(settle[0], 1.)
    rows = []
    for f in FREQS:
        r = response(sc, y0, f, amp)
        if r is not None:
            rows.append(dict(f=f, gain=r[0], phase=r[1], F=r[2], k8=r[3], k9=r[4], k10=r[5]))
    return dict(tissue=tissue, seed=seed, rigid=rigid, dy=dy, amp=amp, rows=rows)


def main():
    args = []
    for tissue in ('soft', 'mid', 'hard'):
        for seed in range(3):
            for rigid in (False, True):
                args.append((tissue, 600+seed, rigid, 0.))
            # working point changes, fixed by the static curve: small offsets of the probe
            for frac in (-.4, .4, 1.):
                sc = scene(tissue, 600+seed, False); s, _ = PP._init(sc, PP.PLATFORM)
                k0 = response(sc, s[2], 20.8, 1e-7)[0]
                args.append((tissue, 600+seed, False, frac*4./k0))
    with ProcessPoolExecutor(max_workers=15) as pool:
        out = [r for r in pool.map(job, args, chunksize=1) if r]
    (HERE/'results'/'contact_response.json').write_text(json.dumps(out, indent=1), encoding='utf-8')
    for tissue in ('soft', 'mid', 'hard'):
        for rigid in (False, True):
            rr = [r for r in out if r['tissue'] == tissue and r['rigid'] == rigid and r['dy'] == 0.]
            print(f'\n{tissue} rigid={rigid} n={len(rr)} F={np.mean([r["rows"][7]["F"] for r in rr]):.2f} k_probe_col8={np.mean([r["rows"][7]["k8"] for r in rr]):.0f}')
            ref = np.array([r['rows'][7]['gain'] for r in rr])
            phc = np.array([r['rows'][7]['phase'] for r in rr]); beta = 2*phc/np.pi
            for j, f in enumerate(FREQS):
                g = np.array([r['rows'][j]['gain'] for r in rr])/ref; ph = np.array([r['rows'][j]['phase'] for r in rr])
                pred = (f/20.8)**beta
                print(f'  {f:5.1f} Hz gain/gain(20.8) {g.mean():.3f}  phase {np.degrees(ph.mean()):+6.1f} deg  power-law pred {pred.mean():.3f}')
            print(f'  gain(20.8) {ref.mean():.0f}, beta from phase {beta.mean():.3f}')
    print('\nstatic nonlinearity: ln k vs ln F at 20.8 Hz')
    for tissue in ('soft', 'mid', 'hard'):
        pts = [(r['rows'][7]['F'], r['rows'][7]['gain']) for r in out if r['tissue'] == tissue and not r['rigid']]
        x, y = np.log(np.array(pts)).T
        print(f'  {tissue}: gamma {np.polyfit(x, y, 1)[0]:.3f}, F range {np.exp(x.min()):.1f}-{np.exp(x.max()):.1f} N')


if __name__ == '__main__':
    main()
