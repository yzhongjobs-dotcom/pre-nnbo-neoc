"""Model-free readings of the Neumann records, to be set against what each mechanism says should happen. Nothing is fitted across protocols here."""
import json
from pathlib import Path

import numpy as np

import neumann_dynamics as D

HERE = Path(__file__).resolve().parent


def slope(x, y):
    return float(np.polyfit(np.log(x), np.log(y), 1)[0])


def main():
    data = list(D.prepare())
    out = {}
    print('DMA per frequency, mean over depth: f, phase deg, tan, local slope of the magnitude, the phase that slope implies, gap')
    print('Causality ties phase to the slope of the magnitude (Bode), locally phase = 90 deg * slope. A time skew between the channels adds 360 f dt to the phase and nothing to the magnitude.')
    dma = [r for r in data if r['kind'] == 'dma']
    out['dma'] = {}
    for tissue in ('white', 'grey'):
        rs = [r for r in dma if r['tissue'] == tissue]
        f = np.array([x['f'] for x in rs[0]['rows']])
        K = np.array([[x['K'] for x in r['rows']] for r in rs])
        ph = np.degrees(np.angle(K))
        grad = np.gradient(np.log(np.abs(K)), np.log(f), axis=1)
        table = []
        for i, ff in enumerate(f):
            local = float(grad[:, i].mean())
            table.append(dict(f=float(ff), phase=float(ph[:, i].mean()), phase_sd=float(ph[:, i].std()), tan=float(np.tan(np.radians(ph[:, i])).mean()), local_slope=local, bode_phase=90.*local, gap=float(ph[:, i].mean()-90.*local)))
            print(f"  {tissue[:5]} {ff:5.1f}  {table[-1]['phase']:6.2f} sd {table[-1]['phase_sd']:.2f}  tan {table[-1]['tan']:.3f}  slope {local:.3f}  90*slope {90.*local:5.1f}  gap {table[-1]['gap']:+5.1f}")
        gap = np.array([x['gap'] for x in table])
        skew = float(np.sum(gap*360.*f)/np.sum((360.*f)**2))
        print(f"  {tissue[:5]} time skew the gap would imply: {skew*1e3:+.3f} ms, {skew*360.*20.8:+.1f} deg at 20.8 Hz")
        out['dma'][tissue] = dict(table=table, skew_ms=skew*1e3)

    print('\nDMA against depth. Hunt-Crossley: tan grows with depth to the power 1. Linear viscoelastic contact: power 0. Stiffness grows with depth to the power n-1')
    out['depth'] = {}
    for tissue in ('white', 'grey'):
        rs = [r for r in dma if r['tissue'] == tissue and r['name'] != 'PSC1501 Indentation_032']
        d0 = np.array([r['rows'][0]['d0'] for r in rs])
        row = {}
        for i, f in enumerate([x['f'] for x in rs[0]['rows']]):
            t = np.array([r['rows'][i]['K'].imag/r['rows'][i]['K'].real for r in rs]); k = np.array([abs(r['rows'][i]['K']) for r in rs])
            row[str(f)] = dict(tan_vs_depth=slope(d0, t), stiffness_vs_depth=slope(d0, k))
            print(f"  {tissue[:5]} {f:5.1f} Hz, depths {np.round(d0,1)}: tan {np.round(t,3)} exponent {row[str(f)]['tan_vs_depth']:+.2f}; |K| exponent {row[str(f)]['stiffness_vs_depth']:.2f}")
        out['depth'][tissue] = row

    print('\nLoading rate: force at one depth against speed. Same site, fastest last, so preconditioning works against the effect')
    out['rate'] = {}
    for tissue in ('white', 'grey'):
        rs = sorted([r for r in data if r['kind'] == 'rate' and r['tissue'] == tissue], key=lambda r: r['speed'])
        v = np.array([r['speed'] for r in rs])
        depth = .85*min(r['depth_peak'] for r in rs)
        F = []; plateau = []
        for r in rs:
            up = r['t'] <= r['t_peak']
            F.append(float(np.interp(depth, r['d'][up], r['f'][up])))
            hold = (r['t'] > r['t_drop']-.5) & (r['t'] < r['t_drop']-.05)
            plateau.append(float(np.mean(r['f'][hold])/np.mean(r['d'][hold])**1.5))
        F = np.array(F); plateau = np.array(plateau)
        lin = np.polyfit(v, F, 1)
        linear_err = float(np.sqrt(np.mean((np.polyval(lin, v)-F)**2))/F.mean()); power = slope(v, F); pw = np.exp(np.polyval(np.polyfit(np.log(v), np.log(F), 1), np.log(v)))
        power_err = float(np.sqrt(np.mean((pw-F)**2))/F.mean())
        print(f"  {tissue[:5]} speeds {np.round(v,0)}\n        F at {depth:.0f} um {np.round(F,3)}\n        slope {power:.3f}; straight line in speed misses by {100*linear_err:.1f}%, power law by {100*power_err:.1f}%; after the 10 s hold, force/d^1.5 max/min {plateau.max()/plateau.min():.3f}, against {F.max()/F.min():.2f} at peak")
        out['rate'][tissue] = dict(speed=v.tolist(), force=F.tolist(), slope=power, linear_error_pct=100*linear_err, power_error_pct=100*power_err, plateau_ratio=float(plateau.max()/plateau.min()))

    print('\nRelaxation, per file: fraction left after the ramp ends at .1, 1, 10, 100, 500 s; local slope between')
    out['relax'] = {}
    pts = np.array([.1, 1., 10., 100., 500.])
    for tissue in ('white', 'grey'):
        rows = []
        for r in [r for r in data if r['kind'] == 'relax' and r['tissue'] == tissue]:
            t = r['t']-r['t_peak']; f = r['f']; q = t > 0.
            f0 = f[q][np.argmax(t[q] > .02)]
            val = np.array([np.mean(f[(t > p*.92) & (t < p*1.08)])/f0 for p in pts])
            rows.append(val)
            print(f"  {tissue[:5]} {r['specimen'][:8]:8s} {r['name'][-3:]} peak {r['peak_force']:5.2f}  {np.round(val,3)}  slopes {np.round(np.diff(np.log(val))/np.diff(np.log(pts)),3)}")
        rows = np.array(rows)
        out['relax'][tissue] = dict(times=pts.tolist(), median=np.median(rows, axis=0).tolist(), low=rows.min(axis=0).tolist(), high=rows.max(axis=0).tolist())

    print('\nCyclic, per file: peak per cycle relative to the first; dissipated fraction of the loop; force on the way back at half depth relative to the way up; fraction of unloading with almost no force')
    out['cyclic'] = {}
    T = 10.036
    for tissue in ('white', 'grey'):
        rows = []
        for r in [r for r in data if r['kind'] == 'cyclic' and r['tissue'] == tissue]:
            t, d, f = r['t'], r['d'], r['f']
            t1 = r['t_peak']
            peaks = []
            for c in range(5):
                tp = t1+c*T; q = (t > tp-1.5) & (t < tp+1.5)
                if q.sum() > 3:
                    peaks.append(float(f[q].max()))
            ratio = np.array(peaks)/peaks[0]
            hyst = []; back = []; lost = []
            for c in range(4):
                tp = t1+c*T; up = (t >= tp-T/2.) & (t <= tp); dn = (t > tp) & (t <= tp+T/2.)
                if up.sum() < 5 or dn.sum() < 5:
                    continue
                o = np.argsort(d[dn])
                e_up = np.trapz(f[up], d[up]); e_dn = -np.trapz(f[dn], d[dn])
                hyst.append(float(1.-e_dn/e_up))
                half = .5*d[up].max()
                back.append(float(np.interp(half, d[dn][o], f[dn][o])/np.interp(half, d[up], f[up])))
                lost.append(float(np.mean(f[dn] < .02*peaks[0])))
            rows.append(dict(name=r['name'], peak_ratio=ratio.tolist(), dissipation=hyst, back_ratio=back, no_force_fraction=lost))
            print(f"  {tissue[:5]} {r['name'][-3:]} peak ratio {np.round(ratio,3)}  dissipated {np.round(hyst,3)}  back/up at half depth {np.round(back,3)}  no force {np.round(lost,2)}")
        out['cyclic'][tissue] = rows
    (HERE/'results'/'neumann_tests.json').write_text(json.dumps(out, indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
