"""Are all scenes valid for the Method? Not a comparison: one repetition of every cell and regime for two forms on two hosts."""
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import method3 as M
import scenes3 as S3

HERE = Path(__file__).resolve().parent
NAMES = ('L',)        # validity of the scenes only, not a comparison


def _job(a):
    regime, (cond, tissue), i = a
    sc = S3.cell_scene3(cond, tissue, regime, i, base=440000)
    out = []
    for name in NAMES:
        log, pol, sec = M.run(name, sc, seed=i+7)
        m = M.metric(log, pol, sc)
        m.update(name=name, regime=regime, cond=cond, tissue=tissue, i=i, sec=sec, force=sc['FD'], probe=sc['probe'], law=sc['contact_law'], slide=sc['slide'], phi_hc=sc['phi_hc'], cover=sc['cover'])
        if log is not None:
            m['d_peak'] = float(np.max(np.abs(np.diff(log[:, 6]))/.008))
        out.append(m)
    return out


def main():
    jobs = [(g, c, 0) for g in S3.REGIMES for c in S3.CELLS3]
    with ProcessPoolExecutor(max_workers=15) as ex:
        rows = [r for b in ex.map(_job, jobs, chunksize=1) for r in b]
    (HERE/'results'/'smoke3.json').write_text(json.dumps(rows, indent=1), encoding='utf-8')
    fail = lambda r: r['failed'] or r['peak'] > 50
    print('rows', len(rows), 'median seconds', np.median([r['sec'] for r in rows]))
    for key in ('regime', 'cond', 'tissue', 'force', 'probe', 'law', 'slide', 'phi_hc'):
        print('\nby', key)
        for v in sorted({r[key] for r in rows}, key=str):
            line = []
            for n in NAMES:
                q = [r for r in rows if r[key] == v and r['name'] == n]
                line.append(f"{n}: L {np.mean([np.log1p(r['J']) for r in q]):.2f} rmse {np.nanmedian([r['rmse'] for r in q]):.2f} fail {sum(fail(r) for r in q)}/{len(q)} long {sum(r['long_loss'] for r in q)}")
            print(f'{str(v):16s}', ' | '.join(line))
    print('\nscenes where all four run away')
    keys = {(r['regime'], r['cond'], r['tissue']) for r in rows}
    n = 0
    for k in sorted(keys):
        q = [r for r in rows if (r['regime'], r['cond'], r['tissue']) == k]
        if all(fail(r) for r in q):
            n += 1; print(k, q[0]['force'], q[0]['probe'], q[0]['law'], q[0]['cover'], [round(r['peak'], 1) for r in q])
    print('count', n, 'of', len(keys))
    print('normal speed peak mm/s: median', np.nanmedian([r.get('d_peak', np.nan) for r in rows])*1000., 'P90', np.nanpercentile([r.get('d_peak', np.nan) for r in rows], 90)*1000.)


if __name__ == '__main__':
    main()
