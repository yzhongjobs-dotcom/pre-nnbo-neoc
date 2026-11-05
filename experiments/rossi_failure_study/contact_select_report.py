"""Tables for the contact-model selection, from the raw rows of results/contact_select.json."""
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
name = sys.argv[1] if len(sys.argv) > 1 else 'contact_select.json'
data = json.loads((HERE/'results'/name).read_text(encoding='utf-8'))
rows = data['rows'] if isinstance(data, dict) else data
REGIMES = ('nominal', 'stop_go', 'reversal', 'curved', 'rate', 'memory')
NAMES = [n for n in ('L', 'W', 'S', 'H', 'Hz', 'AS', 'TR', 'LF', 'WF', 'WA', 'WT', 'HA', 'HT') if any(r['candidate'] == n for r in rows)]
NAMES += ['P5:'+n for n in NAMES if any(r['candidate'] == 'P5:'+n for r in rows)]
L = lambda m: math.log1p(min(m['J'], 1000))
key = lambda r: (r['cond'], r['tissue'], r['i'])
boot = lambda d: np.percentile(d[np.random.default_rng(1).integers(0, len(d), (10000, len(d)))].mean(1), [2.5, 97.5])
fail = lambda m: bool(m['failed']) or m['peak'] > 50

print('n rows', len(rows))
print('\nmean L by regime, then mean over regimes, worst significant difference from L')
print('name     '+' '.join(f'{g:>9s}' for g in REGIMES)+'      mean')
means = {}
for n in NAMES:
    v = [float(np.mean([L(r['metrics']) for r in rows if r['regime'] == g and r['candidate'] == n])) for g in REGIMES]
    means[n] = float(np.mean(v))
    print(f'{n:8s} '+' '.join(f'{x:9.3f}' for x in v)+f' {np.mean(v):9.3f}')
print('\ndifference from L: mean [95%] per regime; * excludes 0')
for n in NAMES:
    if n in ('L', 'P5:L'):
        continue
    base = 'P5:L' if n.startswith('P5:') else 'L'
    cells = []
    for g in REGIMES:
        ref = {key(r): r['metrics'] for r in rows if r['regime'] == g and r['candidate'] == base}
        d = np.array([L(r['metrics'])-L(ref[key(r)]) for r in rows if r['regime'] == g and r['candidate'] == n]); c = boot(d)
        cells.append(f"{d.mean():+.3f}{'*' if c[0] > 0 or c[1] < 0 else ' '}")
    all_ref = {(r['regime'],)+key(r): r['metrics'] for r in rows if r['candidate'] == base}
    d = np.array([L(r['metrics'])-L(all_ref[(r['regime'],)+key(r)]) for r in rows if r['candidate'] == n]); c = boot(d)
    print(f"{n:8s} "+' '.join(f'{x:>8s}' for x in cells)+f"  all {d.mean():+.3f} [{c[0]:+.3f},{c[1]:+.3f}]")
print('\nall regimes pooled: L, RMSE Fn median/P90, RMSE F_axis median, event median, runaways, long losses')
for n in NAMES:
    rr = [r['metrics'] for r in rows if r['candidate'] == n]
    print(f"{n:8s} L {np.mean([L(m) for m in rr]):.3f} Fn {np.nanmedian([m['fn_rmse'] for m in rr]):.3f}/{np.nanpercentile([m['fn_rmse'] for m in rr], 90):.3f} axis {np.nanmedian([m['faxis_rmse'] for m in rr]):.3f} event {np.nanmedian([m['event_rmse'] for m in rr]):.3f} fail {sum(fail(m) for m in rr)} long {sum(m['long_loss'] for m in rr)} call {np.nanmedian([m['call_p50'] for m in rr]):.0f} us")
print('\ntissue groups pooled over regimes, mean L')
GROUPS = (('relaxed', lambda r: r['tissue'] == 'relaxed'), ('soft', lambda r: r['tissue'] == 'soft'), ('mid', lambda r: r['tissue'] == 'mid'), ('hard', lambda r: r['tissue'] == 'hard'), ('mixed', lambda r: r['tissue'] == 'mixed'),
          ('flinch+cough', lambda r: r['cond'] in ('flinch', 'cough')), ('fast', lambda r: r['cond'] == 'shift_fast'), ('smooth', lambda r: r['cond'] in ('breath', 'shift_slow')), ('rest', lambda r: r['cond'] == 'rest'))
print('name     '+' '.join(f'{g[:9]:>9s}' for g, _ in GROUPS))
for n in NAMES:
    print(f'{n:8s} '+' '.join(f"{np.mean([L(r['metrics']) for r in rows if r['candidate'] == n and f(r)]):9.3f}" for _, f in GROUPS))
print('\nranking by mean over regimes', sorted(means, key=means.get))
