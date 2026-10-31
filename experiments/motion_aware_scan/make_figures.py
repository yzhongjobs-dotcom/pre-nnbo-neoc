"""Paper figures (research/figures/*.png).

Fig. 1 failure mechanism on one diagnosis scene (T2m-stiff, sid 113013): contact force and the measurement gain vs the true
       stiffness, earlier version (no acceptance rule, model-reconstructed displacement) vs final method.
Fig. 2 main comparison on the untouched confirmation cohort (cohort_final.json): per-scene RMSE, five families.
Fig. 3 dose-response sweeps (sweep.json): median RMSE and catastrophic scenes vs push ramp time and vs ridge ratio.
Palette: reference categorical slots 1-5 in fixed order (validated; contrast warning -> distinct markers + legends).
Run: python -X utf8 make_figures.py
"""
import json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import common
import cohort_eval as CE
from trace_scene import Traced

H = common.HERE
OUT = os.path.join(H, '..', '..', 'research', 'figures'); os.makedirs(OUT, exist_ok=True)
METH = ['Proposed', 'RossiType', 'WijType', 'DOMPAC', 'YuenType']
LAB = {'Proposed': 'Proposed', 'RossiType': 'Rossi-type', 'WijType': 'Wijayarathne-type (true k)', 'DOMPAC': 'DO-MPAC',
       'YuenType': 'Yuen-type (motion sensor)'}
COL = dict(zip(METH, ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4']))
MRK = dict(zip(METH, ['o', 's', '^', 'D', 'v']))
INK, MUTED, GRID = '#0b0b0b', '#52514e', '#e4e3df'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'figure.facecolor': '#fcfcfb',
                     'axes.facecolor': '#fcfcfb', 'savefig.facecolor': '#fcfcfb', 'lines.linewidth': 1.6})
FAMS = ['S2', 'T3m', 'T2m-soft', 'T2m-stiff', 'T2f-stiff']
FLAB = {'S2': 'S2 static span', 'T3m': 'S3a spine', 'T2m-soft': 'S3b soft push', 'T2m-stiff': 'S3c stiff push',
        'T2f-stiff': 'S3d fast stiff push'}


def fig1():
    sc = [s for s in CE.cohort('diag', 60) if s[0] == 'T2m-stiff' and s[1] == 113013][0]
    runs = {}
    for name, params in (('earlier version', CE.BASE), ('final method', CE.FINAL)):
        pol = Traced(single_model=False, **dict(params))
        log, pol = common.run(pol, sc[2], sc[3], sc[1])
        runs[name] = (log, pol)
    fig, ax = plt.subplots(2, 1, figsize=(6.4, 4.2), sharex=True)
    t0, t1 = 3.5, 6.3
    for name, c in (('earlier version', '#eb6834'), ('final method', '#2a78d6')):
        log, pol = runs[name]; t = log[:, 0]; m = (t >= t0) & (t <= t1)
        ax[0].plot(t[m], log[m, 3], color=c, label=name)
        idx = np.array(sorted(i for i in pol.tr if i < len(t))); tt = t[idx]; mm = (tt >= t0) & (tt <= t1)
        kg = np.array([pol.tr[i][1] for i in idx])
        ax[1].plot(tt[mm], kg[mm]/1e3, color=c, label='measurement gain, '+name)
    log = runs['final method'][0]; t = log[:, 0]; m = (t >= t0) & (t <= t1)
    ax[1].plot(t[m], log[m, 8]/1e3, color=INK, lw=1.0, ls='--', label='true stiffness')
    ax[0].set_yscale('symlog', linthresh=10); ax[0].set_ylabel('contact force [N]')
    ax[0].axhline(4, color=MUTED, lw=0.8, ls=':')
    ax[0].legend(frameon=False, loc='upper left')
    ax[1].set_ylabel('stiffness [kN/m]'); ax[1].set_xlabel('time [s]'); ax[1].set_ylim(15, 50)
    ax[1].legend(frameon=False, loc='upper left', fontsize=8)
    for a in ax: a.axvspan(5.05, 5.35, color=GRID, lw=0)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, 'fig1_failure_mechanism.png'), dpi=200); plt.close(fig)


def fig2():
    rows = json.load(open(os.path.join(H, 'results_v2', 'cohort_final.json')))
    d = {}
    for r in rows: d[(r['cfg'], r['sid'])] = r
    rows = list(d.values())
    fig, ax = plt.subplots(1, 5, figsize=(10, 3.0), sharey=True)
    rng = np.random.default_rng(0)
    for k, fam in enumerate(FAMS):
        for j, m in enumerate(METH):
            v = np.array([r['rmse'] for r in rows if r['cfg'] == m and r['family'] == fam and r['ok']])
            x = j+rng.uniform(-0.18, 0.18, len(v))
            ax[k].scatter(x, v, s=8, color=COL[m], marker=MRK[m], alpha=0.55, lw=0)
            ax[k].plot([j-0.3, j+0.3], [np.median(v)]*2, color=INK, lw=1.4)
        ax[k].set_yscale('log'); ax[k].set_title(FLAB[fam], fontsize=9, color=INK)
        ax[k].set_xticks(range(len(METH))); ax[k].set_xticklabels(['Prop.', 'Rossi', 'Wij.', 'DO-MPAC', 'Yuen'], rotation=45, fontsize=8)
    ax[0].set_ylabel('force RMSE [N] (per scene; bar = median)')
    fig.tight_layout(); fig.savefig(os.path.join(OUT, 'fig2_main_comparison.png'), dpi=200); plt.close(fig)


def fig3():
    rows = json.load(open(os.path.join(H, 'results_v2', 'sweep.json')))
    fig, ax = plt.subplots(2, 2, figsize=(8, 5), sharex='col')
    for col, (kind, xl) in enumerate((('ramp', 'push ramp time [s] (faster to the right)'), ('ratio', 'ridge / base stiffness ratio'))):
        doses = sorted({r['dose'] for r in rows if r['kind'] == kind}, reverse=(kind == 'ramp'))
        xs = np.arange(len(doses))
        for m in METH:
            med, bad = [], []
            for dd in doses:
                R = [r for r in rows if r['cfg'] == m and r['kind'] == kind and r['dose'] == dd]
                ok = [r['rmse'] for r in R if r['ok']]
                med.append(np.median(ok)); bad.append(sum((not r['ok']) or r.get('peak', 0) > 50 for r in R))
            ax[0, col].plot(xs, med, color=COL[m], marker=MRK[m], ms=5, label=LAB[m])
            ax[1, col].plot(xs, bad, color=COL[m], marker=MRK[m], ms=5)
        ax[0, col].set_yscale('log'); ax[0, col].set_ylabel('median force RMSE [N]')
        ax[1, col].set_ylabel('scenes failed or >50 N (of 20)'); ax[1, col].set_xlabel(xl)
        ax[1, col].set_xticks(xs); ax[1, col].set_xticklabels([('%g' % dd) for dd in doses])
    h, l = ax[0, 0].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, fontsize=8, loc='lower center', ncol=5)
    fig.tight_layout(rect=(0, 0.06, 1, 1)); fig.savefig(os.path.join(OUT, 'fig3_dose_response.png'), dpi=200); plt.close(fig)


if __name__ == '__main__':
    import sys
    which = sys.argv[1:] or ['1', '2', '3']
    if '1' in which: fig1()
    if '2' in which: fig2()
    if '3' in which: fig3()
    (None if False else None); print('written to', os.path.abspath(OUT))
