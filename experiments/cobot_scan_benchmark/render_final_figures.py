"""Render the two figures of the final simulation section from final_figure_traces.npz."""
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 9, 'font.family': 'DejaVu Sans'})
D = np.load('completed_tasks/final_figure_traces.npz')
M = [('Fixed_stiff', 'Fixed (stiff-matched)', '#888888'), ('EtR-A-Haninger', 'Haninger', '#8172B2'),
     ('EtR-A', 'Online estimate + rule', '#4C72B0'), ('Candidate2', 'C2', '#55A868'), ('Candidate1_v2', 'C1', '#C44E52')]

def shade(ax, p):
    on, ramp, hold = p[16], p[18], p[19]
    for a in ax:
        a.axvspan(on, on + ramp, color='orange', alpha=.15); a.axvspan(on + ramp + hold, on + 2 * ramp + hold, color='orange', alpha=.15)

if __name__ == '__main__':
    # Fig. 1: T3 scene 0
    fig, ax = plt.subplots(3, 1, figsize=(7, 6), sharex=True, gridspec_kw={'height_ratios': [1, 2, 1.2]})
    for k, lab, c in M:
        L = D[f'T3_combined_0_{k}']; t = L[:, 0]
        if k == 'Fixed_stiff': ax[0].plot(t, L[:, 8] / 1e3, 'k', lw=1)
        ax[1].plot(t, L[:, 3], color=c, lw=1.1 if k == 'Candidate1_v2' else .8, label=lab)
        ax[2].plot(t, L[:, 9], color=c, lw=1.1 if k == 'Candidate1_v2' else .8)
    shade(ax, D['T3_combined_0_p'])
    ax[1].axhspan(3.5, 4.5, color='.9', zorder=0); ax[1].axhline(7, ls='--', color='.5', lw=.8)
    ax[0].set_ylabel('true local\nstiffness (kN/m)'); ax[1].set_ylabel('contact force (N)'); ax[2].set_ylabel('damping D (N s/m)')
    ax[2].set_xlabel('time (s)'); ax[1].legend(ncol=3, loc='upper left', fontsize=7.5, frameon=False); ax[2].set_xlim(1, 10)
    fig.tight_layout(); fig.savefig('completed_tasks/figures/final_fig_T3.png', dpi=200)
    # Fig. 2: T2 soft (0) and stiff (2) push
    fig, ax = plt.subplots(2, 2, figsize=(8, 4.8), sharex=True)
    for j, (i, name) in enumerate(((0, 'soft tissue, 5 mm push'), (2, 'stiff tissue, 5 mm push'))):
        for k, lab, c in M:
            L = D[f'T2_motion_{i}_{k}']; t = L[:, 0]
            ax[0, j].plot(t, L[:, 3], color=c, lw=1.1 if k == 'Candidate1_v2' else .8, label=lab)
            ax[1, j].plot(t, L[:, 2] * 1e3, color=c, lw=1.1 if k == 'Candidate1_v2' else .8)
        shade(ax[:, j], D[f'T2_motion_{i}_p']); ax[0, j].axhspan(3.5, 4.5, color='.9', zorder=0)
        ax[0, j].set_title(name); ax[1, j].set_xlabel('time (s)'); ax[0, j].set_xlim(1, 8)
    ax[0, 0].set_ylabel('contact force (N)'); ax[1, 0].set_ylabel('probe normal position (mm)')
    ax[0, 0].legend(fontsize=7, frameon=False)
    fig.tight_layout(); fig.savefig('completed_tasks/figures/final_fig_T2.png', dpi=200)
    print('ok')
