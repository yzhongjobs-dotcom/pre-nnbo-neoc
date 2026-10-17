"""Exportable plots; seed variation is shown, not hidden by best-seed selection."""
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from contact_env import ROOT

out = ROOT/'results'
test = json.loads((out/'test.json').read_text())
fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
for h, color in [(4, '#1764ab'), (1, '#c65917')]:
    for i, seed in enumerate([11, 22, 33]):
        rows = json.loads((out/f'learning_h{h}_s{seed}.json').read_text())
        axes[0].plot([r['steps']/1000 for r in rows], [r['mean_cost'] for r in rows],
                     color=color, alpha=.65, marker='.', label=f'PPO history {h}' if i==0 else None)
axes[0].set(xlabel='Training decisions (thousands)', ylabel='Validation physical cost J')
axes[0].legend()
labels = ['ID', 'Material OOD', 'Noise + mass OOD']
for offset, name, label, color in [(-.27, 'fixed', 'Fixed target', '#777777'),
         (-.09, 'feedback', 'Simple feedback', '#d89732'),
         (.09, 'preview', 'Causal model preview', '#47a291'),
         (.27, 'prior_preview', 'Known-program schedule', '#9567b2')]:
    vals = [test[g][name]['mean_cost'] for g in test]
    axes[1].bar(np.arange(3)+offset, vals, width=.15, color=color, label=label)
for h, marker, color in [(4, 'o', '#1764ab'), (1, 'x', '#c65917')]:
    for i, seed in enumerate([11, 22, 33]):
        vals = [test[g][f'ppo_h{h}_s{seed}']['mean_cost'] for g in test]
        axes[1].scatter(np.arange(3)+(i-1)*.035, vals, marker=marker, color=color,
                        s=28, label=f'PPO history {h}' if i==0 else None, zorder=3)
axes[1].set(xticks=np.arange(3), xticklabels=labels, ylabel='Test physical cost J')
axes[1].legend(fontsize=7)
fig.savefig(out/'quick_screen.png', dpi=180)
fig.savefig(out/'quick_screen.svg')
