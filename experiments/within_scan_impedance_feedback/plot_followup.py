"""Show both later batches, including reversals and model-transfer scope."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
follow = json.loads((HERE / 'followup_analysis.json').read_text())
coupled = json.loads((HERE / 'coupled_audit.json').read_text())['comparisons']
panels = [
    ('New 24: original plant', [
        ('Learned K/D only', follow['fresh24_gain_learning_vs_continuous_constant']),
        ('Learned K/D + FF', follow['fresh24_joint_learning_vs_continuous_constant'])]),
    ('New 24: coupled plant, frozen policies', [
        ('Learned K/D only', coupled['learned_gain_only_vs_continuous_gain_only']),
        ('Learned K/D + FF', coupled['learned_feedback_vs_continuous_gain_only'])]),
]
fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.6), sharex=True)
for ax, (title, rows) in zip(axes, panels):
    for i, (name, value) in enumerate(rows):
        mid = 100 * value['relative_mean_J_benefit']
        lo, hi = [100 * x for x in value['paired_95pct']]
        ax.errorbar(mid, i, xerr=[[mid-lo], [hi-mid]], fmt='o',
                    color=['#3d70a3', '#b77827'][i], capsize=4)
    ax.set_yticks(range(len(rows)), [x[0] for x in rows])
    ax.set_ylim(-.5, 1.5)
    ax.axvline(0, color='0.5', linewidth=1)
    ax.grid(axis='x', alpha=.2)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel('Mean J benefit vs tuned constant (%)\nPaired 95% bootstrap interval; positive = lower J')
fig.suptitle('Follow-up evidence: no extra FF in the gain-only action family', fontsize=11)
fig.tight_layout()
fig.savefig(HERE / 'comparison_followup.png', dpi=180)
