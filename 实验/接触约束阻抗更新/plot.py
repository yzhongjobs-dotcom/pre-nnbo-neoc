"""Standalone scientific summary; no generated illustrations."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
r=json.loads((HERE/'ranking_results.json').read_text())
names=['fixed','critical','library','teacher','quadratic','rbf_rank','cost_actor_11','cost_actor_22','cost_actor_33']
rows={v['method']:v for v in r['methods']}
labels=['Fixed','Critical damping','Rollout library','Local optimizer','Quadratic ranker','RBF ranker','Cost actor 11','Cost actor 22','Cost actor 33']
fig,axes=plt.subplots(1,3,figsize=(16,5),layout='constrained',gridspec_kw={'width_ratios':[1.3,1,.9]})
y=np.arange(len(names))
for shift,key,color,label in [(-.18,'mean_improvement_over_critical_pct','#4d89b5','Immediate release'),
                              (.18,'measured_release_improvement_over_critical_pct','#dc8e3d','Measured release delay')]:
    axes[0].barh(y+shift,[rows[n][key] for n in names],height=.34,color=color,label=label)
axes[0].set_yticks(y,labels);axes[0].invert_yaxis();axes[0].axvline(0,color='k',lw=.7)
axes[0].set_xlabel('Mean task improvement over tuned critical damping (%)')
axes[0].legend(loc='lower left',fontsize=8)
axes[1].barh(y,[rows[n]['inference_median_ms'] for n in names],color='#4d89b5')
axes[1].set_yticks(y,labels);axes[1].invert_yaxis();axes[1].set_xscale('log')
axes[1].set_xlabel('Median single-query time (ms, log scale)')
anytime=json.loads((HERE/'anytime_results.json').read_text())['methods']
anames=['library_with_immediate_simple','four_step_local_with_immediate_simple','RBF_warm_four_step_with_immediate_simple']
amap={v['method']:v for v in anytime}
axes[2].barh(np.arange(3),[amap[n]['mean_improvement_over_immediate_simple_pct'] for n in anames],color='#538a59')
axes[2].set_yticks(np.arange(3),['Rollout library','Budgeted local','RBF + budgeted local'])
axes[2].invert_yaxis();axes[2].set_xlabel('Task improvement over simple controller (%)')
axes[2].set_title('Stronger online execution:\nimmediate simple + 20ms refinement',fontsize=9)
fig.suptitle('96 known-model contact contexts; right panel is a post-result baseline improvement on the same cases\nForce thresholds inactive; local teacher is not a global optimum',fontsize=11)
fig.savefig(HERE/'comparison.png',dpi=180)


if __name__=='__main__':plt.close('all')
