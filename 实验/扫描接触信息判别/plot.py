from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p=Path(__file__).resolve().parent
r=json.loads((p/'followup.json').read_text())['test']
methods=['robust_PI','fixed_PI','passive_fit','preprobe_fit','embedded_probe_fit','local_J_oracle']
labels=['1 Range-screened PI','2 Mean-tuned PI','3 Passive fit','4 Pre-scan probe','5 Embedded probe','6 Local-J oracle']
colors=plt.cm.tab10(np.arange(6))
fig,ax=plt.subplots(2,3,figsize=(13,7),layout='constrained')
for c,(family,rows) in enumerate(r.items()):
    for i,(name,label,color) in enumerate(zip(methods,labels,colors)):
        v=rows[name]
        ax[0,c].bar(i,v['scan_force_rms'],color=color)
        ax[1,c].scatter(v['J_error2'],v['total_loss'],color=color,s=65,label=label)
        offset=(-12,5) if i==3 else ((5,-12) if i==4 else (5,4))
        ax[1,c].annotate(str(i+1),(v['J_error2'],v['total_loss']),xytext=offset,textcoords='offset points',fontsize=9)
    ax[0,c].set(title=family.replace('_',' '),ylabel='Scan normal-force RMS [N]',xticks=np.arange(6),xticklabels=np.arange(1,7))
    ax[0,c].grid(axis='y',alpha=.2)
    ax[1,c].set(xscale='log',yscale='log',xlabel='Mean squared local-J error',ylabel='Total normalized-wrench loss [N² s]')
    ax[1,c].grid(alpha=.2)
    ax[1,c].margins(x=.12,y=.12)
fig.legend(labels,loc='outside lower center',ncol=3,frameon=False)
fig.suptitle('Independent synthetic screen: 128 paired environments per family\nLower model error need not produce lower task loss; oracle is unavailable knowledge')
fig.savefig(p/'screen.png',dpi=180)
fig.savefig(p/'screen.svg')
