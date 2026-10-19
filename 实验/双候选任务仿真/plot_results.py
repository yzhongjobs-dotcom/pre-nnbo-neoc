"""Standalone scientific figures; paired cases and illustrative trajectories."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path(__file__).resolve().parent
def rows(name):return json.loads((root/(name+'.json')).read_text())['results']
c1=rows('candidate1_test');nets=[rows(f'candidate2_test_seed{s}') for s in [16111,16112,16113]];base=rows('jin_ms_test')
fig,ax=plt.subplots(2,3,figsize=(12,6.5),constrained_layout=True)
colors=['#2b6cb0','#b86520','#24764e'];labels=['C1: value + neighboring update','C2: guided NN','Jin MS task adaptation']
for k,key,title in [(0,'force_rmse_N','Force RMSE (N)'),(1,'tangent_rmse_mm','Tangential RMSE (mm)'),(2,'normal_rmse_mm','Normal RMSE (mm)')]:
    values=np.array([[r[key] for r in c1],np.mean([[r[key] for r in rr] for rr in nets],axis=0),[r[key] for r in base]])
    for j in range(24):ax[0,k].plot(range(3),values[:,j],color='.8',lw=.5,zorder=0)
    for i,c in enumerate(colors):ax[0,k].scatter(np.full(24,i),values[i],s=10,color=c,alpha=.5);ax[0,k].scatter(i,values[i].mean(),s=65,color=c,marker='D',edgecolor='black')
    ax[0,k].set_xticks(range(3),['C1','C2 (3-seed mean)','Jin MS']);ax[0,k].set_ylabel(title);ax[0,k].grid(axis='y',alpha=.2)
for i,name in enumerate(['candidate1_test','candidate2_test_seed16111','jin_ms_test']):
    with np.load(root/(name+'_traces.npz')) as data:a=data['3']
    for k,column,scale in [(0,10,1),(1,1,1000),(2,6,1)]:ax[1,k].plot(a[:,0],a[:,column]*scale,color=colors[i],label=labels[i],lw=1.2)
ax[1,0].axhline(.6,color='black',ls='--',lw=.7,label='Desired force')
for a,title in zip(ax[1],['Force, combined-change case (N)','Tangential error, same case (mm)','Actual normal stiffness (N/m)']):a.set_xlabel('Time (s)');a.set_ylabel(title);a.grid(alpha=.2)
ax[1,0].legend(fontsize=7,loc='best');fig.suptitle('Synthetic continuous scanning: 24 paired cases; lower errors are better',fontsize=12)
fig.savefig(root/'results.png',dpi=220);fig.savefig(root/'results.pdf');plt.close(fig)
