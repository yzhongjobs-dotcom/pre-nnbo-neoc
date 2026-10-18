"""Export a compact comparison figure from retained primary results."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent


def main():
    read=lambda name:json.loads((HERE/name).read_text())
    steady=read('results.json')['results'];updates=read('update_results.json')['updates']
    transient=read('transient_results.json')['cases'];adjoint=read('adjoint_results.json')['cases']
    bounded=read('bounded_results.json')['cases'];delay=read('delay_results.json')['cases']
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    labels=['White force','Colored force','Colored + observation noise']
    families=list(steady);colors=['#0072B2','#D55E00','#009E73']
    ax=axes[0,0];xx=np.arange(4)
    for i,(family,label) in enumerate(zip(families,labels)):
        yy=[steady[family]['methods'][m]['mean_relative_gap']*100 for m in ['fixed_tuned','analytic_schedule','teacher_RBF','actor_11']]
        ax.bar(xx+(i-1)*.25,np.maximum(yy,1e-7),width=.25,label=label,color=colors[i])
    ax.set_yscale('log');ax.set_xticks(xx,['Tuned fixed','Analytic','RBF teacher','Actor seed 11'])
    ax.set_ylabel('Mean cost gap to numerical reference (%)');ax.set_title('Stationary four-gain screen, 96 test cases/family')
    ax.legend(fontsize=7);ax.grid(axis='y',alpha=.2)
    ax=axes[0,1];methods=['warm_Newton','warm_chord','sensitivity_chord']
    for i,(family,label) in enumerate(zip(families,labels)):
        yy=[np.mean([r['time_s'] for r in updates if r['family']==family and r['method']==m])*1000 for m in methods]
        ax.bar(xx[:3]+(i-1)*.25,yy,width=.25,color=colors[i])
    ax.set_xticks(xx[:3],['Warm Newton','Warm chord','Sensitivity + chord']);ax.set_ylabel('Mean update time (ms), including preparation')
    ax.set_title('Same projected-gradient tolerance, 1152 updates/method');ax.grid(axis='y',alpha=.2)
    ax=axes[1,0]
    for method,name,color in [(0,'Hold old','#999999'),(1,'Slew to stationary','#D55E00'),(4,'Two blocks, converged','#0072B2')]:
        yy=[c['methods'][method]['gap_to_reference']*100 for c in transient]
        ax.plot(np.arange(8),yy,'o-',label=name,color=color,ms=4,lw=1)
    yy=[r['history'][-1]['gap_to_12_block_reference']*100 for r in bounded]
    ax.plot(np.arange(8),yy,'s-',label='Two blocks, <=4 steps',color='#009E73',ms=4,lw=1)
    ax.set_xticks(np.arange(8),['W0','W1','W2','W3','N0','N1','N2','N3'])
    ax.set_ylabel('Cost gap to best available 12-block trajectory (%)')
    ax.set_title('0.3 s transients, actual gain-rate constraints');ax.legend(fontsize=7);ax.grid(alpha=.2)
    ax=axes[1,1]
    names=['instantaneous','5ms','20ms','50ms','100ms'];delays=[0,5,20,50,100]
    for i,family in enumerate(['white','colored_observation_noise']):
        rr=[r for r in delay if r['family']==family]
        yy=[100*np.mean([next(v for v in r['releases'] if v['release']==name)['improvement_over_hold'] for r in rr]) for name in names]
        ax.plot(delays,yy,'o-',label='White' if i==0 else 'Colored + observation noise',color=colors[2*i])
    ax.axhline(0,color='gray',lw=.6);ax.set_xlabel('Plan release delay (ms)')
    ax.set_ylabel('Mean cost improvement over holding old gains (%)')
    ax.set_title('Old gains hold during compute; then execute returned plan')
    ax.legend(fontsize=7);ax.grid(alpha=.2)
    fig.suptitle('Optimal impedance update: control value and computation\nLocal synthetic models; no full B0 or robot deployment claim',fontsize=12)
    fig.tight_layout(rect=(0,0,1,.93));fig.savefig(HERE/'comparison.png',dpi=170);plt.close(fig)


if __name__=='__main__':main()
