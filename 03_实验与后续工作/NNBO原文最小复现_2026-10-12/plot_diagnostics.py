"""Python diagnostic figure; simulation trajectories, not measured replicates.

Claim: source-equation interaction runs, but nonzero residual limits optimality.
Panel roles: position tradeoff, applied force, weight trajectory, HJB residual.
Archetype: quantitative grid; no interpolation/downsampling or uncertainty bars.
Scope: tight stated-disturbance bound; alternate bound retained in source CSVs.
"""
from pathlib import Path
import sys
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,'C:/Users/user2/.agents/skills/nature-figure/scripts')
from audit_panel_alignment import require_matplotlib_panel_alignment
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,
                     'pdf.fonttype':42,'svg.fonttype':'none',
                     'axes.spines.right':False,'axes.spines.top':False,
                     'axes.linewidth':.7,'legend.frameon':False})
fig,axes=plt.subplots(2,2,figsize=(7.2,5.6))
fig.subplots_adjust(left=.10,right=.96,bottom=.11,top=.92,wspace=.34,hspace=.48)
colors=['#0072B2','#D55E00']
for i,Q in enumerate([500000,1000000]):
    a=np.load(ROOT/f'Q{Q}_tight_known_disturbance.npz')['data']
    t=a[:,0];color=colors[i];label=f'Q = {Q:,}'
    axes[0,0].plot(t,a[:,1]*1000,color=color,label=label,linewidth=1.2)
    axes[0,1].plot(t,a[:,5],color=color,label=label,linewidth=1.1)
    axes[1,1].semilogy(t,np.abs(a[:,-1]),color=color,label=label,linewidth=1.)
    if i==0:
        axes[0,0].plot(t,a[:,2]*1000,color='#777777',linestyle='--',label='Desired',linewidth=.9)
        axes[0,0].axhline(10,color='#999999',linestyle=':',linewidth=.9,label='Surface')
        for j,c in enumerate(['#0072B2','#D55E00','#009E73']):
            axes[1,0].plot(t,a[:,9+j],color=c,label=f'W{j+1}',linewidth=1.1)
titles=['a  Position tradeoff','b  Applied force (includes exploration)',
        'c  Critic weights (Q = 500,000)','d  Hamiltonian residual magnitude']
ylabels=['Position (mm)','Force (N)','Weight','Absolute residual']
for ax,title,label in zip(axes.flat,titles,ylabels):
    ax.set_title(title,loc='left',fontsize=9,pad=11)
    ax.set_xlabel('Time (s)');ax.set_ylabel(label);ax.set_xlim(0,7)
axes[0,0].legend(loc='upper right',fontsize=7)
axes[0,1].legend(loc='upper right',fontsize=7)
axes[1,0].legend(loc='upper right',fontsize=7)
axes[1,1].legend(loc='upper right',fontsize=7)
fig.canvas.draw()
require_matplotlib_panel_alignment(fig,json_out=str(ROOT/'diagnostics.alignment.json'),
                                  overlay_svg=str(ROOT/'diagnostics.alignment.svg'),
                                  tolerance_pt=1.5,gutter_tolerance_pt=1.5,strict=True)
for fmt in ['pdf','svg','png']:
    fig.savefig(ROOT/f'diagnostics.{fmt}',dpi=220)
plt.close(fig)
(ROOT/'figure_notes.json').write_text(json.dumps({
    'backend':'python','purpose':'diagnostic minimum simulation, not submission figure',
    'data':'full 7001 samples in each Q*_tight_known_disturbance.csv',
    'no_statistical_replicates':True,'no_uncertainty_bars':True,
    'alternate_bound':'uniform_sine_bound retained separately, not hidden or discarded',
    'claim':'finite contact interaction, remaining HJB residual',
    'not_claimed':'exact original figures or HJB optimality'},indent=2),encoding='utf-8')
