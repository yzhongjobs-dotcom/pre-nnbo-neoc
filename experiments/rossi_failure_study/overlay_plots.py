"""Python plots of retained additive-scene results; no simulation or selection.

Figure contract: compare paired force-control loss, then estimate accuracy on
common observations; fixed i=0 examples illustrate behavior without selecting
favorable episodes. All quantitative summaries include all requested methods.
PNG previews, editable SVG/PDF and measured panel geometry are exported.
"""
import json
import math

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import baseline_compare as B
import overlay_compare as O

OUT=B.HERE/'figures'/'overlay'
TITLES=('Breathing + ribs','Breathing + fast shift','Breathing + ribs + fast shift')
DISPLAY={'Huber-D4':'Huber-D4 (current)','D4':'D4','P5':'P5','Rossi-ID':'Rossi-ID',
         'Kalman':'Kalman','MHE':'MHE','Xue-const':'Xue (original)','Xue-chi':'Xue + gate',
         'Rossi-control':'Rossi full','DO-MPAC':'DO-MPAC','Xing-DAHIC':'Xing DAHIC'}
SELECTED=('Huber-D4','Rossi-ID','MHE','Rossi-control')
COLORS=('#197c80','#93734b','#7a6e9b','#8b949b')


def save(fig,name):
    fig.canvas.draw()
    w,h=fig.get_size_inches()*72
    panels=[]
    for j,ax in enumerate(fig.axes):
        p=ax.get_position()
        spec=ax.get_subplotspec()
        panels.append(dict(id=chr(97+j),bbox_pt=[p.x0*w,p.y0*h,p.x1*w,p.y1*h],
                           grid_id='main',row_start=spec.rowspan.start,row_stop=spec.rowspan.stop,
                           col_start=spec.colspan.start,col_stop=spec.colspan.stop))
    widths=[p['bbox_pt'][2]-p['bbox_pt'][0] for p in panels]
    heights=[p['bbox_pt'][3]-p['bbox_pt'][1] for p in panels]
    assert max(widths)-min(widths)<1.5 and max(heights)-min(heights)<1.5
    manifest=dict(schema_version=1,backend='python-matplotlib',figure=dict(width_pt=w,height_pt=h),panels=panels,exemptions=[])
    (OUT/(name+'.layout.json')).write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    for ext in ('svg','pdf','png'):
        fig.savefig(OUT/(name+'.'+ext),dpi=180)
    plt.close(fig)


def axes():
    fig,aa=plt.subplots(1,3,figsize=(12,4.6))
    fig.subplots_adjust(left=.16,right=.98,bottom=.18,top=.85,wspace=.20)
    for j,(ax,title) in enumerate(zip(aa,TITLES)):
        ax.set_title(chr(97+j)+'  '+title,fontsize=9,pad=14,loc='left')
    return fig,aa


def plot(summary,test,module):
    OUT.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.spines.top':False,
                         'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none',
                         'legend.frameon':False,'axes.labelsize':8,'xtick.labelsize':8,'ytick.labelsize':8})
    # All 10 comparators, pairing already computed by original bootstrap.
    fig,aa=axes();methods=[m for m in test['chosen'] if m!='Huber-D4'];y=np.arange(len(methods))
    for ax,c in zip(aa,O.SCENES):
        for j,m in enumerate(methods):
            d=summary['by_scene'][c][m]['delta'];v=d['mean'];lo,hi=d['ci']
            ax.errorbar(v,j,xerr=[[v-lo],[hi-v]],fmt='o',color='#586975',markersize=4,capsize=2,lw=1)
        ax.axvline(0,color='#197c80',ls='--',lw=1)
        ax.set_yticks(y,labels=[DISPLAY[m] for m in methods] if ax is aa[0] else ['']*len(y))
        ax.set_ylim(len(y)-.5,-.5);ax.set_xlabel('Paired loss difference vs current Huber-D4')
    intervals=[summary['by_scene'][c][m]['delta']['ci'] for c in O.SCENES for m in methods]
    limits=(min(v[0] for v in intervals),max(v[1] for v in intervals));span=limits[1]-limits[0]
    for ax in aa:
        ax.set_xlim(limits[0]-.05*span,limits[1]+.05*span)
    fig.text(.16,.045,'Dots: mean difference in ln(1+J). Bars: 95% cluster bootstrap CI (5000 draws; 9 / 3 / 9 physical clusters).',fontsize=7)
    save(fig,'paired_control')

    fig,aa=axes();methods=list(B.ID_METHODS);y=np.arange(len(methods))
    for ax,c in zip(aa,O.SCENES):
        for j,m in enumerate(methods):
            s=summary['module'][f'{c}:{m}:all']
            ax.plot([s['error_med'],s['error_p90']],[j,j],color='#a8afb4',lw=1.3)
            ax.plot(s['error_med'],j,'o',color='#197c80' if m=='Huber-D4' else '#586975',ms=4)
            ax.plot(s['accepted_med'],j,'|',color='#a65442',ms=8)
        ax.set_yticks(y,labels=[DISPLAY[m] for m in methods] if ax is aa[0] else ['']*len(y))
        ax.set_ylim(len(y)-.5,-.5);ax.set_xscale('log');ax.set_xlabel('Absolute log stiffness error')
    values=[summary['module'][f'{c}:{m}:all'][field] for c in O.SCENES for m in methods
            for field in ('error_med','error_p90','accepted_med')]
    values=[v for v in values if v is not None and v>0]
    for ax in aa:
        ax.set_xlim(min(values)/1.2,max(values)*1.2)
    fig.text(.16,.045,'Dots / line ends: finite-error median / P90. Red ticks: accepted-error median. Both fixed log sources; all eligible samples.',fontsize=7)
    save(fig,'identification')

    for quantity in ('force','stiffness'):
        fig,aa=axes()
        for ax,c in zip(aa,O.SCENES):
            tissue='relaxed' if c=='breath_fast' else 'rib'
            cell=(c,tissue,'neural_white',14.,0)
            sc=B.scene(cell,O.BASE,0.)
            for m,col in zip(SELECTED,COLORS):
                if quantity=='force':
                    r=next(r for r in test['rows'] if r['method']==m and B.key(r)==cell)
                    with np.load(B.OUT/r['log_file']) as z:
                        cols=z['cols'].tolist();t=z['log'][:,cols.index(0)];v=z['log'][:,cols.index(18)]
                    ax.plot(t,v,color=col,lw=1.,label=DISPLAY[m],alpha=.9)
                elif m in B.ID_METHODS:
                    r=next(r for r in module['rows'] if r['method']==m and r['source']=='Huber-D4' and B.key(r)==cell)
                    with np.load(B.OUT/r['samples_file']) as z:
                        ax.plot(z['t'],z['estimate']/1000.,color=col,lw=1.,label=DISPLAY[m],alpha=.9)
                        truth=z['truth'].copy();t=z['t'].copy()
            if quantity=='force':
                ax.axhline(8.,color='black',ls='--',lw=.8)
            else:
                ax.plot(t,truth/1000.,color='black',ls='--',lw=1,label='True probe gain')
                ax.set_yscale('log')
            for _,start,end in sc['events']:
                ax.axvspan(start,end,color='#ccd4d9',alpha=.35,zorder=-2)
            ax.set_xlabel('Time (s)');ax.set_xlim(0,10)
            ax.set_ylabel('Axial force (N)' if quantity=='force' else 'Probe gain (kN/m)')
        handles,labels=aa[-1].get_legend_handles_labels()
        fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.56,.98),ncol=len(handles),fontsize=8)
        fig.text(.16,.045,'Fixed examples: neural_white, 14 degrees, seed index 0. Gray shading: fast shift. '+
                 ('Dashed line: 8 N target.' if quantity=='force' else 'All estimates use the same Huber-D4 observations; no full Rossi stiffness shown.'),fontsize=7)
        save(fig,'example_'+quantity)


if __name__=='__main__':
    plot(O.ORIGINAL_LOAD('overlay_summary.json'),O.ORIGINAL_LOAD('overlay_test.json'),O.ORIGINAL_LOAD('overlay_module.json'))
