"""Primary independent-test comparisons, task components and family breakdown."""
from pathlib import Path
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent


def main():
    results=json.loads((HERE/'results.json').read_text());data=np.load(HERE/'paired.npz')
    keys=['fixed_PI','continuous_constant','analytic_inverse','equilibrium_reference','spatial_joint','learned_feedback']
    labels=['Fixed PD + force P','Tuned constant + FF','Analytic inverse K','Equilibrium reference','Predictive gain + FF','Offline feedback']
    out={};rng=np.random.default_rng(8491)
    for family in range(3):
        rows={};mask=np.arange(48)%3==family
        for name in keys:
            v=data['informed_prior_'+name][mask]
            base=data['informed_prior_analytic_inverse'][mask]
            ids=rng.integers(0,len(v),size=(4000,len(v)))
            relative=1-v[ids,0].mean(1)/base[ids,0].mean(1)
            rows[name]={'mean_J':float(v[:,0].mean()),'force_rms_N':float(v[:,1].mean()),
              'normal_pose_rms_mm':float(v[:,10].mean()*1000),'worst_peak_force_N':float(v[:,3].max()),
              'relative_J_benefit_vs_analytic':float(1-v[:,0].mean()/base[:,0].mean()),
              'paired95pct':np.quantile(relative,[.025,.975]).tolist()}
        out[str(family)]=rows
    primary=results['joint']['informed_prior']
    # A shared objective is assessed on both force and normal motion, not force alone.
    fig,axes=plt.subplots(1,3,figsize=(14,4.1),layout='constrained')
    colors=['.6','.65','#86a9a5','#9a9bbf','#e8aa60','#3f7eae']
    ax=axes[0];benefits=[];lo=[];hi=[]
    for name in keys:
        r=primary[name];val=r['benefit_vs_fixed_ratio_of_mean_J']*100
        ci=np.array(r['paired_benefit_vs_fixed_95pct'])*100
        benefits.append(val);lo.append(max(0,val-ci[0]));hi.append(max(0,ci[1]-val))
    ax.barh(np.arange(len(keys)),benefits,color=colors)
    ax.errorbar(benefits,np.arange(len(keys)),xerr=[lo,hi],fmt='none',ecolor='black',capsize=2)
    ax.set_yticks(np.arange(len(keys)),labels);ax.axvline(0,color='black',lw=.6)
    ax.set_xlabel('Joint-cost improvement vs fixed (%)');ax.set_title('48 independent scans; paired 95% intervals')
    ax.invert_yaxis()
    ax=axes[1]
    for name,label,color in zip(keys,labels,colors):
        r=primary[name];ax.scatter(r['normal_pose_rms_mm'],r['force_rms_N'],s=40,label=label,color=color)
    ax.set_xlabel('Mean normal-position RMS (mm)');ax.set_ylabel('Mean force RMS (N)')
    ax.set_title('Both task components; target force 3 N');ax.grid(alpha=.2)
    ax=axes[2]
    for name,label,color in zip(keys,labels,colors):
        t=results['timing'][name];r=primary[name]
        ax.scatter(t['median_ms'],r['mean_J'],s=40,label=label,color=color)
        ax.plot([t['median_ms'],t['p95_ms']],[r['mean_J']]*2,lw=1,color=color)
    ax.set_xscale('log');ax.set_xlabel('Response + decision time (ms)');ax.set_ylabel('Mean joint cost')
    ax.set_title('Prototype timing; full I/O excluded');ax.grid(alpha=.2)
    fig.savefig(HERE/'comparison.png',dpi=180);plt.close(fig)
    payload={'families':out,'scope':'Per-environment bootstrap, conditional on one fitted policy; no physical generalization claim',
      'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      'source_results_sha256':hashlib.sha256((HERE/'results.json').read_bytes()).hexdigest()}
    (HERE/'analysis.json').write_text(json.dumps(payload,indent=2)+'\n')
    for family,rows in out.items():
        print('family',family,'learning vs analytic',rows['learned_feedback']['relative_J_benefit_vs_analytic'],rows['learned_feedback']['paired95pct'])
    print('all methods maxima',{name:primary[name]['worst_peak_force_N'] for name in keys})
if __name__=='__main__':main()
