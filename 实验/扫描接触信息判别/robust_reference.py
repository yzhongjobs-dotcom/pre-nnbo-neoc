"""Follow-up: strengthen fixed baseline after rare held-out oscillations.

This is explicitly exploratory after inspecting the first held-out set. The new
confirmatory set is drawn once after freezing the selection rule below.
Frozen scalar stability screening does not certify the nonlinear coupled loop.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import screen as s

HERE = Path(__file__).resolve().parent


def radius(bw, beta, scale, tau):
    d, dt = np.exp(-s.DT/tau), s.DT
    # [force, force-velocity, previous force integral, previous force, force n-2]
    integral = np.array([0.,0.,1.,0.,-dt])
    command = scale*(beta*bw*bw*integral-np.array([0.,0.,0.,0.,bw]))
    mat = np.zeros((5,5))
    mat[0] = np.array([1.,tau*(1-d),0.,0.,0.])+(dt-tau*(1-d))*command
    mat[1] = np.array([0.,d,0.,0.,0.])+(1-d)*command
    mat[2],mat[3],mat[4] = integral,np.array([1.,0.,0.,0.,0.]),np.array([0.,0.,0.,1.,0.])
    if beta == 0:
        # Unused integral carries an irrelevant unit eigenvalue.
        mat = mat[np.ix_([0,1,3,4],[0,1,3,4])]
    return float(np.abs(np.linalg.eigvals(mat)).max())


def summary(rr, ref):
    diff = s.score(rr)-s.score(ref)
    rng = np.random.default_rng(9182)
    boot = diff[rng.integers(0,len(diff),size=(10000,len(diff)))].mean(axis=1)
    return dict(total_loss=float(s.score(rr).mean()), median_loss=float(np.median(s.score(rr))),
        p95_loss=float(np.quantile(s.score(rr),.95)),
        scan_force_rms=float(np.sqrt(rr[:,4].mean()/s.SCAN)),
        scan_moment_rms=float(np.sqrt(rr[:,5].mean()/s.SCAN)),
        mean_peak=float(rr[:,2].mean()), worst_peak=float(rr[:,2].max()),
        band_fraction_illustrative=float(rr[:,6].mean()),
        J_error2=float(rr[:,7].mean()),total_time=float(rr[0,8]),
        paired_loss_difference_vs_robust=float(diff.mean()),
        paired_percent_vs_robust=100*float(s.score(rr).mean()/s.score(ref).mean()-1),
        paired_bootstrap_95ci_difference=np.quantile(boot,[.025,.975]).tolist())


def main():
    first = json.loads((HERE/'results.json').read_text())
    # Pointwise stiffness bound implies J <= diag(Kmax, Kmax/3) in Loewner order.
    scale_max = 12e6*np.exp(.65)*(2*s.A*s.L)/1500
    scales = np.geomspace(.04,scale_max,200)
    screened=[]
    for cfg in first['calibration']['fixed_PI']:
        rho=max(radius(cfg['bandwidth'],cfg['integral'],g,tau)
                for g in scales for tau in (.020,.030))
        screened.append(dict(**cfg,max_frozen_scalar_radius=rho))
    admissible=[c for c in screened if c['max_frozen_scalar_radius'] < 1-1e-8]
    robust=min(admissible,key=lambda c:c['mean_total_loss'])
    frozen=dict(original_selected=first['selected'],robust_PI=robust,scales=[float(scales[0]),float(scales[-1])],
        selection="lowest original calibration mean among stable frozen scalar grid; exploratory response to first test tails",
        additional_test_seeds=[9100,9101,9102],additional_noise_seeds=[9200,9201,9202],
        episodes_per_family=128,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (HERE/'followup_protocol.json').write_text(json.dumps(frozen,indent=2)+'\n')
    results, arrays={},{}
    for f,label in enumerate(('uniform_static','heterogeneous_static','heterogeneous_moving')):
        ps=s.environments(9100+f,128,f)
        xs,nn=s.initial_poses(ps),s.noise_for(9200+f,128)
        ref=s.batch(ps,xs,nn,0,robust['bandwidth'],robust['integral'])
        rr=[ref]
        for m,name in enumerate(s.NAMES):
            cfg=first['selected'][name]
            rr.append(s.batch(ps,xs,nn,m,cfg['bandwidth'],cfg['integral']))
        rr=np.array(rr)
        names=['robust_PI']+s.NAMES
        results[label]={name:summary(r,ref) for name,r in zip(names,rr)}
        arrays[label+'_metrics']=rr
        arrays[label+'_parameters']=ps
        arrays[label+'_initial_poses']=xs
        print(label,results[label],flush=True)
    np.savez_compressed(HERE/'followup_paired.npz',**arrays)
    (HERE/'followup.json').write_text(json.dumps(dict(protocol=frozen,screened=screened,test=results),indent=2)+'\n')


if __name__ == '__main__':
    main()
