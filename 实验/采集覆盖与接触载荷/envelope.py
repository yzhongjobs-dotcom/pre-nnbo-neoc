"""Quasistatic whole-path contact-coverage/load screen, not an US simulator.

Winkler derivation and independent shear-coupled elastic-foundation audit.
Known geometry is shared with the strong geometric/search baselines.
No claim of zero preparation, clinical force limit or dynamic safety.
"""
from pathlib import Path
import hashlib
import json
import platform
import time
import numpy as np
from scipy.linalg import solve
from scipy.optimize import linprog, minimize_scalar

HERE = Path(__file__).resolve().parent
A, L, R, MM = .015, .010, .050, .001
MARGIN, LIMIT, ANGLE = .00005, 6., .25
METHODS = ['full_profile_LP', 'two_read_reduction', 'geometry_chord',
           'direct_scalar_search', 'best_zero_moment', 'noisy_two_read']


def scene(seed, n, family):
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        out.append([np.exp(rng.uniform(np.log(.6e6),np.log(6e6))),
                    np.exp(rng.uniform(np.log(.6e6),np.log(6e6))),
                    rng.uniform(.025,.055),rng.uniform(.003,.010),
                    rng.uniform(-.05,.05),rng.uniform(-.001,.001),
                    rng.uniform(.0002,.0005),rng.uniform(-np.pi,np.pi),
                    rng.uniform(.001,.002) if family else 0.,
                    rng.uniform(.004,.010),rng.uniform(.025,.055)])
    return np.array(out)


def material(xx,p):
    blend = .5*(1+np.tanh((xx-p[2])/p[3]))
    return np.exp((1-blend)*np.log(p[0])+blend*np.log(p[1]))


def geometry(xx,p):
    return (p[5]+p[4]*(xx-.04)+p[6]*np.cos(2*np.pi*xx/.08+p[7])
            +p[8]*np.exp(-.5*((xx-p[10])/p[9])**2))


def system(p,s,n=61,shear_length=0.):
    # Full physical elastic domain is three times the probe width. Condense
    # the unconstrained tissue outside the probe (not a clamped edge at probe).
    du = 2*A/(n-1)
    uu = np.arange(-(3*(n-1)//2),3*(n-1)//2+1)*du
    weights = np.full(len(uu),L*du)
    weights[[0,-1]] *= .5
    ground = weights*material(s+uu,p)
    k = np.diag(ground)
    shear = 3e6*L*shear_length**2/du
    for i in range(len(uu)-1):
        k[i,i] += shear;k[i+1,i+1] += shear
        k[i,i+1] -= shear;k[i+1,i] -= shear
    idx = np.flatnonzero(np.abs(uu) <= A+1e-12)
    outside = np.setdiff1d(np.arange(len(uu)),idx)
    eff = k[np.ix_(idx,idx)]
    if shear_length:
        eff = eff-k[np.ix_(idx,outside)]@solve(k[np.ix_(outside,outside)],
             k[np.ix_(outside,idx)],assume_a='pos')
    else:
        # Trapezoidal footprint weights, including half cells at probe ends.
        eff = eff.copy();eff[0,0]*=.5;eff[-1,-1]*=.5
    u = uu[idx]
    g = geometry(s+u,p)+u*u/(2*R)
    return u,g,eff


def maps(u,g,k):
    a = np.column_stack((np.ones(len(u)),u/A))
    b = MM*k@a
    off = k@g
    # Non-negative normal forces AND nominal geometric indentation margin.
    ub = np.vstack((-a,-b))
    rhs = np.r_[-(g+MARGIN)/MM,-off]
    return a,b,off,ub,rhs


def lp(u,g,k,c=None,zero_moment=False,pressure=True):
    a,b,off,ub,rhs=maps(u,g,k)
    coeff=b.sum(axis=0) if c is None else np.array([1.,c])
    ae,be=None,None
    if zero_moment:
        ae=(u@b)[None,:];be=np.array([u@off])
    if not pressure:ub,rhs=-a,-(g+MARGIN)/MM
    sol=linprog(coeff,A_ub=ub,b_ub=rhs,A_eq=ae,b_eq=be,
        bounds=[(-10.,30.),(-A*ANGLE/MM,A*ANGLE/MM)],method='highs')
    if not sol.success:return None
    return sol.x


def fixed_angle(theta,u,g,k):
    _,b,off,_,_=maps(u,g,k)
    rr=u/A
    roll=A*theta/MM
    z=max(np.max((g+MARGIN)/MM-rr*roll),np.max((off-b[:,1]*roll)/b[:,0]))
    return np.array([z,roll])


def metrics(q,u,g,k):
    if q is None:return [np.nan,np.nan,np.nan,np.nan]
    d=MM*(q[0]+q[1]*u/A)-g
    f=k@d
    # f>=0 characterizes whether the full-contact candidate is an equilibrium.
    # Negative entries invalidate this assumed active set, not negative force.
    return [float(f.sum()),float(u@f),float(f.min()),float(d.min())]


def evaluate(p,s,rng,n=61,shear_length=0.,noise_sd=(.03,.0003),delta=.0003,averages=16):
    u,g,k=system(p,s,n,shear_length)
    a,b,off,_,_=maps(u,g,k)
    q0=fixed_angle(0.,u,g,k)
    # A common already-established full-contact seed; its cost is NOT zero.
    q0[0] += .1
    m0=metrics(q0,u,g,k)
    q1=q0+np.array([delta/MM,0.])
    m1=metrics(q1,u,g,k)
    c_exact=(m1[1]-m0[1])/(A*(m1[0]-m0[0]))
    full=lp(u,g,k)
    # Reduction uses geometry + differential wrench only. With shear coupling,
    # using untransformed geometry is deliberately a model-mismatch test.
    reduced=lp(u,g,k,c_exact,pressure=False)
    theta=float(np.clip((g[-1]-g[0])/(2*A),-ANGLE,ANGLE))
    chord=fixed_angle(theta,u,g,k)
    count=0
    def fun(th):
        nonlocal count
        count+=1
        return metrics(fixed_angle(th,u,g,k),u,g,k)[0]
    opt=minimize_scalar(fun,bounds=(-ANGLE,ANGLE),method='bounded',
                        options={'xatol':1e-9,'maxiter':100})
    search=fixed_angle(opt.x,u,g,k)
    zero=lp(u,g,k,zero_moment=True)
    noise=rng.normal(size=(2,2))*np.array(noise_sd)/np.sqrt(averages)
    dy=np.array([m1[0]-m0[0],(m1[1]-m0[1])/A])+np.array([noise[1,0]-noise[0,0],(noise[1,1]-noise[0,1])/A])
    # Synthetic SNR safeguard. Unsafe calibration is not performed; fallback
    # chooses c=0. This does not repair an unsafe seed or prove safe engagement.
    performed=m1[0] <= LIMIT and m0[0] <= LIMIT
    ce=np.clip(dy[1]/dy[0],-.99,.99) if performed and dy[0]>.03 else 0.
    noisy=lp(u,g,k,float(ce),pressure=False)
    qs=[full,reduced,chord,search,zero,noisy]
    mm=np.array([metrics(q,u,g,k) for q in qs])
    return mm,dict(search_calls=count,calibration_seed_force=m0[0],
        calibration_step_force=m1[0],calibration_performed=performed,c_exact=float(c_exact),
        c_noisy=float(ce),noise_sd=list(noise_sd),delta=delta,averages=averages,
        geometry_convex=bool(np.min(np.diff(g,n=2))>=-1e-12)),np.array([q if q is not None else [np.nan,np.nan] for q in qs])


def aggregate(raw):
    # One row = one scene/position, methods paired. Whole-scan aggregation is
    # separate, so adjacent positions are NOT treated as independent subjects.
    out={}
    for m,name in enumerate(METHODS):
        ff=raw[:,m,0]
        feasible=np.isfinite(ff)&(raw[:,m,2]>=-1e-8)&(raw[:,m,3]>=MARGIN-1e-9)
        out[name]=dict(full_contact_feasible_fraction=float(feasible.mean()),
            within_6N_and_contact_fraction=float((feasible&(ff<=LIMIT+1e-8)).mean()),
            mean_force_on_feasible=float(np.mean(ff[feasible])) if feasible.any() else None,
            mean_excess_force_on_common_feasible=float(np.mean((ff-raw[:,0,0])[feasible])) if feasible.any() else None,
            max_excess_force_on_common_feasible=float(np.max((ff-raw[:,0,0])[feasible])) if feasible.any() else None,
            max_abs_force_difference_to_full=float(np.nanmax(np.abs(ff-raw[:,0,0]))))
    return out


def main():
    start=time.perf_counter();rng=np.random.default_rng(5120)
    protocol=dict(methods=METHODS,model='elastic foundation, known geometry, planar quasistatic pose; no ultrasound',
        a=A,length=L,radius=R,geometric_indentation_margin=MARGIN,
        illustrative_load_ceiling=LIMIT,angle_bound=ANGLE,force_noise_sd=.03,moment_noise_sd=.0003,
        calibration_displacement=.0003,averages_per_read=16,scene_seeds=[1701,1702],
        scenes_per_family=48,path_positions=np.linspace(0,.08,17).tolist(),
        coupled_scene_seed=1901,coupled_scenes=12,coupled_path_positions=np.linspace(0,.08,9).tolist(),
        coupled_lengths=[.0005,.0015,.003],
        note='Contact coverage is mechanical support, not diagnostic image quality. Two reads require a common full-contact seed and same location. Search calls assume an ideal complete-contact depth oracle, so are optimistic.',
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (HERE/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    results,arrays={},{}
    for family,label in enumerate(('smooth_surface','localized_relief')):
        ps=scene(1701+family,48,family)
        raw,meta,poses=[],[],[]
        for p in ps:
            for s in protocol['path_positions']:
                mm,md,qq=evaluate(p,s,rng);raw.append(mm);meta.append(md);poses.append(qq)
        raw=np.array(raw);arrays[label+'_metrics']=raw;arrays[label+'_poses']=np.array(poses)
        arrays[label+'_parameters']=ps
        # Group into 48 whole paths, report path-wise all-position feasibility.
        shaped=raw.reshape(48,17,len(METHODS),4)
        complete=(np.isfinite(shaped[:,:,:,0])&(shaped[:,:,:,0]<=LIMIT+1e-8)
            &(shaped[:,:,:,2]>=-1e-8)&(shaped[:,:,:,3]>=MARGIN-1e-9)).all(axis=1)
        results[label]=dict(points=aggregate(raw),whole_path_feasible_fraction={name:float(complete[:,m].mean()) for m,name in enumerate(METHODS)},
            mean_search_calls=float(np.mean([x['search_calls'] for x in meta])),
            calibration_performed_fraction=float(np.mean([x['calibration_performed'] for x in meta])),
            calibration_seed_over_6N_fraction=float(np.mean([x['calibration_seed_force']>LIMIT for x in meta])),
            geometry_convex_fraction=float(np.mean([x['geometry_convex'] for x in meta])),
            max_seed_force=float(max(x['calibration_seed_force'] for x in meta)))
        arrays[label+'_meta']=np.array([[md['c_exact'],md['c_noisy'],md['calibration_seed_force'],md['calibration_step_force'],md['geometry_convex']] for md in meta])
        print(label,results[label],flush=True)
    ps=scene(1901,12,1)
    for ell in protocol['coupled_lengths']:
        raw=[]
        for p in ps:
            for s in protocol['coupled_path_positions']:
                mm,_,_=evaluate(p,s,rng,shear_length=ell);raw.append(mm)
        raw=np.array(raw);label='coupled_'+str(ell)
        arrays[label+'_metrics']=raw
        results[label]=aggregate(raw)
        print(label,results[label],flush=True)
    arrays['coupled_parameters']=ps
    # Explicit deterministic zero-moment/full-contact incompatibility example.
    alpha=.8 # right/left material ratio 9, larger than the threshold 5.
    results['zero_moment_example']=dict(ratio=9,alpha=alpha,flat_face_full_support_min_indentation_factor=1-1.5*alpha,
        scope='Flat probe, uniform planar surface, two constant halves, zero torque; full-support equations imply a negative edge indentation. Not all force-based controllers impose zero torque.')
    np.savez_compressed(HERE/'paired.npz',**arrays)
    (HERE/'results.json').write_text(json.dumps(dict(protocol=protocol,results=results,
        runtime_s=time.perf_counter()-start,python=platform.python_version(),numpy=np.__version__),indent=2)+'\n')


if __name__=='__main__':main()
