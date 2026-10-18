"""Whole-scan capability headroom under a contact/load envelope.

All methods here know the future geometry AND mechanical coefficients. This is
an optimistic planning screen, NOT an online unknown-tissue result. Velocities
only: acceleration, latency, image acquisition and robot dynamics are absent.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
import envelope as e

HERE=Path(__file__).resolve().parent
VMAX, VZMAX = .008, .020


def prepare(p,npath=81,ncontact=41):
    ss=np.linspace(0,.08,npath)
    rows=[]
    for s in ss:
        u,g,k=e.system(p,s,ncontact)
        a,b,off,ub,rhs=e.maps(u,g,k)
        # Winkler pressure is redundant when indentation margin is positive.
        rows.append(dict(a=a,g=g,u=u,k=k,c=b.sum(axis=0),offset=float(off.sum()),
            minimum=e.lp(u,g,k),chord=e.fixed_angle(float(np.clip((g[-1]-g[0])/(2*e.A),-e.ANGLE,e.ANGLE)),u,g,k)))
    return ss,rows


def joint_plan(ss,rows,cap,omega):
    n=len(ss);nv=2*n+n-1
    ii=[];jj=[];vv=[];rhs=[]
    def add(entries,b):
        r=len(rhs)
        for k,v in entries:ii.append(r);jj.append(k);vv.append(v)
        rhs.append(b)
    for i,row in enumerate(rows):
        for a,g in zip(row['a'],row['g']):
            add([(2*i,-a[0]),(2*i+1,-a[1])],-(g+e.MARGIN)/e.MM)
        add([(2*i,row['c'][0]),(2*i+1,row['c'][1])],cap+row['offset'])
    for i in range(n-1):
        for axis,speed in [(0,VZMAX/e.MM),(1,e.A*omega/e.MM)]:
            add([(2*(i+1)+axis,1.),(2*i+axis,-1.),(2*n+i,-speed)],0.)
            add([(2*(i+1)+axis,-1.),(2*i+axis,1.),(2*n+i,-speed)],0.)
    aub=coo_matrix((vv,(ii,jj)),shape=(len(rhs),nv)).tocsr()
    bounds=[]
    for i in range(n):bounds.extend([(-10.,30.),(-e.A*e.ANGLE/e.MM,e.A*e.ANGLE/e.MM)])
    # A common already-established minimum-load starting pose; no free reset.
    bounds[0]=(rows[0]['minimum'][0],rows[0]['minimum'][0])
    bounds[1]=(rows[0]['minimum'][1],rows[0]['minimum'][1])
    bounds.extend([(ds/VMAX,None) for ds in np.diff(ss)])
    cost=np.r_[np.zeros(2*n),np.ones(n-1)]
    sol=linprog(cost,A_ub=aub,b_ub=rhs,bounds=bounds,method='highs')
    if not sol.success:return None,None
    # Lexicographic second objective: among minimum-time paths, reduce spatial
    # force integral; no gain preference or invented target impedance.
    from scipy.sparse import vstack,csr_matrix
    secondary=np.r_[np.concatenate([r['c'] for r in rows]),np.zeros(n-1)]
    rr=linprog(secondary,A_ub=vstack([aub,csr_matrix(cost[None,:])]),
        b_ub=np.r_[rhs,sol.fun+1e-7],bounds=bounds,method='highs')
    if rr.success:sol=rr
    return sol.x[:2*n].reshape(n,2),sol.x[2*n:]


def rate_time(q,ss,seed,omega):
    move=np.abs(np.diff(q,axis=0))*e.MM
    dt=np.maximum.reduce([np.diff(ss)/VMAX,move[:,0]/VZMAX,move[:,1]/(e.A*omega)])
    initial=np.abs(q[0]-seed)*e.MM
    return float(dt.sum()+max(initial[0]/VZMAX,initial[1]/(e.A*omega))),dt


def feasible(q,rows,cap):
    if q is None:return False
    for qq,row in zip(q,rows):
        m=e.metrics(qq,row['u'],row['g'],row['k'])
        if m[0]>cap+1e-7 or m[2]<-1e-8 or m[3]<e.MARGIN-1e-9:return False
    return True


def fixed_angle_best(ss,rows,cap,omega):
    # Find feasible constant-angle interval exactly from all affine constraints.
    lower,upper=-e.ANGLE,e.ANGLE
    for row in rows:
        c0,c1=row['c'];r=row['u']/e.A
        for g,ri in zip((row['g']+e.MARGIN)/e.MM,r):
            aa=(c1-c0*ri)*e.A/e.MM
            bb=cap+row['offset']-c0*g
            if abs(aa)<1e-12:
                if bb<0:return None,None
            elif aa>0:upper=min(upper,bb/aa)
            else:lower=max(lower,bb/aa)
    if lower>upper:return None,None
    # Generous offline strong reference; includes 201 angles in the feasible
    # interval rather than a weak arbitrary fixed orientation.
    best=None
    for theta in np.linspace(lower,upper,201):
        q=np.array([e.fixed_angle(theta,r['u'],r['g'],r['k']) for r in rows])
        tt,dt=rate_time(q,ss,rows[0]['minimum'],omega)
        if best is None or tt<best[0]:best=(tt,q,dt)
    return best[1],best[2]


def dense_check(p,ss,q,cap,ncontact=161):
    min_d,peak=1e9,-1e9
    for i in range(len(ss)-1):
        for fraction in (.25,.5,.75):
            s=ss[i]+fraction*(ss[i+1]-ss[i]);qq=q[i]+fraction*(q[i+1]-q[i])
            u,g,k=e.system(p,s,ncontact)
            m=e.metrics(qq,u,g,k)
            min_d=min(min_d,m[3]);peak=max(peak,m[0])
    return dict(minimum_intersample_penetration=float(min_d),
        peak_intersample_force=float(peak),ceiling_excess=float(max(0.,peak-cap)),
        note='Checks true mechanical touch >=0, not preservation of the 50um design margin or acoustic coupling.')


def main():
    protocol=dict(seed_by_family=[2301,2302],scenes_per_family=12,
        path_points=81,contact_nodes=41,load_design_ceilings=[1.5,2.,3.,6.],
        roll_rate_bounds=[.1,.25,.5],scan_speed_bound=VMAX,normal_speed_bound=VZMAX,
        methods=['joint_pose_progress_oracle','pointwise_min_load_then_rate_scale','chord_then_rate_scale','best_constant_roll'],
        information='All methods know future geometry and the mechanics. Full-contact equilibrium seed is supplied. No claim of online/no-preparation operation.',
        first_phase='min total scan time',second_phase='min spatial force sum at fixed minimum time',
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        envelope_sha256=hashlib.sha256((HERE/'envelope.py').read_bytes()).hexdigest())
    (HERE/'path_protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    rec=[];traces={};dense=[]
    for family,label in enumerate(('smooth_surface','localized_relief')):
        for case,p in enumerate(e.scene(2301+family,12,family)):
            ss,rows=prepare(p)
            point=np.array([r['minimum'] for r in rows]);chord=np.array([r['chord'] for r in rows])
            for cap in protocol['load_design_ceilings']:
                for omega in protocol['roll_rate_bounds']:
                    joint,dt=joint_plan(ss,rows,cap,omega)
                    const,dconst=fixed_angle_best(ss,rows,cap,omega)
                    qs=[joint,point,chord,const]
                    times=[];valid=[]
                    for q in qs:
                        ok=feasible(q,rows,cap);valid.append(ok)
                        times.append(rate_time(q,ss,rows[0]['minimum'],omega)[0] if ok else None)
                    # Joint rate_time at its chosen poses equals minimum time;
                    # do not accidentally use secondary LP's force objective.
                    rec.append(dict(family=label,case=case,cap=cap,omega=omega,
                        time=times,grid_feasible=valid))
                    if case==0 and cap==3. and omega==.25:
                        for m,q in enumerate(qs):
                            if q is not None:traces[label+'_'+str(m)]=q
                        traces[label+'_s']=ss
                        if joint is not None:dense.append(dict(family=label,**dense_check(p,ss,joint,cap)))
            print(label,case,flush=True)
    summary=[]
    for family in ('smooth_surface','localized_relief'):
        for cap in protocol['load_design_ceilings']:
            for omega in protocol['roll_rate_bounds']:
                rows=[r for r in rec if r['family']==family and r['cap']==cap and r['omega']==omega]
                ok=np.array([r['grid_feasible'] for r in rows])
                tt=np.array([[v if v is not None else np.nan for v in r['time']] for r in rows])
                comparison=[]
                for m in range(1,4):
                    common=ok[:,0]&ok[:,m]
                    comparison.append(dict(baseline=protocol['methods'][m],common_cases=int(common.sum()),
                        mean_time_s_joint=float(tt[common,0].mean()) if common.any() else None,
                        mean_time_s_baseline=float(tt[common,m].mean()) if common.any() else None,
                        mean_fractional_time_reduction=float(np.mean(1-tt[common,0]/tt[common,m])) if common.any() else None))
                summary.append(dict(family=family,cap=cap,omega=omega,
                    feasible_count=ok.sum(axis=0).tolist(),comparisons=comparison))
    np.savez_compressed(HERE/'path_traces.npz',**traces)
    (HERE/'path_results.json').write_text(json.dumps(dict(protocol=protocol,summary=summary,
        cases=rec,dense_checks=dense),indent=2)+'\n')


if __name__=='__main__':main()
