"""Independent solver, integration, replay and offline-foresight diagnostics."""
import json,time
import numpy as np
from scipy.optimize import minimize,brentq
import scan_validation as v

HERE=v.HERE
def write(name,obj):
    (HERE/name).write_text(json.dumps(obj,indent=2),encoding='utf8')

def replay(kind,speed,gains,targets,dt):
    duration=v.CFG['scan_length']/speed;steps=round(duration/dt)
    y=v.initial(speed,gains)[None,:];cost=0.;rows=[]
    for j in range(steps):
        t=j*dt;index=min(len(targets)-1,int((t+1e-10)/duration*len(targets)))
        goal=targets[index][None,:]
        _,_,_,out=v.dynamics(y,t,speed,goal,None,kind)
        fn,e,delta,command,u,vt=out
        rows.append([t,*y[0],fn[0],*e[0],delta[0]])
        y,l,p=v.rk4(y,t,dt,speed,goal,None,kind);cost+=l[0]
    a=np.array(rows)
    return dict(cost=float(cost),force_rmse_N=float(np.sqrt(np.mean((a[:,10]-v.CFG['force_target'])**2))),
                track_rmse_m=float(np.sqrt(np.mean(a[:,11]**2))),normal_rmse_m=float(np.sqrt(np.mean(a[:,12]**2))),
                force_peak_N=float(a[:,10].max()),min_energy_J=float(a[:,9].min()),
                min_contact_m=float(a[:,13].min())),a

def offline_value(q,state,speed,kind,dt=.008,blocks=4):
    # Full true spatial field is intentionally supplied only to this OFFLINE diagnostic.
    q=np.atleast_2d(q);y=np.repeat(state[None,:],len(q),axis=0)
    targets=v.LOW+q.reshape(len(q),blocks,4)*v.SPAN
    duration=v.CFG['scan_length']/speed;steps=round(duration/dt);dt=duration/steps
    cost=np.zeros(len(q))
    for j in range(steps):
        goal=targets[:,min(blocks-1,j*blocks//steps)]
        y,l,p=v.rk4(y,j*dt,dt,speed,goal,None,kind);cost+=l+p
    return cost

def oracle(kind,speed,gains):
    state=v.initial(speed,gains);q=np.tile((gains-v.LOW)/v.SPAN,4)
    def obj(q):
        ep=2e-5;eye=np.eye(len(q))*ep
        costs=offline_value(np.vstack((q,q+eye,q-eye)),state,speed,kind)
        return costs[0],(costs[1:1+len(q)]-costs[1+len(q):])/(2*ep)
    starts=[q,q.copy()]
    starts[1][1::4]=.1;starts[1][3::4]=.2
    out=[];tic=time.perf_counter()
    incumbent=float(offline_value(q,state,speed,kind)[0])
    for start in starts:
        res=minimize(obj,start,jac=True,method='L-BFGS-B',bounds=[(0,1)]*len(q),
                     options=dict(maxiter=45,ftol=1e-8,gtol=1e-5,maxls=25))
        out.append(dict(q=res.x.tolist(),cost=float(res.fun),success=bool(res.success),
                        iterations=int(res.nit),message=str(res.message)))
    best=min(out+[dict(q=q.tolist(),cost=incumbent,success=True,iterations=0,message='fixed incumbent')],key=lambda r:r['cost'])
    targets=v.LOW+np.array(best['q']).reshape(4,4)*v.SPAN
    metric,trace=replay(kind,speed,gains,targets,.002)
    fine,fine_trace=replay(kind,speed,gains,targets,.001)
    return dict(starts=out,selected=best,fixed_incumbent_objective=incumbent,
                metrics=metric,fine_metrics=fine,seconds=time.perf_counter()-tic,
                information='OFFLINE full true future material, not a causal comparison')

def main():
    best=json.loads((HERE/'tuning.json').read_text(encoding='utf8'))['selected'];gains=np.array(best['fixed']['gains'])
    # Exact scalar root independently checks the contact-table inverse.
    errors=[]
    for d in np.linspace(.00015,.0008,21):
        def psi(a):
            s=np.sqrt(v.CFG['lref']**2-a*a)
            return (v.CFG['lref']**2-v.CFG['lref']*s+a*a*np.log((v.CFG['lref']+s)/a))/(2*v.CFG['radius'])
        a=brentq(lambda a:psi(a)-d,1e-7,.01,xtol=1e-14)
        exact=v.CFG['length']*45000*np.pi*a*a/(4*v.CFG['radius'])
        errors.append(abs(v.contact(d,45000)-exact)/exact)
    geometry=dict(max_relative_force_table_error=float(max(errors)),
        max_surface_slope=v.CFG['surface_amplitude']*2*np.pi/v.CFG['surface_wavelength'],
        minimum_relative_curvature=1-v.CFG['radius']*v.CFG['surface_amplitude']*(2*np.pi/v.CFG['surface_wavelength'])**2,
        model='local parabolic curved-surface cylinder contact; small indentation/footprint approximation')
    write('geometry_audit.json',geometry);print('geometry audit',geometry,flush=True)
    # Replay existing target schedules: isolating integration from optimizer differences.
    integration={}
    for name in ['ramp_slow','ramp_fast','cycle_slow','cycle_fast']:
        p=HERE/(name+'_mpc_solves.json')
        while not p.exists():time.sleep(.2)
        solves=json.loads(p.read_text(encoding='utf8'))
        targets=np.array([v.LOW+np.array(s['q'])[:4]*v.SPAN for s in solves])
        kind=name.split('_')[0];speed=.003 if name.endswith('slow') else .006
        coarse,a=replay(kind,speed,gains,targets,.002);fine,b=replay(kind,speed,gains,targets,.001)
        original=np.loadtxt(HERE/(name+'_mpc.csv'),delimiter=',',skiprows=1)
        integration[name]=dict(coarse=coarse,fine=fine,
            cost_relative_difference=abs(fine['cost']/coarse['cost']-1),
            coarse_replay_state_max_abs=float(np.max(abs(a[:,:14]-original[:,:14]))))
    write('integration_audit.json',integration);print('integration audit complete',flush=True)
    independent=[]
    for name in ['ramp_fast','cycle_slow']:
        solves=json.loads((HERE/(name+'_mpc_solves.json')).read_text(encoding='utf8'))
        for idx in [len(solves)//4,len(solves)//2,3*len(solves)//4]:
            s=solves[idx];state=np.array(s['state']);q=np.array(s['q']);env=tuple(s['env'])
            kind=name.split('_')[0];speed=.003 if name.endswith('slow') else .006;t=s['t']
            scalar=lambda p:float(v.shooting(p,state,t,speed,env,kind,dt=.002)[0])
            res=minimize(scalar,q,method='SLSQP',bounds=[(0,1)]*len(q),
                         options=dict(maxiter=60,ftol=1e-8,eps=1e-5))
            independent.append(dict(case=name,t=t,original_fine_objective=scalar(q),
                independent_fine_objective=float(res.fun),success=bool(res.success),message=str(res.message),
                improvement_fraction=float(1-res.fun/scalar(q))))
            write('independent_solver_audit.json',independent)
    print('independent solver audit complete',flush=True)
    oracles={}
    for kind,speed in [('ramp',.006),('cycle',.003)]:
        name=kind+('_slow' if speed==.003 else '_fast')
        oracles[name]=oracle(kind,speed,gains);write('offline_oracle_audit.json',oracles)
        print('offline oracle',name,oracles[name]['metrics'],flush=True)
    # Exact deterministic rerun of one held-out causal optimized closed loop.
    r,a,solves=v.run('ramp',.006,'mpc',best['fixed'],save='reproduction_ramp_fast_mpc')
    old=np.loadtxt(HERE/'ramp_fast_mpc.csv',delimiter=',',skiprows=1)
    reproduction=dict(case='ramp_fast_mpc',trajectory_exact=bool(np.array_equal(a,old)),
                      max_absolute_difference=float(np.max(abs(a-old))),metrics=r,
                      timing_note='Solver timings measured again; timing equality is not required.')
    write('reproducibility.json',reproduction);print('reproduction',reproduction,flush=True)
    assert reproduction['trajectory_exact']
    assert max(i['cost_relative_difference'] for i in integration.values())<.005
    assert max(errors)<1e-5

if __name__=='__main__':main()
