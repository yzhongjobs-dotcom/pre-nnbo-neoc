"""Mechanism, selected-configuration and causal replay checks; no test tuning."""
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.linalg import expm

import id_baselines as I
import ctrl_baselines as C
from comparators import nominal_mpc
import physics_plant3 as P3
import scenes3 as S3


def checks(chosen=None):
    rng=np.random.default_rng(46);n=64;dt=.008;R=.04**2
    Y=rng.normal(0,.0002,n);zz=2300*Y+1.2+rng.normal(0,.04,n)
    x0=np.array([30.,0.]);p0=np.diag([300.**2,100.**2]);x=x0.copy();P=p0.copy()
    H=np.c_[1000*Y,np.ones(n)]
    for h,z in zip(H,zz):x,P,_=I.kf_step(x,P,h,z,R,np.eye(2),np.zeros((2,2)))
    J=np.linalg.inv(p0)+H.T@H/R;target=np.linalg.solve(J,np.linalg.solve(p0,x0)+H.T@zz/R)
    err=float(np.max(np.abs(x-target)));assert err<1e-7
    tr,var,sse=I.mhe_solve(Y,zz,R,qk=0.,qd=0.)
    ls=np.linalg.lstsq(H,zz,rcond=None)[0];err2=float(np.max(abs(tr[-1]-ls)));assert err2<1e-9
    # Finite-window arrival + history equals the recursive posterior exactly.
    qk,qd=3000.,.3;tr,var,_=I.mhe_solve(Y,zz,R,x0,p0,qk,qd,dt)
    x=x0.copy();P=p0.copy();Q=np.diag([qk*qk*dt/1e6,qd*qd*dt])
    for j,(h,z) in enumerate(zip(H,zz)):
        x,P,_=I.kf_step(x,P,h,z,R,np.eye(2),np.zeros((2,2)) if j==0 else Q)
    err3=float(np.max(abs(tr[-1]-x)));varerr=abs(var-P[0,0]*1e6)
    assert err3<1e-7 and varerr<1e-5
    th=np.array([6.,1.2,.1]);phi=np.array([1.,-4.,.02]);lf=2.;de=.03
    ref=th+phi*(lf-th@phi)/(1+phi@phi)
    assert np.allclose(I.xue_update(th,phi,lf,de),ref)
    assert np.array_equal(I.xue_update(np.zeros(3),phi,.0625,.03125),np.zeros(3))
    # Independently condense step responses and compare (27), not another
    # invocation of the implementation under test.
    do=(chosen or {}).get('DO-MPAC',dict(Np=15,q_over_r=1e6,eso_bandwidth=8.))
    horizon=do['Np'];ratio=do['q_over_r']
    g,poles=nominal_mpc(dt,horizon,ratio);A=np.array([[1.,dt],[0.,1.]]);B=np.array([dt*dt/2,dt]);c=np.array([1.,0.])
    Psi=np.array([c@np.linalg.matrix_power(A,j+1) for j in range(horizon)])
    Theta=np.array([[c@np.linalg.matrix_power(A,i-j)@B if i>=j else 0. for j in range(horizon)] for i in range(horizon)])
    gain=np.linalg.solve(Theta.T@Theta+1/ratio*np.eye(horizon),Theta.T@Psi)[0]
    err4=float(np.max(abs(g-gain)));assert err4<1e-8 and max(abs(poles))<1
    w=do['eso_bandwidth'];Ac=np.array([[-3*w,1.,0.],[-3*w*w,0.,1.],[-w**3,0.,0.]])
    poleerr=float(np.max(abs(np.linalg.eigvals(Ac)+w)));assert poleerr<.001
    # Nominal Rossi full controller convergence: static linear spring and
    # the same 2-pole actuator/LPF/delay, no surface or stiffness truth fed in.
    p=C.RossiController(**(chosen or {}).get('Rossi-control',dict(K_I=5.,A_E=2e-5)));p.force_set=8.;p.reset(P3.PLATFORM)
    y=.002;vel=0.;f=2.;sent=y;k=1000.;history=[]
    for j in range(1500):
        row=np.full(26,np.nan);row[[0,2,5,11]]=[j*dt,y,f,sent];history.append(row)
        cmd=p(SimpleNamespace(n=len(history),log=np.asarray(history)))
        cmd=float(np.clip(cmd,sent-.25*dt,sent+.25*dt))
        # Match P3's finalized historical command column after the callback.
        history[-1][11]=cmd
        # With ndel=1 P3 shifts [older,sent] to [sent,new] before the
        # interval: the actuator sees sent, not the command from two ticks ago.
        for substep in range(16):
            if substep%2==0:f+=p.afc*(k*y-f)
            a=(sent-y-.019*vel)/((.019-P3.ACTUATOR_T2)*P3.ACTUATOR_T2)
            vel+=.0005*a;y+=.0005*vel
        sent=cmd
    staticerr=abs(f-8.)
    if chosen is None:assert staticerr<.05
    # The observation object masks truth, including all columns called out in
    # the task; ForceSensing sees the same masked object in every policy.
    obs=P3.Obs(np.ones((2,26)),2)
    assert np.isnan(obs.log[:,[3,4,6,9,18,20,21,22]]).all()
    return dict(kalman_constant_rls_max_error=err,mhe_constant_ls_max_error=err2,
                mhe_finite_vs_recursive_max_error=err3,mhe_endpoint_variance_error=varerr,
                xue_printed_law=True,dompac_gain_max_error=err4,dompac_nominal_poles=poles.tolist(),
                eso_poles= np.linalg.eigvals(Ac).real.tolist(),eso_pole_error=poleerr,
                rossi_static_force_error=staticerr,rossi_static_converged=bool(staticerr<.05),
                rossi_static_tail_rmse=float(np.sqrt(np.mean((np.asarray(history)[-250:,5]-8.)**2))),
                rossi_static_khat=p.k_hat,truth_mask=True)


def source_replay_check():
    """Replay the source policy; clipped requests must reproduce sent commands."""
    import baseline_compare as B
    data=B.load('baseline_test.json');results={}
    for method in ('Huber-D4','D4'):
        source=next(r for r in data['rows'] if r['method']==method and r['cond']=='combo')
        z=np.load(B.OUT/source['log_file']);log=np.full((len(z['log']),26),np.nan);log[:,z['cols']]=z['log']
        observed=z['observed']
        sc=B.scene(B.key(source),B.TEST_BASE,data['start_gap']);pol=B.build(method,data['chosen'][method],sc)
        pol.reset(dict(P3.PLATFORM,tau_a=sc['tau_a_hat']));diff=[]
        for n in range(1,len(log)+1):
            u=pol(B.source_observation(log,observed,n));sent=observed[n-1,B.OBS_COLS.index(11)];limit=P3.PLATFORM['VMAX']*P3.PLATFORM['Tc']
            diff.append(float(np.clip(u,sent-limit,sent+limit)-log[n-1,11]))
        err=float(np.max(np.abs(diff)));assert err<1e-10,(method,err)
        results[method]=err
    return results


if __name__=='__main__':
    import sys
    dest=Path(__file__).resolve().parent/'results'/'baseline_checks.json'
    if '--replay' in sys.argv:
        out=json.loads(dest.read_text(encoding='utf-8'));out['source_replay_command_error']=source_replay_check()
    elif '--selected' in sys.argv:
        chosen=json.loads((dest.parent/'baseline_dev.json').read_text(encoding='utf-8'))['chosen']
        out=json.loads(dest.read_text(encoding='utf-8'));out['selected_configuration_checks']=checks(chosen)
    else:
        previous=json.loads(dest.read_text(encoding='utf-8')) if dest.exists() else {}
        out=checks()
        for k in ('default_method_command_difference','source_replay_command_error','selected_configuration_checks'):
            if k in previous:out[k]=previous[k]
    dest.write_text(json.dumps(out,indent=2,default=str),encoding='utf-8');print(out)
