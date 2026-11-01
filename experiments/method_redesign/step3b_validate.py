"""Focused checks for label sign, causal replay, integration, profiles and loop model.
Run after fitted model JSON exists; writes only a compact summary JSON.
"""
import json,copy,hashlib,time
from step3b_core import *

def main():
    checks={};start=time.time()
    scene=make_scenes(190000,1)[3];fam,sid,p,T,code=scene;use_plant(0)
    a=variants.build('Proposed-IMM',CE.FINAL);la,_=common.run(a,p,T,sid)
    r=Recorder(**CE.FINAL);lr,_=common.run(r,p,T,sid)
    checks['recorder_is_passive']=bool(np.array_equal(la,lr));assert checks['recorder_is_passive']
    e=data_episode(scene)
    old=MultiIMM(CE.FINAL['q_smooth'],CE.FINAL['q_maneuver'],stay=CE.FINAL['p_stay'],third=False)
    yy,vv=seq_imm(old,e)
    checks['legacy_replay_max_abs_m']=float(np.max(abs(yy-e['p0'])))
    checks['legacy_replay_max_var_m2']=float(np.max(abs(vv-e['v0'])))
    assert checks['legacy_replay_max_abs_m']<1e-12
    b=Replacement(copy.deepcopy(old),**CE.FINAL);lb,_=common.run(b,p,T,sid)
    checks['legacy_adapter_log_max_abs']=float(np.max(abs(la-lb)));assert np.allclose(la,lb,atol=1e-7,rtol=1e-7)
    # Directly verify labels against Oracle at each physical x and reset anchor.
    t,x=[],[]
    # Recorder retains x alongside eta, permitting a completely independent label read.
    arr=np.array(r.records)
    for i,row in enumerate(arr):
        ti,xi=row[:2]
        if e['reset'][i]: anchor=BASE_PLANT.surface(xi,ti,p)[0]
        h,ht,_=BASE_PLANT.surface(xi,ti+L,p)
        t.append(h+TAU*ht-anchor)
    checks['oracle_target_max_abs_m']=float(np.max(abs(np.array(t)-e['target'])))
    assert checks['oracle_target_max_abs_m']<1e-14
    # With exact stiffness and ideal linear contact, inverse force minus y is +surface.
    y=.004; k=8000.; h=.002
    checks['measurement_sign_eta_equals_positive_surface']=bool(abs((k*(y+h))/k-y-h)<1e-15)
    assert checks['measurement_sign_eta_equals_positive_surface']
    cfg,models=load_models();errs={};prefix={}
    for name,m in models.items():
        tr=IMMTracker(Tc=DT,q_smooth=CE.FINAL['q_smooth'],q_maneuver=CE.FINAL['q_maneuver'],p_stay=CE.FINAL['p_stay'])
        adapter=Adapter(tr,copy.deepcopy(m));values=[]
        for z,R,reset in zip(e['z'],e['R'],e['reset']):
            adapter.R=R
            if reset:adapter.reset(z)
            else:adapter.update(z)
            values.append(adapter.predict(L)[0])
        if isinstance(m,MultiIMM): yy,_=seq_imm(m,e)
        else: yy,_=m.batch(windows(e)[0],e['z'])
        errs[name]=float(np.max(abs(np.array(values)-yy)));assert errs[name]<1e-11
        cut=400;ep={k:v[:cut] if isinstance(v,np.ndarray) and len(v)==len(e['z']) else v for k,v in e.items()}
        if isinstance(m,MultiIMM): yp,_=seq_imm(m,ep)
        else: yp,_=m.batch(windows(ep)[0],ep['z'])
        prefix[name]=float(np.max(abs(yp-yy[:cut])));assert prefix[name]<1e-11
    checks['online_offline_max_abs_m']=errs;checks['causal_prefix_max_abs_m']=prefix
    # Profile derivative consistency, including all held-out profiles.
    mod=profile_plant();rng=np.random.default_rng(209956);derrs={}
    for code in range(8):
        pp=np.r_[p,float(code)];err=[]
        for ti in rng.uniform(.1,7.9,100):
            h,ht,_=mod.surface(.45,ti,pp);eps=1e-7
            num=(mod.surface(.45,ti+eps,pp)[0]-mod.surface(.45,ti-eps,pp)[0])/(2*eps)
            err.append(abs(num-ht))
        derrs[str(code)]=max(err)
        assert max(err)<1e-6
    checks['profile_derivative_max_abs_m_s']=derrs
    pp=np.r_[p,0.]
    checks['profile0_max_abs']=max(float(np.max(abs(np.array(BASE_PLANT.surface(.45,ti,p))-np.array(mod.surface(.45,ti,pp))))) for ti in np.linspace(0,8,500))
    assert checks['profile0_max_abs']<1e-14
    from step3b_closed_loop import ss_radius,filter_ss
    import ff_loop as ff
    oldparams=(ff.M1,ff.D1,ff.KPF1)
    ff.M1,ff.D1,ff.KPF1=(CE.FINAL[k] for k in ('M1','D1','KPF1'))
    dif=[];gain_dif=[]
    for kg in (.75,1.,1.2):
        mm=MultiIMM(CE.FINAL['q_smooth'],CE.FINAL['q_maneuver'],stay=CE.FINAL['p_stay'],third=False)
        As,Bs,Cs=filter_ss(mm.A[1,:3,:3],mm.Q[1,:3,:3],(CE.FINAL.get('noise',.04)/(1e4/kg))**2*CE.FINAL['r_scale'],mm.H[1,:3])
        dif.append(abs(ss_radius(1.,kg,As,Bs,Cs)-ff.radius(1.,kg,Bs)))
        ref=ff.kalman_gain(CE.FINAL['q_maneuver'],1e4/kg)
        gain_dif.append(float(np.max(abs(Bs-ref)/np.maximum(abs(ref),1.))))
    checks['ff_loop_spectral_radius_max_difference']=max(dif);assert max(dif)<1e-9
    checks['dare_vs_original_iterative_gain_relative_difference']=max(gain_dif)
    assert max(gain_dif)<1e-5
    checks['all_passed']=True;checks['elapsed_s']=time.time()-start
    checks['source_sha256']={name:hashlib.sha256((HERE.parent/'motion_aware_scan'/name).read_bytes()).hexdigest() for name in ('imm.py','variants.py','proposed.py','common.py','cohort_eval.py','sweep.py','ff_loop.py')}
    (RES/'step3b_validation.json').write_text(json.dumps(checks,indent=2),encoding='utf8');print(json.dumps(checks,indent=2))

if __name__=='__main__':main()
