"""Open-loop protocol. Run with --workers 8 (default). No final seeds allowed.
Transient trajectories are stored only in --cache outside delivery results.
"""
import argparse, json, pickle, platform, time
from concurrent.futures import ProcessPoolExecutor
from copy import deepcopy
import numpy as np
from scipy.optimize import minimize
from scipy.linalg import cho_factor,cho_solve
from step3b_core import *


def collect(scenes,cache,workers):
    if cache.exists():
        data=pickle.loads(cache.read_bytes())
        assert [(e['family'],e['sid'],e['code']) for e in data]==[(s[0],s[1],s[4]) for s in scenes]
        return data
    t=time.perf_counter(); data=[]
    data=file_pool('step3b_core','data_episode',scenes,cache.parent/'workers',workers)
    print('data',cache.stem,len(data),'seconds',round(time.perf_counter()-t),flush=True)
    cache.write_bytes(pickle.dumps(data,protocol=5))
    return data


def selection_score(pred,data):
    # equal family weight; scene is the independent unit, no loss/contact exclusion
    values={f:[] for f in FAMS}
    for e,(y,v) in zip(data,pred):
        mask=e['t']>=1
        values[e['family']].append(float(np.mean((y[mask]-e['target'][mask])**2)))
    return float(np.mean([np.mean(v) for v in values.values()]))


def select_imm(data,third=True):
    rng=np.random.default_rng(209950 if third else 209951); trials=[]
    best=None
    for i in range(61):
        if i==0:
            kw=dict(qs=CE.FINAL['q_smooth'],qm=CE.FINAL['q_maneuver'],stay=CE.FINAL['p_stay'])
            if third: kw['qj']=100.
        else:
            kw=dict(qs=10**rng.uniform(-9,-3),qm=10**rng.uniform(-3,1),stay=float(rng.uniform(.9,.999)))
            if third: kw['qj']=10**rng.uniform(-1,5)
        model=MultiIMM(**kw,third=third)
        pred=[seq_imm(model,e) for e in data]
        score=selection_score(pred,data)
        trials.append(dict(params=kw,mse_m2=score))
        if best is None or score<best[0]: best=(score,kw)
        if i%10==0: print('IMM',third,i,'best um',np.sqrt(best[0])*1e6,flush=True)
    return best[1],trials


def fit_gp(train):
    # Unsupervised covariance likelihood, same causal training observations.
    # One window per scene, balanced across families/profiles, no validation labels.
    rng=np.random.default_rng(209952); xx=[]
    for e in train:
        X,age=windows(e); ids=np.flatnonzero((e['t']>=1)&(age>=W))
        if len(ids):
            j=int(rng.choice(ids)); x=X[j]; xx.append(x-x.mean())
    X=np.array(xx)
    # metres; covariance models displacement variation, not absolute reset level
    bounds=[(np.log(1e-5),np.log(.02)),(np.log(.015),np.log(2.)),
            (np.log(1e-6),np.log(.02)),(np.log(1.5),np.log(7.)),
            (np.log(.1),np.log(4.)),(np.log(1e-7),np.log(.002))]
    t=np.arange(W)*DT; d=t[:,None]-t[None,:]
    count=0
    def objective(theta):
        nonlocal count
        count+=1
        a,ell,ap,period,lp,noise=np.exp(theta)
        K=a*a*np.exp(-d*d/(2*ell*ell))+ap*ap*np.exp(-2*np.sin(np.pi*d/period)**2/lp**2)+noise*noise*np.eye(W)
        cf=cho_factor(K); sol=cho_solve(cf,X.T)
        return float(.5*(np.mean(np.sum(X.T*sol,axis=0))+2*np.log(np.diag(cf[0])).sum()))
    x0=np.log([.003,.15,.001,4.,1.,.00001])
    # Powell maxfev gives a strict 61 training objective evaluations, no restarts.
    fit=minimize(objective,x0,method='Powell',bounds=bounds,options=dict(maxfev=61,xtol=.01,ftol=.001))
    return fit.x.tolist(),dict(evaluations=count,success=bool(fit.success),message=str(fit.message),train_windows=len(X),nll=float(fit.fun))


def fit_net(train,val):
    rng=np.random.default_rng(209953)
    Xs=[]; Ys=[]
    for e in train:
        X,age=windows(e); mask=(e['t']>=.5)
        Xs.append(X[mask]*1000); Ys.append((e['target'][mask]-e['z'][mask])*1000)
    X=np.concatenate(Xs);Y=np.concatenate(Ys)
    net=Net(); opt=dict(t=0,m=[np.zeros_like(x) for x in net.w],v=[np.zeros_like(x) for x in net.w])
    VX=[];VY=[]
    for e in val:
        x,_=windows(e);m=e['t']>=1;VX.append(x[m]*1000); VY.append((e['target'][m]-e['z'][m])*1000)
    vx=np.concatenate(VX);vy=np.concatenate(VY)
    best=(float('inf'),None,0); trace=[]; stale=0
    for epoch in range(61):
        # Fixed, equal update budget per epoch, resampled training windows.
        ids=rng.integers(0,len(Y),16384)
        losses=[]
        for start in range(0,len(ids),256):
            ix=ids[start:start+256]; losses.append(net.step(X[ix],Y[ix],opt,3e-4))
        _,out=net.forward(vx);lv=np.clip(out[:,1],-14,8)
        vnll=float(np.mean(.5*(lv+(out[:,0]-vy)**2*np.exp(-lv))))
        pred=[net.batch(windows(e)[0],e['z']) for e in val]
        score=selection_score(pred,val)
        trace.append(dict(epoch=epoch+1,training_nll=float(np.mean(losses)),validation_nll=vnll,validation_rmse_um=np.sqrt(score)*1e6))
        if vnll<best[0]-1e-5: best=(vnll,deepcopy(net.w),epoch+1);stale=0
        else: stale+=1
        if epoch%10==0: print('net',epoch+1,trace[-1],flush=True)
        if stale>=12: break
    net.w=best[1]
    return dict(hidden=48,weights=[w.tolist() for w in net.w],epoch=best[2],parameters=sum(w.size for w in net.w)),trace


def eval_pred(name,model,e):
    if name=='P0':
        # Recompute uncertainty for old transient caches as well: forecast
        # process variance is needed, even though FINAL returns only a mean.
        baseline=MultiIMM(CE.FINAL['q_smooth'],CE.FINAL['q_maneuver'],stay=CE.FINAL['p_stay'],third=False)
        _,v=seq_imm(baseline,e)
        return e['p0'],v
    if isinstance(model,MultiIMM): return seq_imm(model,e)
    return model.batch(windows(e)[0],e['z'])


def metrics(err,var):
    if len(err)==0: return dict(n_ticks=0)
    return dict(n_ticks=len(err),rmse_um=float(np.sqrt(np.mean(err**2))*1e6),p95_um=float(np.percentile(abs(err),95)*1e6),
                coverage_1=float(np.mean(abs(err)<=np.sqrt(np.maximum(var,0)))),coverage_2=float(np.mean(abs(err)<=2*np.sqrt(np.maximum(var,0)))))


def evaluate(models,data):
    out={}
    for name,model in models.items():
        groups={};scenes=[]
        for e in data:
            y,v=eval_pred(name,model,e);mask=e['t']>=1;err=y-e['target']
            keys=[('all','all'),(e['family'],'all')]+[('all',ph) for ph in ('breathing','pre50ms','ramp','hold','retract')]+[(e['family'],ph) for ph in ('breathing','pre50ms','ramp','hold','retract')]
            for fam,ph in keys:
                m=mask if ph=='all' else mask&(e['phase']==ph)
                if np.any(m):
                    g=groups.setdefault(fam+'|'+ph,[[],[]]);g[0].append(err[m]);g[1].append(v[m])
            de=np.diff(y)-np.diff(e['target']);dm=mask[1:]&~e['reset'][1:]
            scenes.append(dict(family=e['family'],sid=e['sid'],code=e['code'],**metrics(err[mask],v[mask]),increment_rmse_um=float(np.sqrt(np.mean(de[dm]**2))*1e6)))
        out[name]=dict(groups={k:metrics(np.concatenate(x),np.concatenate(v)) for k,(x,v) in groups.items()},scenes=scenes)
    return out


def benchmark(models,data):
    # Serial hot timing, wrapper/update/prediction included; no process contention.
    out={}
    e=data[0]
    for name,model in models.items():
        samples=[]
        legacy=None
        if name=='P0': legacy=IMMTracker(Tc=DT,q_smooth=CE.FINAL['q_smooth'],q_maneuver=CE.FINAL['q_maneuver'],p_stay=CE.FINAL['p_stay'])
        buffer=[]
        for z,R,reset in zip(e['z'],e['R'],e['reset']):
            t=time.perf_counter_ns()
            if legacy:
                legacy.R=R
                if reset: legacy.reset(z)
                else: legacy.update(z)
                s,sd=legacy.predict(L);value=s+TAU*sd
            elif isinstance(model,MultiIMM):
                if reset: model.reset(z,R)
                else: model.update(z,R)
                value=model.output()
            else:
                if reset: buffer=[]
                buffer.append(z);buffer=buffer[-W:]
                x=np.pad(np.array(buffer),(max(0,W-len(buffer)),0),mode='edge')[-W:]
                value=model.batch((x-z)[None,:],np.array([z]))
            samples.append((time.perf_counter_ns()-t)/1e6)
        out[name]=dict(p50_ms=float(np.percentile(samples[100:],50)),p99_ms=float(np.percentile(samples[100:],99)),n=len(samples)-100)
    return out


def write_md(out):
    lines=['# Step 3b: causal motion prediction','',
           'All prediction errors below are in μm. Predictions are evaluated on the SAME FINAL closed-loop measurements; no truth enters predictor inputs. The target sign is eta=surface+reset constant=-s_eq+constant. Oracle uses future physical surface at the current observed x, not future x.', '',
           '| Model | ID RMSE | ID P95 | OOD RMSE | OOD P95 | ID 1σ/2σ | P50/P99 ms |', '|---|---:|---:|---:|---:|---|---|']
    for name in out['evaluation']['id']:
        a=out['evaluation']['id'][name]['groups']['all|all'];b=out['evaluation']['ood'][name]['groups']['all|all'];t=out['timing'][name]
        lines.append(f"| {name} | {a['rmse_um']:.2f} | {a['p95_um']:.2f} | {b['rmse_um']:.2f} | {b['p95_um']:.2f} | {a['coverage_1']:.3f}/{a['coverage_2']:.3f} | {t['p50_ms']:.3f}/{t['p99_ms']:.3f} |")
    for split in ('id','ood'):
        lines+=['',f'## {split}: family and stage errors','', '| Group | '+ ' | '.join(out['evaluation'][split])+' |','|---|'+'---:|'*len(out['evaluation'][split])]
        keys=out['evaluation'][split]['P0']['groups'].keys()
        for key in keys:
            if key=='all|all' or key.startswith('all|') or key.endswith('|all'):
                cells=[f"{out['evaluation'][split][n]['groups'][key]['rmse_um']:.2f}/{out['evaluation'][split][n]['groups'][key]['p95_um']:.2f}" for n in out['evaluation'][split]]
                lines.append('| '+key.replace('|',' / ')+' | '+' | '.join(cells)+' |')
    lines+=['','## Protocol and boundaries','',out['protocol'], '',
            'Selection is frozen BEFORE diagnostic evaluation. P0r is a budget-matched retuned two-CA sensitivity reference (not a new architecture); P0 stays FINAL. P1 has stationary / CA / constant-jerk modes. P2 fits SE+periodic covariance on training observations, profiles the constant mean, and predicts the exact linear functional with derivative covariance. P3 is a 63-sample temporal convolution and 1x1 Gaussian head, 3170 parameters, trained with Gaussian NLL. Early stopping is by validation NLL; closed-loop candidate selection is validation scene MSE, never diagnostic/OOD outcomes.', '',
            'Calibration is tick coverage, descriptive only. P0/P1 covariances model state uncertainty but omit contact-model and stiffness bias; none is calibrated on evaluation. Reset offsets and gain error can accumulate: absolute physical-target error includes these, so scene increment RMSE is also retained in JSON. Equal causal information does not imply equal supervision: P3 uses training physical-target labels, P2 uses measurement likelihood, and IMM uses validation labels for selection. Computational budgets are explicit, rather than pretending those training operations are equivalent.', '',
            'Training includes 100 scenes/family with 20 each of smoothstep, cosine, linear ramp, minimum jerk, and smoothstep with varied breathing; nonzero profiles include a second respiratory harmonic. ID is the original cohort, 20/family. OOD has trapezoidal velocity, two-stage pushes and phase-continuous frequency jumps, 20/family/profile, with matching noise/control structure. Stage labels are temporal event intervals; breathing means outside event intervals and is not necessarily free of stiffness transitions in S2/T3m. All ticks after 1 s are retained, including post-loss resets; incomplete episodes are counted in metadata.', '',
            'No variance scaling of feed-forward is applied. Multiplying an absolute prediction by a time-varying confidence gain can create offset-dependent command jumps; a risk tradeoff and state-consistent increment rule would need separate design. Hot serial timings include prediction/update, exclude stiffness estimation, plant and cold compilation, and are not real-time certification.']
    (RES/'step3b_predictors.md').write_text('\n'.join(lines)+'\n',encoding='utf8')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=8);parser.add_argument('--cache',type=Path,required=True);args=parser.parse_args();args.cache.mkdir(exist_ok=True,parents=True)
    started=time.time()
    train=collect(make_scenes(200000,100,(0,1,2,3,4)),args.cache/'train.pkl',args.workers)
    val=collect(make_scenes(190000,10),args.cache/'val.pkl',args.workers)
    select=val
    p1,trial1=select_imm(select)
    p0r,trial0=select_imm(select,third=False)
    p2,gptrace=fit_gp(train);print('GP',gptrace,flush=True)
    p3,nntrace=fit_net(train,select)
    selected={'P1':p1,'P0r':p0r,'P2':p2,'P3':p3}
    (RES/'step3b_predictors.json').write_text(json.dumps(dict(selected=selected)),encoding='utf8')
    _,models=load_models();models={'P0':None,**models,'P0r':MultiIMM(**p0r,third=False)}
    scores={n:selection_score([eval_pred(n,m,e) for e in select],select) for n,m in models.items()}
    # Always select the two best new candidates, even if P0 is better.
    chosen=sorted(('P1','P2','P3'),key=lambda n:scores[n])[:2]
    print('frozen closed-loop candidates',chosen,'validation',scores,flush=True)
    diag=collect(make_scenes(110000,20),args.cache/'diag.pkl',args.workers)
    ood=[]
    for code in (5,6,7):
        ood+=collect(make_scenes(115000+(code-5)*5000,20,(code,)),args.cache/f'ood{code}.pkl',args.workers)
    # Note last OOD base125000+4000 stays below130000.
    out=dict(selected=selected,validation_scores_m2=scores,closed_loop_candidates=chosen,
             search=dict(P1=trial1,P0r=trial0,P2=gptrace,P3=nntrace),
             evaluation=dict(id=evaluate(models,diag),ood=evaluate(models,ood)),timing=benchmark(models,diag),
             ood_by_profile={str(code):evaluate(models,[e for e in ood if e['code']==code]) for code in (5,6,7)},
             dataset={n:dict(n_scenes=len(ds),incomplete=sum(not e['full_length'] for e in ds),n_ticks=sum(len(e['z']) for e in ds)) for n,ds in [('train',train),('validation',select),('id',diag),('ood',ood)]},
             protocol='Training seeds 200000+ (100/family); selection 190000+ (exactly10/family); ID110000+ (20/family), OOD115000+/120000+/125000+ (20/family/profile, max seed129019). No final130000+ cohort accessed. The OOD125000 block is separate from FINAL and stops before130000. Each new model has at most61 validation looks (GP one); P0r61, P1 61, NN≤61epochs with12epoch patience and16384 sampled windows/epoch. All models see the same causal eta and measurement-variance sequences; GP/NN use63 samples, IMM recursively summarizes history. No evaluation tuning.',
             machine=dict(python=sys.version,platform=platform.platform(),numpy=np.__version__,threads=1),elapsed_s=time.time()-started)
    (RES/'step3b_predictors.json').write_text(json.dumps(out,indent=2,allow_nan=False),encoding='utf8');write_md(out)
    print('completed',out['elapsed_s'],flush=True)

if __name__=='__main__': main()
