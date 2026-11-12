"""Fixed trajectory-correction closed loops. New files only, no tuning or plant options.
python -X utf8 trajectory_loop.py --project-root REPO --output-root MIRROR
Equivalence precedes episodes; exact Huber regression precedes candidates.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[_key] = '1'
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
import sys
sys.dont_write_bytecode = True
import argparse
import json
import math
import tempfile
import time
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# Discover project before imports, including Windows worker imports of this file.
if '--project-root' in sys.argv:
    os.environ['TRAJECTORY_LOOP_PROJECT'] = sys.argv[sys.argv.index('--project-root')+1]
PROJECT = Path(os.environ.get('TRAJECTORY_LOOP_PROJECT', Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(PROJECT/'experiments/rossi_failure_study'))
# Never modify existing JIT caches in the repository.
if 'TRAJECTORY_LOOP_CACHE' not in os.environ:
    os.environ['TRAJECTORY_LOOP_CACHE'] = tempfile.mkdtemp(prefix='trajectory_loop_numba_')
os.environ['NUMBA_CACHE_DIR'] = os.environ['TRAJECTORY_LOOP_CACHE']

import numpy as np
import corrected_comparison as CC
import memory_segment_eval as SE
import overlay_compare as O
import dual_channel_replay as D
import dual_channel_v3_replay as V3
P3, B, M = SE.P3, O.B, SE.M
CASES = tuple(range(24))+tuple(range(36,42))
FAST = (36,37,39,40)
GAINS = {36:12.,37:8.,38:3.,39:1.5,40:3.}
SEGMENTS = ('plateau','slow','rise','fall','reversal','event')
METRICS = ('rmse','peak','over_15','loss_time','osc','ripple')


class OnlineFilter:
    """Truth-free, tick-for-tick port of dual_channel_v3_replay.run."""
    def __init__(self):
        dt = V3.TC
        self.tau = 11.5*dt
        w = 2*math.pi*.25
        self.A = np.eye(5)
        self.A[0,1] = dt
        cw, sw = math.cos(w*dt), math.sin(w*dt)
        self.A[2,2], self.A[2,3], self.A[3,2], self.A[3,3] = cw, sw, -sw, cw
        self.ct, self.st = math.cos(w*self.tau), math.sin(w*self.tau)
        self.Js = np.array([[dt**3/3,dt**2/2],[dt**2/2,dt]])
        self.resets = 0
        self.reset()

    def reset(self):
        self.previous_n = None
        self.y_sums, self.f_sums, self.fs = [0.], [0.], []
        self.active, self.last = False, -1000
        self.out = self.sd = self.g = np.nan

    def observe(self, n, y, force):
        if self.previous_n is not None and n != self.previous_n+1:
            self.resets += 1
            self.reset()
        self.previous_n = n
        # Same sequential cumsum and subtraction as D.trailing.
        self.y_sums.append(self.y_sums[-1]+float(y))
        self.f_sums.append(self.f_sums[-1]+float(force))
        j = len(self.fs)
        lo_y, lo_f = max(0,j-23), max(0,j-5)
        yw = (self.y_sums[j+1]-self.y_sums[lo_y])/(j+1-lo_y)
        fs = (self.f_sums[j+1]-self.f_sums[lo_f])/(j+1-lo_f)
        self.fs.append(fs)
        return j, yw, fs

    def step(self, n, y, force, k_w, sigma):
        j, yw, fs = self.observe(n,y,force)
        finite = np.isfinite(k_w)
        rel = sigma/k_w if finite else np.inf
        valid = finite and rel <= 1.
        m = float(np.log(np.clip(k_w,200.,3e5))) if finite else np.nan
        R = 24*rel**2
        if finite:
            self.last = j
        if fs < 1. or j-self.last > 25 or not np.isfinite(fs+y):
            self.active = False
            return self.out,self.sd,self.g
        if not self.active:
            if not valid:
                return self.out,self.sd,self.g
            C0 = fs/math.exp(m)
            self.y0 = y
            b0 = 1./3.
            self.x = np.array([C0,0.,0.,0.,b0])
            self.P = np.diag([C0*C0*R/24,(C0*.5)**2,(C0*.5)**2,(C0*.5)**2,(b0*1.1)**2])
            self.active = True
        else:
            x, P = self.x,self.P
            c2 = max(x[0]+x[2]+x[4]*(y-self.y0),1e-6)**2
            x = self.A@x
            P = self.A@P@self.A.T
            P[:2,:2] += 1.0*c2*self.Js
            P[2,2] += 1.0*c2*V3.TC
            P[3,3] += 1.0*c2*V3.TC
            P[4,4] += 1e-7*x[4]*x[4]
            if valid and j >= 24 and self.fs[j-12] >= 1.:
                Cobs = self.fs[j-12]/math.exp(m)
                H = np.array([1.,-self.tau,self.ct,-self.st,yw-self.y0])
                Rc = Cobs*Cobs*R
                PH = P@H
                S = float(H@PH)+Rc
                K = PH/S
                x = x+K*(Cobs-H@x)
                I = np.eye(5)-np.outer(K,H)
                P = I@P@I.T+Rc*np.outer(K,K)
            self.x,self.P = x,P
        H0 = np.array([1.,0.,1.,0.,y-self.y0])
        C = max(float(H0@self.x),1e-7)
        self.out = math.log(fs)-math.log(C)
        self.sd = math.sqrt(max(float(H0@self.P@H0),0.))/C
        self.g = 1./self.x[4] if self.x[4] > 0 else np.nan
        return self.out,self.sd,self.g


def equivalence(project):
    checks = []
    for index in (0,36,39):
        path = project/'experiments/rossi_failure_study/results/huber_windows'/f'case_{index:02d}.npz'
        with np.load(path,allow_pickle=False) as source:
            log,raw,Y = source['log'],source['raw'],source['Y']
        data = D.preprocess(log[:,0],log[:,2],log[:,4],raw,Y)
        case = dict(data=data,y=log[:,2],yw=D.trailing(log[:,2],24))
        reference = V3.run(case,(1.0,1.0,1e-7),.25)[0]
        rel = np.sqrt(data[4]/24)
        filt = OnlineFilter()
        actual = np.array([filt.step(j+1,log[j,2],log[j,4],raw[j],rel[j]*raw[j])[0]
                           for j in range(len(raw))])
        same = np.array_equal(np.isfinite(reference),np.isfinite(actual))
        finite = np.isfinite(reference)&np.isfinite(actual)
        delta = float(np.max(np.abs(reference[finite]-actual[finite])))
        check = dict(index=index,max_abs_difference=delta,same_finite_pattern=same,
                     n=len(raw),met=bool(same and delta < 1e-9))
        checks.append(check)
        print('Equivalence '+json.dumps(check),flush=True)
        if not check['met']:
            raise RuntimeError('Online/replay equivalence failed; no closed loops run')
    return checks


class ReturnedEstimate:
    def __init__(self,original,arm,index):
        self.original,self.arm,self.index = original,arm,index
        self.filter = DelagCorrection(force_ratio_factor=arm == 'delag') if arm in ('delag', 'delag2') else OnlineFilter()
        self.t,self.n,self.k,self.kw,self.g = [],[],[],[],[]
        self.beta = []
        self.force_ratio = []

    def __call__(self,o,n,policy):
        k_w,sigma,chi = self.original(o,n,policy)
        k_c,gain = k_w,np.nan
        y,force = float(o.log[n-1,2]),float(o.log[n-1,5])
        beta = np.nan
        if self.arm in ('delag', 'delag2'):
            k_c,beta = self.filter.step(n,y,force,k_w,sigma)
            self.force_ratio.append(self.filter.last_force_ratio)
        elif self.arm == 'traj':
            out,_,gain = self.filter.step(n,y,force,k_w,sigma)
            if np.isfinite(k_w) and np.isfinite(out):
                k_c = float(np.clip(np.exp(out),200.,3e5))
        elif self.arm == 'traj_priv':
            j,yw,fs = self.filter.observe(n,y,force)
            gain = GAINS[self.index]
            if np.isfinite(k_w) and j >= 12:
                Cw = self.filter.fs[j-12]/k_w
                denominator = Cw+(y-yw)/gain
                if denominator < .3*Cw:
                    denominator = .3*Cw
                k_c = float(fs/denominator)
        self.t.append(float(o.log[n-1,0]))
        self.n.append(n)
        self.k.append(k_c)
        self.kw.append(k_w)
        self.g.append(gain)
        self.beta.append(beta)
        if (not np.isfinite(k_w) or self.arm == 'huber'
                or (self.arm in ('delag', 'delag2') and (k_c == k_w or not np.isfinite(sigma)))):
            return k_w,sigma,chi
        return k_c,sigma*k_c/k_w,chi


class DelagCorrection:
    """Window level plus bounded lag correction; return the pre-update gain."""
    def __init__(self, force_ratio_factor=True):
        self.force_ratio_factor = force_ratio_factor
        self.resets = 0
        self.reset()

    def reset(self):
        self.previous_n = None
        self.y_sums,self.f_sums = [0.],[0.]
        self.forces,self.cw,self.y24,self.rel = [],[],[],[]
        self.Sxx = self.Sxy = self.Sxn = self.beta = 0.

    def step(self,n,y,force,k_w,sigma):
        if self.previous_n is not None and n != self.previous_n+1:
            self.resets += 1
            self.reset()
        self.previous_n = n
        j = len(self.forces)
        self.forces.append(force)
        self.y_sums.append(self.y_sums[-1]+y)
        self.f_sums.append(self.f_sums[-1]+force)
        lo6,lo24 = max(0,j-5),max(0,j-23)
        y6 = (self.y_sums[j+1]-self.y_sums[lo6])/(j+1-lo6)
        y24 = (self.y_sums[j+1]-self.y_sums[lo24])/(j+1-lo24)
        F6 = (self.f_sums[j+1]-self.f_sums[lo6])/(j+1-lo6)
        F24 = (self.f_sums[j+1]-self.f_sums[lo24])/(j+1-lo24)
        self.last_force_ratio = F6/F24 if F24 != 0. else np.nan
        finite = np.isfinite(k_w)
        Cw = F24/k_w if finite else np.nan
        rel = sigma/k_w if finite else np.nan
        self.cw.append(Cw)
        self.y24.append(y24)
        self.rel.append(rel)
        beta_used,k_c = self.beta,k_w
        if finite and np.isfinite(sigma) and j >= 23 and F24 >= 1.:
            x = y6-y24
            r = np.clip(1.+beta_used*x/Cw,.5,2.)
            if self.force_ratio_factor:
                k_c = float(np.clip(k_w*(F6/F24)/r,200.,3e5))
            elif beta_used != 0.:
                k_c = float(np.clip(k_w/r,200.,3e5))
        if (j >= 24 and np.isfinite(Cw) and np.isfinite(self.cw[j-24])
                and rel <= 1. and self.rel[j-24] <= 1.
                and all(f >= 1. for f in self.forces[j-24:j+1])):
            dx = y24-self.y24[j-24]
            dC = Cw-self.cw[j-24]
            v = Cw**2*rel**2+self.cw[j-24]**2*self.rel[j-24]**2
            lam = math.exp(-V3.TC/4.)
            self.Sxx = lam*self.Sxx+dx*dx
            self.Sxy = lam*self.Sxy+dx*dC
            self.Sxn = lam*self.Sxn+dx*dx*v
            self.beta = float(np.clip((self.Sxy-math.sqrt(24*self.Sxn))/(self.Sxx+2e-5),0.,1.))
        return k_c,beta_used


def job(task):
    arm,index,output = task
    started = time.perf_counter()
    case = CC.cases()[index]
    O.setup()
    sc = SE.MC.scene_for(tuple(case['cell']),CC.BASE,case['variant'])
    pol = B.build('Huber-D4',O.chosen()['Huber-D4'],sc)
    CC.own_mechanism(pol)
    recorder = ReturnedEstimate(pol.estimator.estimate,arm,index)
    pol.estimator.estimate = recorder
    log = P3.episode(sc,pol,case['seed'])
    metric = M.metric(log,pol,sc)
    if log is None:
        raise RuntimeError(f'Episode returned no log: {index} {arm}')
    dest = Path(output)/f'case_{index:02d}_{arm}.npz'
    arrays = dict(t=log[:,0],fax=log[:,18],k_eff=log[:,26],y=log[:,2],
                  t_est=np.asarray(recorder.t),k_ret=np.asarray(recorder.k),
                  n_est=np.asarray(recorder.n),trace=np.asarray(pol.trace))
    if arm == 'traj':
        arrays.update(k_w=np.asarray(recorder.kw),g=np.asarray(recorder.g))
    elif arm in ('delag', 'delag2'):
        arrays.update(k_w=np.asarray(recorder.kw),beta=np.asarray(recorder.beta),
                      force_ratio=np.asarray(recorder.force_ratio))
    with dest.open('xb') as stream:
        np.savez_compressed(stream,**arrays)
    row = dict(index=index,group=case['group'],arm=arm,seed=case['seed'],FD=float(sc['FD']),
                L=float(np.log1p(metric['J'])),**{k:float(metric[k]) for k in METRICS},
                J=float(metric['J']),event_info=sc['event_info'],
                filter_resets=recorder.filter.resets,n_est=len(recorder.t),
                trace_shape=list(arrays['trace'].shape),file=dest.name,runtime_s=time.perf_counter()-started)
    if arm == 'delag2':
        row['beta_zero_identity'] = beta_zero_identity(recorder)
    return row


def beta_zero_identity(recorder):
    zero = np.asarray(recorder.beta) == 0.
    returned = np.asarray(recorder.k, dtype=np.float64)[zero].view(np.uint64)
    window = np.asarray(recorder.kw, dtype=np.float64)[zero].view(np.uint64)
    mismatches = int(np.count_nonzero(returned != window))
    assert mismatches == 0, 'delag2 beta=0 must return k_w bit for bit'
    return dict(n=int(zero.sum()), mismatches=mismatches, bitwise_equal=True)


def subprocess_job(task, spool):
    """Independent process per job; eight concurrent slots, file IPC only."""
    arm,index,output = task
    result = spool/f'{index:02d}_{arm}.json'
    logfile = spool/f'{index:02d}_{arm}.stdout'
    command = [sys.executable,'-X','utf8',str(Path(__file__).resolve()),
               '--project-root',str(PROJECT),'--worker-arm',arm,
               '--worker-index',str(index),'--worker-output',output,
               '--worker-result',str(result)]
    with logfile.open('x',encoding='utf-8') as stream:
        done = subprocess.run(command,stdin=subprocess.DEVNULL,stdout=stream,stderr=stream,
                              creationflags=subprocess.CREATE_NO_WINDOW if os.name == 'nt' else 0)
    if done.returncode:
        raise RuntimeError(f'Worker {index} {arm} failed: '+logfile.read_text(encoding='utf-8'))
    return json.loads(result.read_text(encoding='utf-8'))


def scored_arrays(row,output):
    with np.load(output/row['file'],allow_pickle=False) as f:
        t,fax,kt = f['t'],f['fax'],f['k_eff']
        te,kr = f['t_est'],f['k_ret']
    proxy = np.zeros((len(t),10))
    proxy[:,0],proxy[:,7],proxy[:,9] = t,fax,kt
    _,masks = D.masks_for(proxy,row)
    ids = np.rint(te/V3.TC).astype(int)
    assert np.all((ids >= 0)&(ids < len(t)))
    assert np.max(np.abs(t[ids]-te),initial=0.) < 1e-9
    truth = np.interp(te,t,kt)
    at_force = np.interp(te,t,fax)
    score = (te > 1.5)&(truth > 0)&(at_force > .5)
    e = np.full(len(te),np.nan)
    finite = np.isfinite(kr)&(kr > 0)&(truth > 0)
    e[finite] = np.log(kr[finite]/truth[finite])
    segments = {s:masks[s][ids]&score for s in SEGMENTS}
    # Missing or skipped calls must break sustained windows, never bridge gaps.
    full_e = np.full(len(t),np.nan)
    full_e[ids] = e
    sustained = D.sustained(full_e,masks['outside_event'])
    force = {s:fax[masks[s]]-row['FD'] for s in ('fast','plateau')}
    return e,segments,sustained,force


def error_stats(e,mask):
    good = mask&np.isfinite(e)
    return dict(D.stats(e[good]),n_scored=int(mask.sum()),n_missing=int((mask&~good).sum()))


def summarize(rows,output):
    cache = {}
    for row in rows:
        e,masks,sustained,force = scored_arrays(row,output)
        cache[row['arm'],row['index']] = e,masks,sustained,force
        row['identification'] = {s:error_stats(e,masks[s]) for s in SEGMENTS}
        row['sustained_underestimate'] = sustained
        row['force_rms'] = {s:dict(n=len(a),rms=D.rms(a) if len(a) else None) for s,a in force.items()}
    def aggregate(subset,identification=False):
        result = {}
        for arm in ('huber','traj','traj_priv','delag','delag2'):
            rr = [r for r in subset if r['arm'] == arm]
            if not rr:
                continue
            data = [cache[arm,r['index']] for r in rr]
            rec = dict(cases=[r['index'] for r in rr],mean_L=float(np.mean([r['L'] for r in rr])),
                       median_L=float(np.median([r['L'] for r in rr])),force_rms={})
            for s in ('fast','plateau'):
                a = np.concatenate([d[3][s] for d in data])
                rec['force_rms'][s] = dict(n=len(a),rms=D.rms(a) if len(a) else None)
            if identification:
                e = np.concatenate([d[0] for d in data])
                rec['segments'] = {s:error_stats(e,np.concatenate([d[1][s] for d in data])) for s in SEGMENTS}
                vals = [d[2] for d in data if d[2] is not None]
                rec['sustained'] = dict(median=float(np.median(vals)) if vals else None,
                                        worst=max(vals) if vals else None)
            result[arm] = rec
        return result
    groups = {g:aggregate([r for r in rows if r['group'] == g]) for g in D.BASE_GROUPS}
    subsets = [('A',tuple(range(6))),('A_slow_thin',tuple(range(6,12))),
               ('B',tuple(range(12,18))),('A+B',tuple(range(18,24))),
               ('A_large_fast',FAST),('A_large_slow',(38,41))]
    ident = {g:aggregate([r for r in rows if r['index'] in ii],True) for g,ii in subsets}
    return groups,ident


def predictions(rows,groups,identification):
    fast = identification['A_large_fast']
    h,c = fast['huber']['segments'],fast['traj']['segments']
    reductions = {s:1-c[s]['median_abs']/h[s]['median_abs'] for s in ('fall','reversal')}
    p1 = dict(met=bool(abs(c['rise']['mean']) <= .06 and c['fall']['mean'] < .12
                      and all(v >= .25 for v in reductions.values())),
              rise_mean={a:fast[a]['segments']['rise']['mean'] for a in ('huber','traj')},
              fall_mean={a:fast[a]['segments']['fall']['mean'] for a in ('huber','traj')},
              median_abs_reduction=reductions,
              per_case={str(i):{a:next(r for r in rows if r['index'] == i and r['arm'] == a)['identification']
                               for a in ('huber','traj')} for i in FAST})
    get = {(r['arm'],r['index']):r for r in rows}
    large = groups['A_large']
    deltas = {str(i):get['traj',i]['L']-get['huber',i]['L'] for i in range(36,42)}
    p2 = dict(met=bool(large['traj']['mean_L'] < large['huber']['mean_L']
                      and large['traj']['mean_L'] <= .22 and max(deltas.values()) <= .05),
              huber_mean=large['huber']['mean_L'],traj_mean=large['traj']['mean_L'],
              per_case_delta=deltas,worst_delta=max(deltas.values()))
    normal = {g:groups[g]['traj']['mean_L']-groups[g]['huber']['mean_L'] for g in ('A','A_slow_thin')}
    p3 = dict(met=bool(all(abs(v) <= .02 for v in normal.values())),mean_L_delta=normal)
    weak = {g:groups[g]['traj']['mean_L']-groups[g]['huber']['mean_L'] for g in ('B','A+B')}
    p4 = dict(met=bool(all(v > 0 for v in weak.values())),mean_L_delta=weak)
    priv = {str(i):get['traj_priv',i]['L']-get['traj',i]['L'] for i in GAINS}
    p5 = dict(met=bool(all(v <= 0 for v in priv.values())),compared_cases=list(GAINS),
              privileged_minus_learned_L=priv,
              paired_mean_L={a:float(np.mean([get[a,i]['L'] for i in GAINS])) for a in ('traj','traj_priv')},
              interpretation='At least as good: no larger L on each prescribed paired case; no privileged run on 41.')
    return dict(P1=p1,P2=p2,P3=p3,P4=p4,P5=p5)


def round2_predictions(rows,groups,identification,beta_stats):
    get = {(r['arm'],r['index']):r for r in rows}
    delta = {i:get['delag',i]['L']-get['huber',i]['L'] for i in CASES}
    large = groups['A_large']['delag']['mean_L']
    q1 = dict(met=bool(large <= .22 and max(delta[i] for i in range(36,42)) <= .02),
              mean_L=large,per_case_delta={str(i):delta[i] for i in range(36,42)},
              worst_delta=max(delta[i] for i in range(36,42)))
    normal = {g:groups[g]['delag']['mean_L']-groups[g]['huber']['mean_L']
              for g in ('A','A_slow_thin','B','A+B')}
    mean_ok = {g:abs(d) <= .02 for g,d in normal.items()}
    exception = normal['A_slow_thin'] < -.02 and delta[7] < 0
    q2 = dict(met=bool(max(delta[i] for i in range(24)) <= .05
                       and all(ok or (g == 'A_slow_thin' and exception) for g,ok in mean_ok.items())),
              mean_L_delta=normal,within_002=mean_ok,slow_thin_improvement_exception=bool(exception),
              case7_delta=delta[7],worst_delta=max(delta[i] for i in range(24)),
              per_case_delta={str(i):delta[i] for i in range(24)})
    fast = identification['A_large_fast']['delag']['segments']
    q3 = dict(met=bool(abs(fast['rise']['mean']) <= .06 and fast['fall']['mean'] < .12),
              rise_mean=fast['rise']['mean'],fall_mean=fast['fall']['mean'])
    fractions = {str(i):beta_stats[str(i)]['fraction_above_005'] for i in FAST}
    b_medians = {str(i):beta_stats[str(i)]['median'] for i in range(12,18)}
    q4 = dict(met=bool(all(v > .5 for v in fractions.values()) and all(v < .1 for v in b_medians.values())),
              large_fast_fraction_above_005=fractions,B_median=b_medians,
              convention='For each case, estimator ticks t > 3; most means strictly over half; beta is the gain used.')
    return dict(Q1=q1,Q2=q2,Q3=q3,Q4=q4)


def run_round2(project,output_root):
    base = output_root/'experiments/rossi_failure_study/results'
    output,dest = base/'trajectory_loop',base/'trajectory_loop.json'
    source = project/'experiments/rossi_failure_study/results/trajectory_loop.json'
    result = json.loads(source.read_text(encoding='utf-8'))
    if any(r['arm'] == 'delag' for r in result['rows']):
        raise FileExistsError('Round 2 rows already exist')
    output.mkdir(parents=True,exist_ok=True)
    if any((output/f'case_{i:02d}_delag.npz').exists() for i in CASES):
        raise FileExistsError('Refusing to overwrite delag traces')
    # In a mirror, preserve all round-1 traces for the shared scoring routine.
    if output_root != project:
        import shutil
        for row in result['rows']:
            target = output/row['file']
            if not target.exists():
                shutil.copy2(source.parent/'trajectory_loop'/row['file'],target)
    specification = dict(
        Q1='A_large mean L <= .22; every delag-huber case delta <= .02.',
        Q2='Cases 0-23 delta <= .05; group mean delta within +/-.02, allowing A_slow_thin improvement through case 7.',
        Q3='Pooled A_large fast rise mean log error within +/-.06; fall mean < +.12.',
        Q4='Each A_large fast case: beta > .05 for >50% of t>3 ticks; each B case median beta < .1.')
    print('Round 2 predictions BEFORE RUN '+json.dumps(specification),flush=True)
    started = time.perf_counter()
    spool = Path(tempfile.mkdtemp(prefix='trajectory_loop_round2_workers_'))
    rows = []
    with ThreadPoolExecutor(8) as pool:
        pending = [pool.submit(subprocess_job,('delag',i,str(output)),spool) for i in CASES]
        for future in as_completed(pending):
            row = future.result()
            rows.append(row)
            print(f"delag {row['index']:02d}: L={row['L']:.9f}, peak={row['peak']:.4f}",flush=True)
    # Compute shared summaries without changing any stored row's fields.
    import copy
    groups,identification = summarize(copy.deepcopy(result['rows'])+rows,output)
    beta_stats = {}
    for row in rows:
        with np.load(output/row['file'],allow_pickle=False) as f:
            b = f['beta'][f['t_est'] > 3.]
        if not len(b) or not np.all(np.isfinite(b)):
            raise ValueError('Missing beta ticks after 3 s')
        beta_stats[str(row['index'])] = dict(n=len(b),median=float(np.median(b)),
            final=float(b[-1]),fraction_above_005=float(np.mean(b > .05)))
    for comparison in result['comparison']:
        row = next(r for r in rows if r['index'] == comparison['index'])
        comparison['L']['delag'],comparison['peak']['delag'] = row['L'],row['peak']
    result['rows'].extend(rows)
    result['rows'].sort(key=lambda r:(r['index'],r['arm']))
    result['groups'],result['identification'] = groups,identification
    result['round2'] = dict(status='complete',workers=8,cases=list(CASES),new_episodes=len(rows),
        reused_rows=len(result['rows'])-len(rows),runtime_s=time.perf_counter()-started,
        prespecified_predictions=specification,
        parameters=dict(short_window=6,long_window=24,gain_lag=24,Tc=V3.TC,
            decay_time=4.,noise_overlap=24,ridge=2e-5,r_bounds=[.5,2.],k_bounds=[200.,3e5]),
        beta=beta_stats,predictions=round2_predictions(result['rows'],groups,identification,beta_stats),
        conventions=dict(beta='Gain used before current tick update; all ticks included in t>3 statistics.',
            history='Every estimator call, including invalid windows; nonconsecutive n resets all delag histories and gain sums.',
            contact='Gain update requires F>=1 at all 25 ticks j-24 through j, inclusive.',
            decay='Sums decay only on eligible gain updates, as prescribed.',
            scoring='Unchanged round-1 masks, returned-stiffness log error, and pooled force RMS.'))
    # Verify the stored rows were not altered by shared scoring.
    original = json.loads(source.read_text(encoding='utf-8'))
    assert [r for r in result['rows'] if r['arm'] != 'delag'] == original['rows']
    with dest.open('w',encoding='utf-8') as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write('\n')
    print('Round 2 verdicts '+json.dumps(result['round2']['predictions']),flush=True)
    print(f"Written {dest}; 30 new episodes; elapsed {result['round2']['runtime_s']:.1f}s",flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root',type=Path,default=PROJECT)
    parser.add_argument('--output-root',type=Path)
    parser.add_argument('--round2',action='store_true')
    parser.add_argument('--worker-arm',choices=('huber','traj','traj_priv','delag','delag2'))
    parser.add_argument('--worker-index',type=int)
    parser.add_argument('--worker-output',type=Path)
    parser.add_argument('--worker-result',type=Path)
    args = parser.parse_args()
    if args.worker_arm:
        row = job((args.worker_arm,args.worker_index,str(args.worker_output)))
        with args.worker_result.open('x',encoding='utf-8') as stream:
            json.dump(row,stream,allow_nan=False)
        return
    project = args.project_root.resolve()
    output_root = (args.output_root or project).resolve()
    if args.round2:
        run_round2(project,output_root)
        return
    base = output_root/'experiments/rossi_failure_study/results'
    output,dest = base/'trajectory_loop',base/'trajectory_loop.json'
    if (output.exists() and any(output.iterdir())) or dest.exists():
        raise FileExistsError('Refusing to overwrite trajectory_loop outputs')
    output.mkdir(parents=True,exist_ok=True)
    started = time.perf_counter()
    checks = equivalence(project)
    refs = json.loads((project/'experiments/rossi_failure_study/results/completion_closed_loop.json').read_text(encoding='utf-8'))['rows']
    stored = {(r['method'],r['index']):r for r in refs}
    rows,regression = [],[]
    spool = Path(tempfile.mkdtemp(prefix='trajectory_loop_workers_'))
    with ThreadPoolExecutor(8) as pool:
        pending = [pool.submit(subprocess_job,('huber',i,str(output)),spool) for i in CASES]
        for future in as_completed(pending):
            row = future.result()
            delta = abs(row['L']-stored['Huber-D4_nogate',row['index']]['L'])
            regression.append(dict(index=row['index'],abs_dL=delta))
            rows.append(row)
            print(f"huber {row['index']:02d}: L={row['L']:.9f}, |dL|={delta}",flush=True)
        worst = max(r['abs_dL'] for r in regression)
        if worst != 0.:
            with dest.open('x',encoding='utf-8') as stream:
                json.dump(dict(status='STOPPED: harness regression failed',equivalence=checks,
                               regression=regression,max_abs_dL=worst,rows=rows),stream,indent=2)
            raise RuntimeError(f'Nonzero Huber regression {worst}; no candidate loops run')
        print('Harness passed: max |dL| = 0; starting fixed candidate runs',flush=True)
        tasks = [('traj',i,str(output)) for i in CASES]+[('traj_priv',i,str(output)) for i in GAINS]
        pending = [pool.submit(subprocess_job,task,spool) for task in tasks]
        for future in as_completed(pending):
            row = future.result()
            rows.append(row)
            print(f"{row['arm']} {row['index']:02d}: L={row['L']:.9f}, peak={row['peak']:.4f}",flush=True)
    rows.sort(key=lambda r:(r['index'],r['arm']))
    groups,identification = summarize(rows,output)
    comparison = []
    for index in CASES:
        rr = [r for r in rows if r['index'] == index]
        comparison.append(dict(index=index,group=rr[0]['group'],L={r['arm']:r['L'] for r in rr},
                               peak={r['arm']:r['peak'] for r in rr},
                               stored={name:{k:stored[name,index][k] for k in ('L','peak')}
                                       for name in ('Huber-D4_imm','oracle_keff_off')}))
    result = dict(status='complete',workers=8,execution_backend='Eight concurrent independent subprocess workers; file IPC (Windows denied multiprocessing pipe)',cases=list(CASES),
                  filter_parameters=dict(qs=1.,qo=1.,qbeta=1e-7,f0=.25),privileged_gains=GAINS,
                  equivalence=checks,harness_regression=dict(max_abs_dL=worst,per_case=regression),
                  rows=rows,comparison=comparison,groups=groups,identification=identification,
                  predictions=predictions(rows,groups,identification),runtime_s=time.perf_counter()-started,
                  conventions=dict(group_errors='Pooled finite estimator returns; missing outputs counted, never held or imputed.',
                      segments='D.masks_for on each own full control-grid k_eff and fax, sampled at estimator ticks; events separate.',
                      force_rms='Pooled full control-grid force samples using same masks; fast=rise|fall|reversal.',
                      sustained='D.sustained: 13 contiguous control ticks outside events; missing/skipped calls break windows.',
                      privileged_startup='Retain raw before 12 filter ticks; then prescribed correction without additional estimator clipping; controller retains its own gain bounds.',
                      filter_history='Estimator-call history; discontinuous n resets filter and trailing histories.',
                      privileged_scope='36-40 only: aggregate excludes 41; slow privileged identification contains 38 only.'))
    with dest.open('x',encoding='utf-8') as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write('\n')
    print('Predictions '+json.dumps(result['predictions']),flush=True)
    print(f"Written {dest}; {len(rows)} episodes; elapsed {result['runtime_s']:.1f}s",flush=True)


if __name__ == '__main__':
    main()

