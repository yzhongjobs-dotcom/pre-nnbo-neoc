"""Compile the SAME prediction derivatives for a fair implementation-cost check.

Generated C/binaries live outside the repository. This is a timing diagnostic,
not a new controller. No comparison to embedded solvers without a matched run.
"""
import argparse,json,subprocess,time,os
from pathlib import Path
import numpy as np
import casadi as ca
from optimization import Problem,context,FastNeighbor,ReusedOptimizer
from core import m,plant,local_contact,encode
from experiment import ROOT,BASE


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--path',default='/tmp/impedance_ocp_codegen');ap.add_argument('--rates',action='store_true');ap.add_argument('--true-plant',action='store_true');ap.add_argument('--bisection',action='store_true');args=ap.parse_args();out=Path(args.path);out.mkdir(exist_ok=True)
    pr=Problem(rate_actions=args.rates,gains_only=args.rates,true_plant=args.true_plant,dt=.002 if args.rates else .004,closed_projection=not args.bisection);old=os.getcwd();os.chdir(out)
    try:
        cg=ca.CodeGenerator('prediction.c')
        # Full Hessian/cross generated ~15 MB C and exceeded the current
        # compiler memory limit at -O3. Compile the shared online value/gradient
        # kernel; keep occasional curvature preparation as graph evaluation.
        for f in (pr.value,pr.vg):cg.add(f)
        cg.generate();tic=time.perf_counter();subprocess.run(['gcc','-O3','-fPIC','-shared','prediction.c','-o','prediction.so','-lm'],check=True);compile_s=time.perf_counter()-tic
    finally:os.chdir(old)
    cfg=np.array(json.loads((BASE/'tuned_fixed.json').read_text())['cfg']);p=m.sample(1,9011)[0];cfg[7]=p[9];y=plant.initial(p,cfg[:4],65);zr=y[1]
    q=np.array([np.log(p[0]),0,0,0,0,np.log(p[12]),p[10]]);anchor=np.zeros(2);sensed=local_contact(0.,y,q,anchor,p[9])[0]
    c=context(0,y,zr,p if args.true_plant else q,anchor,sensed,cfg);a=np.zeros(10) if args.rates else np.tile(encode(np.r_[cfg[:4],0]),2);rows={};difference={}
    for name in ('value','vg'):
        f=getattr(pr,name);native=ca.external(f.name(),str(out/'prediction.so'));v1=f(a,c);v2=native(a,c)
        if not isinstance(v1,tuple):v1=(v1,);v2=(v2,)
        difference[name]=max(float(np.max(abs(np.array(x)-np.array(z)))) for x,z in zip(v1,v2))
        assert difference[name]<1e-8,(name,difference[name])
        times=[]
        for ff in (f,native):
            tic=time.perf_counter()
            for _ in range(100):ff(a,c)
            times.append((time.perf_counter()-tic)*10.)
        rows[name]={'graph_ms':times[0],'compiled_ms':times[1]};setattr(pr,name,native)
    record={'compile_seconds':compile_s,'timings':rows,'max_differences':difference,'path':str(out),'casadi':ca.__version__,'compiled':['value','vg'],'full_hessian_compile':'initial -O3 attempt exceeded compiler memory; not a runtime result'}
    name='compiled_training_cost_check.json' if args.true_plant else ('compiled_rate_cost_check.json' if args.rates else 'compiled_cost_check.json')
    (ROOT/name).write_text(json.dumps(record,indent=2));print(record,flush=True)
    (out/'manifest.json').write_text(json.dumps({'rate_actions':args.rates,'true_plant':args.true_plant,'prediction_version':'sampled_force_v2' if args.rates else 'target_legacy','dt':pr.dt,'projection':'bisection' if args.bisection else 'closed_when_rate'},indent=2))


if __name__=='__main__':main()
