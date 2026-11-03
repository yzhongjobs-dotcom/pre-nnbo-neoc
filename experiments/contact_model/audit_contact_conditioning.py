"""Sensitivity of profiled contact fitting on IDENTICAL frozen measurements.

No reference stiffness used. Compares machine-scale outer differences with
central differences at a declared scale; the scale sweep diagnoses derivatives,
not a search for better task RMSE. Closed-loop claims require a separate run.
"""
import json
import argparse
from pathlib import Path
import numpy as np
import dynamic_contact as DC
from audit_contact_roots import old_inverse

HERE=Path(__file__).resolve().parent


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--stable-only',action='store_true');args=ap.parse_args()
    a=np.load(HERE/'results'/'finite_patch_input_S2_110000.npz')
    signals=DC.Signals();windows=[]
    for j,row in enumerate(a['measured'][:128]):
        signals.tau_servo=float(a['baseline_tau'][j]);signals.append(row)
        if j+1 in (64,80,96,112,128):
            windows.append((signals.window(),float(a['baseline_kraw'][j])))
    original=DC.least_squares;rows=[]
    for derivative in (('stable_absolute',) if args.stable_only else ('machine_forward','central_1e-4')):
        def solve(*args,**kw):
            if derivative=='central_1e-4':kw.update(jac='3-point',diff_step=1e-4)
            return original(*args,**kw)
        DC.least_squares=solve
        for roots in ('fixed','early'):
            fitter=DC.ContactFit(spatial=True,stable_jacobian=derivative=='stable_absolute')
            if roots=='fixed':
                fitter.model.inverse=lambda x,f,t:old_inverse(fitter.model,x,f,t)
            for obs,hint in windows:
                fit=fitter.fit(obs,hint)
                theta=fit['theta'];T=fitter.memory.operator(16)
                def residual(th):return fitter.surface_profile(th,obs,T)[0]
                scale=np.array([1.,1.,.01,.01]);js=[]
                for h in (1e-2,1e-3,1e-4,1e-5,1e-6,1e-7):
                    J=np.column_stack([(residual(theta+np.eye(4)[i]*scale*h)-
                                        residual(theta-np.eye(4)[i]*scale*h))/(2*h)
                                       for i in range(4)])
                    js.append(J)
                ref=js[2];norm=max(float(np.linalg.norm(ref)),1e-12)
                rows.append(dict(derivative=derivative,roots=roots,time=obs['tend'],
                                 k=fit['k'],risk=fit['risk'],success=fit['success'],theta=theta.tolist(),
                                 depth=fit['depth'],residual_rms=fit['residual_rms'],
                                 jacobian_step_sweep_relative=[float(np.linalg.norm(J-ref)/norm) for J in js]))
        print(derivative,'done',flush=True)
    DC.least_squares=original
    result=dict(protocol=dict(input='finite_patch_input_S2_110000.npz',truth_used=False,
                              steps=[1e-2,1e-3,1e-4,1e-5,1e-6,1e-7]),rows=rows)
    label='contact_conditioning_stable_audit.json' if args.stable_only else 'contact_conditioning_audit.json'
    (HERE/'results'/label).write_text(json.dumps(result,indent=2),encoding='utf8')
    for r in rows:
        print(r['derivative'],r['roots'],round(r['time'],3),round(r['k'],2),round(r['risk'],3),
              np.round(r['jacobian_step_sweep_relative'],5),flush=True)


if __name__=='__main__':main()
