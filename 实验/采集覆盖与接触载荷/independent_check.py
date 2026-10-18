"""Unilateral-contact equilibrium and mesh audit of the envelope screen."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.optimize import minimize
import envelope as e

HERE = Path(__file__).resolve().parent


def unilateral(k,d):
    # min 1/2 w'Kw subject w>=d. Here w is downward tissue displacement;
    # separation is w-d, forces Kw>=0, and their product must vanish.
    x=d/e.MM
    out=minimize(lambda v: (.5*v@k@v,k@v),x,jac=True,
        bounds=[(v,None) for v in x],method='L-BFGS-B',
        options=dict(ftol=1e-15,gtol=1e-10,maxiter=5000,maxls=100))
    w=out.x*e.MM;f=k@w;gap=w-d
    residual=max(float(max(0.,-f.min())),float(np.max(np.abs(f*(gap>1e-8)))))
    return dict(contact_node_fraction=float(np.mean(gap<=1e-8)),
        max_separation_m=float(gap.max()),force_N=float(f.sum()),
        force_kkt_residual_N=residual,
        minimum_gap_m=float(gap.min()),success=bool(out.success))


def main():
    raw=np.load(HERE/'paired.npz')
    rng=np.random.default_rng(7301);records=[]
    for ell in [.0005,.0015,.003]:
        for case,p in enumerate(raw['coupled_parameters']):
            for s in np.linspace(0,.08,9):
                _,_,qq=e.evaluate(p,s,rng,shear_length=ell)
                u,g,k=e.system(p,s,shear_length=ell)
                d=e.MM*(qq[1,0]+qq[1,1]*u/e.A)-g
                records.append(dict(ell=ell,case=case,s=float(s),**unilateral(k,d)))
    groups=[]
    for ell in [.0005,.0015,.003]:
        rr=[r for r in records if r['ell']==ell]
        groups.append(dict(ell=ell,cases=len(rr),
            mean_actual_contact_node_fraction=float(np.mean([r['contact_node_fraction'] for r in rr])),
            minimum_actual_contact_node_fraction=min(r['contact_node_fraction'] for r in rr),
            max_separation_m=max(r['max_separation_m'] for r in rr),
            max_force_kkt_residual_N=max(r['force_kkt_residual_N'] for r in rr)))
    # Preselected first three scenes, three locations, two models. Compare
    # optimal load on increasingly dense footprints, not only root residuals.
    meshes=[]
    for case,p in enumerate(raw['coupled_parameters'][:3]):
        for s in [0.,.04,.08]:
            for ell in [0.,.003]:
                forces=[]
                for n in [61,121,241]:
                    u,g,k=e.system(p,s,n,ell);q=e.lp(u,g,k)
                    forces.append(e.metrics(q,u,g,k)[0])
                meshes.append(dict(case=case,s=s,ell=ell,forces=forces,
                    relative_change_121_to_241=abs(forces[-1]/forces[-2]-1)))
    # Independent coordinate-descent complementarity solution, one fixed
    # largest-coupling case. It does not call the optimizer above.
    p=raw['coupled_parameters'][0];u,g,k=e.system(p,.04,shear_length=.003)
    _,_,qq=e.evaluate(p,.04,rng,shear_length=.003)
    d=e.MM*(qq[1,0]+qq[1,1]*u/e.A)-g;w=d.copy()
    for it in range(20000):
        old=w.copy()
        for i in range(len(w)):
            w[i]=max(d[i],w[i]-(k[i]@w)/k[i,i])
        if np.max(np.abs(old-w))<1e-13:break
    independent=unilateral(k,d)
    comparison=dict(coordinate_descent_iterations=it+1,
        coordinate_descent_force_N=float(np.sum(k@w)),
        optimizer_force_N=independent['force_N'],
        force_difference_N=abs(float(np.sum(k@w))-independent['force_N']))
    result=dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        envelope_sha256=hashlib.sha256((HERE/'envelope.py').read_bytes()).hexdigest(),
        note='Actual unilateral equilibrium, not a count of negative forces in an assumed all-active state. Contact-node fraction is a grid measure, not acoustic image quality.',
        groups=groups,mesh_cases=meshes,independent_solver=comparison,
        max_mesh_relative_change_121_to_241=max(m['relative_change_121_to_241'] for m in meshes),
        contact_cases=records)
    (HERE/'independent_checks.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(groups=groups,mesh_max=result['max_mesh_relative_change_121_to_241'],
                         solver=comparison),indent=2))


if __name__=='__main__':main()
