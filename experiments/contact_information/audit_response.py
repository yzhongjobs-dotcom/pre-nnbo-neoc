"""Same-input response prediction and control-functional information audit."""
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import memory_contact as MC
import structured_memory as SM
import local_contact as LC
from scipy.optimize import lsq_linear
from dynamic_compliance import compliance_information

HERE=Path(__file__).resolve().parent


def matched_free(measured,tau):
    """Same FOUR nuisance parameters as shared CA: B0,B1,B2,initial exp.

This checks whether shared-physics gains are just fewer nuisance parameters.
It is a matched-capacity diagnostic, not a new proposed contact model.
"""
    r=MC.fit(measured,model='memoryc_040',tau=tau)
    if r is None:return None
    n=len(measured);j=np.arange(n-33,n-1)
    B,U,_=MC.features(measured[j,0],measured[j,2],measured[j+1,2],r['te'],r['ye'],r['tw'],
                      tau,r['tau_m'],r['origin'])
    C=B[:,[0,1,2,5]];X=LC.residualize(C,U)
    z=measured[j+1,3]-math.exp(-2*math.pi*30*.008)*measured[j,3]
    zr=LC.residualize(C,z)
    transform=np.array([[0.,1.,0.],[-1.,1.,0.],[0.,0.,1.]])
    A=X@transform;scale=np.linalg.norm(A,axis=0)
    solution=lsq_linear(A/scale,zr,bounds=(np.array([200,200,0])*scale,
                                                    np.array([3e5,3e5,3e5])*scale),tol=1e-12)
    theta=transform@(solution.x/scale);e=zr-X@theta;dof=len(j)-7
    chi=float(e@e/(dof*.04**2));norm=np.linalg.norm(X,axis=0)
    _,s,v=np.linalg.svd(X/norm,full_matrices=False)
    cov=.04**2*max(1.,chi)*((v.T/s**2)@v)/norm[:,None]/norm[None,:]
    beta=np.zeros(6);beta[[0,1,2,5]]=np.linalg.lstsq(C,z-U@theta,rcond=1e-12)[0]
    omega=2*math.pi*20.8*r['tau_m'];functional=np.array([1.,0.,omega**2/(1+omega**2)])
    gain=float(functional@theta);sd=math.sqrt(float(functional@cov@functional))
    r.update(coefficients=theta,cov=cov,beta=beta,chi2=chi,dof=dof,k=gain,sigma=sd,
             risk=sd/gain*(1+math.sqrt(dof*max(chi-1,0))),physical=bool(solution.success),
             min_gain=float(min(theta[0],theta[0]-theta[1])))
    return r


def predict(r,measured,n,h,tau,memory):
    j=np.arange(n-1,n-1+h)
    if memory:
        B,U,_=MC.features(measured[j,0],measured[j,2],measured[j+1,2],r['te'],r['ye'],r['tw'],
                          tau,r['tau_m'],r['origin'],r['qend'])
    else:
        B,U=LC.features(measured[j,0],measured[j,2],measured[j+1,2],r['te'],r['ye'],r['tw'],tau)
        U=U[:,:2]
    decay=math.exp(-2*math.pi*30*.008);w=decay**np.arange(h-1,-1,-1)
    return float(decay**h*measured[n-1,3]+w@(B@r['beta']+U@r['coefficients']))


def main():
    target=HERE/'results/shared_response_same_input_matched_v2.json'
    if target.exists():raise ValueError('preserve output')
    result=dict(protocol=dict(conditioning='future encoder supplied, future force hidden; not online surface prediction',
        stride=16,windows_independent=False,input_hashes={},no_truth=True),rows=[],summary=[])
    for family in ('T3m','T2m-stiff'):
        path=next((HERE/'results').glob(f'lean_input_{family}_*.npz'))
        result['protocol']['input_hashes'][path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        d=np.load(path);m=d['measured'];taus=d['causal_tau']
        for n in range(192,len(m)-4,16):
            tau=float(taus[n-1])
            for model in ('lean','free','free_matched','shared'):
                r=(LC.fit(m[:n],nw=32,rate=False,return_aux=True,tau=tau) if model=='lean' else
                   MC.fit(m[:n],model='memoryc_040',tau=tau) if model=='free' else
                   matched_free(m[:n],tau) if model=='free_matched' else SM.fit(m[:n],tau=tau))
                if r is None:continue
                row=dict(family=family,model=model,end=n,chi2=r['chi2'],k=r['k'],probe_risk=r['risk'])
                if model!='lean':row['inverse_risk']=compliance_information(r,.031)['risk']
                for h in (1,4):row[f'h{h}_error']=predict(r,m,n,h,tau,model!='lean')-m[n-1+h,3]
                result['rows'].append(row)
        for model in ('lean','free','free_matched','shared'):
            rows=[r for r in result['rows'] if r['family']==family and r['model']==model]
            summary=dict(family=family,model=model,windows=len(rows),
                rmse={str(h):float(np.sqrt(np.mean([r[f'h{h}_error']**2 for r in rows]))) for h in (1,4)},
                median_probe_risk=float(np.median([r['probe_risk'] for r in rows])))
            if model!='lean':summary['median_inverse_risk']=float(np.median([r['inverse_risk'] for r in rows]))
            result['summary'].append(summary);print(summary)
    target.write_text(json.dumps(result,indent=2),encoding='utf8')


if __name__=='__main__':main()
