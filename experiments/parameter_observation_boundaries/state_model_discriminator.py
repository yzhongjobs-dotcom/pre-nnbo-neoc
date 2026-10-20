"""Original Hertz-footprint field: physical variation vs information revision.

The scalar Bayesian coefficient audit uses known geometry and boundary plus an
ideal normal first moment. It is NOT the full MHE or raw wrist-moment pipeline.
"""
from pathlib import Path
import json
import numpy as np

ROOT=Path(__file__).resolve().parent
A=.003;R=.03;L=.02;E0=30000.;EF=45000.;B=0.;V=.0006;H=.05
CH=L/(2*R);I0=np.pi*A*A/2


def features(s):
    r=np.clip(B-s,-A,A);g=np.sqrt(max(A*A-r*r,0.))
    If=.5*(A*A*np.arccos(r/A)-r*g)
    I1=(A*A-r*r)**1.5/3
    return np.array([CH*If,CH*I1])


def posterior(path,with_moment,seed=43):
    rng=np.random.default_rng(seed);noise=np.array([.003,2e-6]);prior_variance=10000.**2
    variance=prior_variance;mean=E0;rows=[]
    for t,s,eff in path:
        d=features(s);measurement=d*(EF-E0)+rng.normal(size=2)*noise
        idx=[0,1] if with_moment else [0]
        information=np.sum((d[idx]/noise[idx])**2)
        numerator=(mean-E0)/variance+np.sum(d[idx]*measurement[idx]/noise[idx]**2)
        variance=1/(1/variance+information);mean=E0+variance*numerator
        rows.append([t,s,eff,mean,np.sqrt(variance),*measurement])
    return np.array(rows)


def main():
    out=ROOT/'state_model_discriminator';out.mkdir(exist_ok=True)
    times=np.arange(0,A/V+H/2,H);positions=-A+V*times
    effective=np.array([E0+(EF-E0)*features(s)[0]/(CH*I0) for s in positions])
    path=np.column_stack((times,positions,effective))
    curves={};seed_records=[]
    for mode,with_moment in [('F_history',False),('F_plus_ideal_moment',True)]:
        a=posterior(path,with_moment);curves[mode]=a
        np.savetxt(out/(mode+'.csv'),a,delimiter=',',
          header='t,s,effective_E,estimated_Ef,posterior_std_Ef,delta_F,ideal_m1',comments='')
        for seed in range(43,63):
            z=posterior(path,with_moment,seed)
            seed_records.append(dict(mode=mode,seed=seed,
                maximum_descriptor_increment=float(np.max(abs(np.diff(z[:,3])))),
                maximum_physical_increment=float(np.max(abs(np.diff(z[:,2])))),
                descriptor_relative_increment=float(np.max(abs(np.diff(z[:,3])))/E0),
                final_error=float(z[-1,3]-EF)))
    alpha=features(positions[1])[0]/(CH*I0)
    # Exact noiseless observation identifies the endpoint coefficient whenever
    # alpha>0. With a uniform-model prior, this first revision is 50% even though
    # the presently encountered effective material changes only ~0.03%.
    recovered=E0+((E0+(EF-E0)*alpha)-E0)/alpha
    result=dict(model=dict(R=R,L=L,a=A,b=B,scan_velocity=V,period=H,current_E=E0,true_front_E=EF,
                    noise_F=.003,noise_ideal_moment=2e-6,prior_std_front_E=10000.,known_geometry_and_boundary=True),
       physical_path=dict(samples=len(path),maximum_relative_physical_increment=float(np.max(abs(np.diff(effective))))/E0,
                  first_relative_physical_increment=float((effective[1]-effective[0])/E0)),
       exact_information_counterexample=dict(first_overlap_fraction=alpha,recovered_front_E=recovered,
           first_relative_descriptor_revision=(recovered-E0)/E0,
           conclusion='Slow current effective stiffness does not imply small revisions of the unknown future material descriptor.'),
       correct_spatial_model=dict(descriptor=[E0,EF,B],descriptor_revisions_due_to_scanning_alone=0,
                  conclusion='A known fixed spatial contact field changes its current output with scan position; that state change alone does not require a new HJB root.'),
       noisy_information_audit=seed_records,
       scope='Model-internal identifiability and information-update counterexample. Known boundary and ideal first moment; not full MHE, raw wrench validation or control performance.')
    (out/'results.json').write_text(json.dumps(result,indent=2),encoding='utf8')
    print(json.dumps({k:v for k,v in result.items() if k!='noisy_information_audit'},indent=2),flush=True)
    assert abs(recovered-EF)<1e-5 and result['physical_path']['maximum_relative_physical_increment']<.004


if __name__=='__main__':main()
