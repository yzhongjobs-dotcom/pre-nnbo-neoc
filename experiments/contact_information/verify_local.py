"""Algebraic identification and signal constraints before closed-loop use."""
import json
from pathlib import Path
import numpy as np
from local_contact import features,fit,design,optimize_probe,periodic_encoder,residualize


def main():
    probe=optimize_probe(n=32,rate=True);wave=np.array(probe['wave'])
    assert probe['peak']<=1+1e-12 and probe['rms']<=1/np.sqrt(2)+1e-12
    assert probe['max_slew']<=2*np.sin(np.pi*20.8*.008)+1e-12
    y=np.tile(periodic_encoder(wave),4)*.0001+.002;t=np.arange(len(y))*.008
    rows=np.c_[t,.02*t,y,np.zeros(len(t)),np.zeros((len(t),3))]
    B,U=features(t[-33:-1],y[-33:-1],y[-32:],t[-1],y[-1],.256)
    theta=np.array([12000.,2000.,30.]);beta=np.array([4.,.1,.02,-.03,.01])
    z=U@theta+B@beta;decay=np.exp(-2*np.pi*30*.008)
    rows[-33,3]=4.
    for j,force in enumerate(z):rows[-32+j,3]=decay*rows[-33+j,3]+force
    r=fit(rows,nw=32,rate=True)
    np.testing.assert_allclose(r['coefficients'],theta,rtol=1e-9,atol=1e-6)
    constant=fit(rows,nw=32,constant=True)
    # A longer constant-gain fit need not identify the END gain of varying
    # contact. Keep its signed bias, rather than treating lower variance as gain.
    constant_bias=float(constant['k']-theta[0])
    assert abs(constant_bias)>100.
    Bb,Uu,zz=design(rows,32);Ur=residualize(Bb,Uu)
    gain_row=np.linalg.pinv(Ur)[0];exact_variance=float(.04**2*(gain_row@gain_row))
    np.testing.assert_allclose(r['cov'][0,0],exact_variance,rtol=1e-9)
    rng=np.random.default_rng(73651);samples=[];nominal=[];linear_differences=[]
    for _ in range(250):
        noisy=rows.copy();perturb=rng.normal(0.,.04,32)
        for j,noise in enumerate(perturb):noisy[-32+j,3]+=noise+(decay*(noisy[-33+j,3]-rows[-33+j,3]))
        f=fit(noisy,nw=32);samples.append(f['k']);nominal.append(f['cov'][0,0]/max(f['chi2'],1.))
        linear_differences.append(float(f['k']-(theta[0]+gain_row@perturb)))
    variance=float(np.var(samples,ddof=1));predicted=float(np.mean(nominal))
    print('linear check',max(abs(np.array(linear_differences))),'empirical/predicted',variance/predicted)
    assert max(abs(np.array(linear_differences)))<1e-6
    # Exact covariance equality is the implementation check; increase this
    # cheap vectorised Gaussian sequence to reduce finite Monte-Carlo scatter.
    mc=rng.normal(0.,.04,(10000,32))@gain_row
    mc_variance=float(np.var(mc,ddof=1));assert abs(mc_variance/exact_variance-1)<.05
    # Causal prefix: arbitrary future changes cannot change a prefix fit.
    longer=np.r_[rows,np.full((10,7),1e9)]
    np.testing.assert_array_equal(fit(longer[:len(rows)],nw=32)['coefficients'],r['coefficients'])
    result=dict(exact_coefficients=r['coefficients'].tolist(),expected=theta.tolist(),
                long_constant_end_bias=constant_bias,noise_variance_empirical=variance,
                noise_variance_predicted=predicted,
                exact_gain_row_variance=exact_variance,vectorized_mc_variance=mc_variance,
                causal_prefix=True,probe_constraints=True,probe=probe,
                nominal_relative_sd_at_original_budget=.04/(.16740639598337909*np.sqrt(probe['information'])))
    target=Path(__file__).resolve().parent/'results';target.mkdir(exist_ok=True)
    (target/'local_verification.json').write_text(json.dumps(result,indent=2),encoding='utf8')
    print(json.dumps({k:v for k,v in result.items() if k!='probe'},indent=2))


if __name__=='__main__':main()
