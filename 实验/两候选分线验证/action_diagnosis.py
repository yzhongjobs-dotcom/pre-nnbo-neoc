"""A finite target-gain dead zone versus the original B0 rate action."""
import json
import numpy as np
from experiment import ROOT,BASE
from core import encode,cost_gradient
from rate_core import rate_gradient


def main():
    d=np.load(ROOT/'training72_contexts.npz');cfg=np.array(json.loads((BASE/'tuned_fixed.json').read_text())['cfg']);cfg[8]=25.;rows=[]
    for i in np.random.default_rng(1921).choice(len(d['x']),size=32,replace=False):
        cc=cfg.copy();cc[7]=d['p'][i,9];q=np.r_[d['q'][i],d['sensed'][i]];target=cfg[:4].copy();target[1]=800.
        aa=np.tile(encode(np.r_[target,0.]),2);rate=np.zeros(10);rate[[1,6]]=.99
        _,g=cost_gradient(d['t'][i],d['y'][i],d['zr'][i],q,d['anchor'][i],cc[7],cc,aa,d['p'][i],True)
        _,gr=rate_gradient(d['t'][i],d['y'][i],d['zr'][i],q,d['anchor'][i],cc[7],cc,rate,d['p'][i])
        rows.append({'context':int(i),'target_Kn_gradient_norm':float(np.linalg.norm(g[[1,6]])),'rate_Kn_gradient_norm':float(np.linalg.norm(gr[[1,6]]))})
    out={'rows':rows,'target_zero_fraction':float(np.mean([r['target_Kn_gradient_norm']<1e-9 for r in rows])),'rate_zero_fraction':float(np.mean([r['rate_Kn_gradient_norm']<1e-9 for r in rows])),'scope':'32 related training contexts, fixed initial gains; diagnostic, not new independent control trials'}
    (ROOT/'action_diagnosis.json').write_text(json.dumps(out,indent=2));print({k:v for k,v in out.items() if k!='rows'})


if __name__=='__main__':main()
