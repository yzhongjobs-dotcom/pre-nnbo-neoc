"""Replay the adapter's critic; distinguish inverse-range from solver failure.

Source Eq52 with the declared nonzero-force shift and quadratic basis gives
S(d)=f-Fd+(A*d+b)/(B*d**n), d>0. When A>0,b<0,n>=1 its positive term
has a finite maximum. A sufficiently negative force error has NO real root.
This diagnoses our adaptation, not the original zero-target-force NNBO.
"""
import json
import numpy as np
from experiment import ROOT
from baselines import NNBO


def main():
    d=json.loads((ROOT/'candidate1_test_w1.json').read_text());traces=np.load(ROOT/'candidate1_test_w1_traces.npz');cfg=np.array(d['cfg']);rows=[]
    for saved in d['results']:
        if saved['method']!='nnbo_adapted':continue
        i=saved['case'];p=np.array(d['parameters'][i]);cfg[7]=p[9];tr=traces[f'{i}_nnbo_adapted'];zr=tr[0,2]
        controller=NNBO(np.array(saved['hc']),zr,cfg,omega=d['config']['omega']);outside=0;outside_failures=0
        for row in tr:
            previous=controller.failures;controller.decide(row[0],row[1:11],row[28:31],dt=.01)
            A=controller.w[0];b=A*(controller.surface-controller.r)+.5*controller.w[1]*controller.r
            if A<=0 or b>=0 or controller.n<1:continue
            if controller.n-1<1e-8:maximum=A/controller.b
            else:
                peak=-controller.n*b/((controller.n-1)*A)
                maximum=np.exp(np.log(-b/(controller.n-1))-np.log(controller.b)-controller.n*np.log(peak))
            impossible=row[28]-3+maximum < -1e-8
            outside+=impossible;outside_failures+=impossible and controller.failures>previous
        np.testing.assert_allclose(controller.w,saved['final_critic'],rtol=1e-12,atol=1e-12)
        assert controller.failures==saved['reference_failures']
        rows.append({'case':i,'reference_failures':controller.failures,'outside_inverse_range':int(outside),'outside_and_failed':int(outside_failures)})
    result={'rows':rows,'total_failures':sum(r['reference_failures'] for r in rows),'proven_outside_range':sum(r['outside_inverse_range'] for r in rows),'outside_and_failed':sum(r['outside_and_failed'] for r in rows),'scope':'Our nonzero-force quadratic-critic adaptation only. Original paper f>=0, Fd=0 has a different inversion range.'}
    (ROOT/'nnbo_inverse_range_diagnosis.json').write_text(json.dumps(result,indent=2));print({k:v for k,v in result.items() if k!='rows'})


if __name__=='__main__':main()
