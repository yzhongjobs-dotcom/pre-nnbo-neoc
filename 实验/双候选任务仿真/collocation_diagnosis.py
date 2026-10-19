"""Reproduce the rejected global polynomial-critic development attempt."""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
import physics as p
from critic import Critic

if __name__=='__main__':
    rows=[]
    for T in [.05,.1,.2,.3]:
        p.HORIZON=T
        eta=np.array([1e5,1.4e5,.15,.2,2.,3.,.025,0.,0.,0.,.001,.01]);x,eta=p.initial(eta);c=Critic(eta);c.setup(eta)
        res=least_squares(lambda w:c.residual(w)[0],np.zeros(66),jac=lambda w:c.residual(w)[1],max_nfev=500,x_scale='jac')
        rows.append({'horizon':T,'evaluations':res.nfev,'maximum_residual':float(np.max(abs(res.fun))),'success':bool(res.success)})
    p.HORIZON=.3
    Path(__file__).with_name('collocation_diagnosis.json').write_text(json.dumps({'scope':'No acceptable root found by this basis/node/solver attempt; not a proof that B0 has no root. Not used in primary experiments.','results':rows},indent=2))
