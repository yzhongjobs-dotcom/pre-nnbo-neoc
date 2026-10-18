"""Numerical refinement of the largest follow-up regression; no policy retuning."""
from pathlib import Path
import json,hashlib
import numpy as np
import model as m
from run import noises
HERE=Path(__file__).resolve().parent
fit=json.loads((HERE/'pure_gain_fit.json').read_text());data=np.load(HERE/'pure_gain_paired.npz')
ps=data['test_params'];ns=noises(ps,8502);cfg=np.array(fit['cfg'])
reg=data['learned_gain_only'][:,0]-data['continuous_gain_only'][:,0]
worst=int(np.argmax(reg));cases=sorted(set([0,2,worst]));rows=[]
latency=json.loads((HERE/'results.json').read_text())['timing']['learned_feedback']['p95_ms']/1000
for case in cases:
 for name in ['continuous_gain_only','learned_gain_only']:
  theta=np.array(fit['policies'][name])
  fine=m.episode(ps[case],ns[case],12,cfg,policy=theta,geometry_prior=True,latency=latency,dt=.00025,nq=129)[0]
  old=data[name][case]
  row={'case':case,'family':case%3,'method':name,'original_J':float(old[0]),'fine_J':float(fine[0]),
       'relative_difference':float(fine[0]/old[0]-1)}
  rows.append(row);print(row,flush=True)
out={'rows':rows,'worst_case':worst,'scope':'Post-hoc numerical diagnostic on unchanged policies and test cases; not additional independent data',
  'source_sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['model.py','audit_followup.py','pure_gain_fit.json']}}
(HERE/'followup_checks.json').write_text(json.dumps(out,indent=2)+'\n')
