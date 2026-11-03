"""Checkpoint integrity, delayed causality, and version-compatible replay."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[key]='1'
import json
import hashlib
from pathlib import Path
import numpy as np
from run_joint import common,CE
from response_policy import ResponsePolicy
from dynamic_compliance import reconstruct_increment

HERE=Path(__file__).resolve().parent
OUT=HERE/'results'


def main():
    inputs=('shared_response_v1','shared_response_history_v2','shared_response_uncertainty_v3')
    summary=[];sources=[]
    for stem in inputs:
        d=json.loads((OUT/(stem+'.json')).read_text())
        assert len(d['rows'])==d['protocol']['expected_rows']
        for r in d['rows']:
            assert r['exact_common_startup']
            assert hashlib.sha256((OUT/r['measurements']).read_bytes()).hexdigest()==r['measurements_sha256']
            assert all(f['time']-f['fit_time']>=.032-1e-10 for f in r['fits'])
            summary.append(dict(batch=stem,family=r['family'],model=r['model'],rmse=r['metrics']['rmse'],
                baseline=r['baseline_metrics']['rmse'],peak=r['metrics']['peak'],used=r['used'],
                maximum_recorded_fit_ms=max(f['fit_seconds']*1000 for f in r['fits'])))
            if (r['family']=='T3m' and r['model']=='shared_dynamic') or r['model']=='shared_history_uncertain':
                sources.append(r)
    replay=[]
    for r in sources:
        m=np.load(OUT/r['measurements'])['measured'];name=r['model']
        pol=ResponsePolicy(response_model='shared',conversion='dynamic',history=name.startswith('shared_history'),
                           uncertainty=name=='shared_history_uncertain',single_model=False,**CE.FINAL)
        pol.reset(dict(common.P.PLATFORMS['PC125']))
        maximum=0.
        for j in range(len(m)):
            # Only the measured prefix exists. Truth and future columns are NaN.
            log=np.full((j+1,16),np.nan);log[:,[0,1,2,5,11,12,13]]=m[:j+1]
            log[j,11]=m[j-1,4] if j else m[0,2]
            command=pol(common.P.Obs(log,j+1))
            maximum=max(maximum,abs(command-m[j,4]))
        assert maximum==0.
        replay.append(dict(family=r['family'],model=name,max_command_error=maximum))
    # Read only complete same-input summaries; overlapping windows are not
    # independent trials and are not used for a significance claim.
    same=json.loads((OUT/'shared_response_same_input_matched_v2.json').read_text())
    for p,h in same['protocol']['input_hashes'].items():
        assert hashlib.sha256((OUT/p).read_bytes()).hexdigest()==h
    physical=json.loads((OUT/'shared_response_verification.json').read_text())
    result=dict(complete_closed_loop_rows=len(summary),rows=summary,
        compatible_prefix_replay=replay,independent_physics_checks=physical,
        same_input_prediction=same['summary'],
        variance_limits='conditional local model; initial-state bias and model mismatch not certified',
        computation_limits='32-ms publication lag simulated; callback is a synchronous prototype, not certified asynchronous deployment',
        no_new_candidate_adopted=True)
    target=OUT/'shared_response_checkpoint.json'
    target.write_text(json.dumps(result,indent=2),encoding='utf8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
