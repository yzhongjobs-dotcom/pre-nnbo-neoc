"""Train fixed gains, freeze a protocol, then run new mechanism cases."""
from pathlib import Path
import hashlib
import json
import time
import numpy as np
import mechanism as m

ROOT=Path(__file__).resolve().parent


def main():
    # Test contexts generated only AFTER recording the completed protocol.
    train=m.sample(16,11501)
    ol=np.mean([m.open_loop_values(c,m.prepare(c)) for c in train],axis=0)
    fixed=list(map(int,np.unravel_index(np.argmin(ol),ol.shape)))
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'mechanism.py',ROOT/'evaluate.py']}
    protocol=dict(train=train,fixed=fixed,n_test=64,test_seed=11601,noise_seed=11602,
                  planning_qmc_power=12,planning_qmc_seed=11001,
                  n_mc_per_mode=4096,gains=m.GAINS.tolist(),modes=m.MODES.tolist(),
                  methods=['trained_fixed_two_stage','frozen_belief_plan','future_encoder_plan','future_encoder_force_plan'],
                  source_sha256=hashes,
                  scope='New local two-decision mechanism cases, not B0/continuous scan or field-prior superiority',
                  latency='Not injected into plant; score and per-context construction measured separately; not deployment performance')
    (ROOT/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    # JIT all branches on a TRAINING case before collecting any timing.
    m.run_case(train[0],fixed,nmc=8,nodes_weights=m.sobol_nodes(5))
    test=m.sample(protocol['n_test'],protocol['test_seed'])
    rows=[];nodes=m.sobol_nodes(12)
    start=time.perf_counter()
    for i,ctx in enumerate(test):
        J,actions,values,timing,setup=m.run_case(ctx,fixed,nmc=4096,seed=11602+i,nodes_weights=nodes)
        rows.append(dict(context=ctx,J=J.tolist(),actions=actions,planning_ms=timing,setup_ms=setup,
                         action_scores=[v.tolist() for v in values]))
        print(i,np.round(J,7).tolist(),actions,flush=True)
    values=np.array([r['J'] for r in rows]);rel=(values[:,1:2]-values[:,2:])/values[:,1:2]
    rng=np.random.default_rng(11603);indices=rng.integers(0,len(rows),(10000,len(rows)))
    ci=np.quantile(rel[indices].mean(axis=1),[.025,.975],axis=0)
    result=dict(protocol=protocol,rows=rows,mean_J=values.mean(axis=0).tolist(),
                relative_gain_vs_frozen=rel.mean(axis=0).tolist(),relative_gain_vs_frozen_CI=ci.tolist(),
                different_first_actions_vs_frozen=[sum(r['actions'][j]!=r['actions'][1] for r in rows) for j in [2,3]],
                planning_median_ms=np.median([r['planning_ms'] for r in rows],axis=0).tolist(),
                preparation_median_ms=float(np.median([r['setup_ms'] for r in rows])),
                elapsed_s=time.perf_counter()-start,
                limitations='Known hypothesis set, exact encoders, no gain-rate/energy/delay/contact-loss model, no full field baseline')
    (ROOT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez_compressed(ROOT/'paired.npz',J=values)
    print('SUMMARY',json.dumps({k:v for k,v in result.items() if k not in ['protocol','rows']},indent=2),flush=True)


if __name__=='__main__': main()
