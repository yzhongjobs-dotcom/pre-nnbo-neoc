"""Post-test numerical sensitivity; never counted as new test cases."""
import json
from pathlib import Path
import numpy as np
import mechanism as m

ROOT=Path(__file__).resolve().parent


def main():
    result=json.loads((ROOT/'results.json').read_text())
    rows=result['rows'];J=np.array([r['J'] for r in rows])
    selected=sorted(set([0,1,int(np.argmin((J[:,1]-J[:,3])/J[:,1])),
                         int(np.argmax((J[:,1]-J[:,3])/J[:,1]))]))
    out=[]
    for i in selected:
        r=rows[i];ctx=r['context']
        costs,actions,values,_,_=m.run_case(ctx,result['protocol']['fixed'],nmc=32768,
            seed=11602+i,nodes_weights=m.sobol_nodes(14,seed=11701))
        out.append(dict(case=i,original_J=r['J'],refined_J=costs.tolist(),
                        original_actions=r['actions'],refined_actions=actions))
        print(out[-1],flush=True)
    data=dict(selection='First two plus worst/best full-Bayes relative gain, selected AFTER test',
              qmc_power=14,qmc_scramble_seed=11701,mc_per_mode=32768,rows=out,
              scope='Numerical sensitivity only; no new independent performance examples')
    (ROOT/'integration_audit.json').write_text(json.dumps(data,indent=2)+'\n')


if __name__=='__main__':main()
