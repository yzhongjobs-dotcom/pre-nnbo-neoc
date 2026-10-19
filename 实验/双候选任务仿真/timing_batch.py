"""Run publication-delay comparisons sequentially on an otherwise idle CPU."""
import argparse,json,subprocess,sys
from pathlib import Path

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--prefix',default='');args=ap.parse_args()
    root=Path(__file__).resolve().parent
    weights=json.loads((root/'baseline_axis_tuning.json').read_text())['selected']
    jobs=[('candidate1',None,'timed_candidate1.json'),('warm_value',None,'timed_candidate1_without_neighbor.json'),
          ('full_value',None,'timed_full_value.json'),('candidate2','candidate2_seed16111.pt','timed_candidate2.json'),('jin_ms',None,'timed_jin_ms.json')]
    for method,model,out in jobs:
        command=[sys.executable,str(root/'run.py'),'--method',method,'--n','8','--seed','16301','--timed','--output',args.prefix+out]
        if model:command+=['--model',model]
        if method=='jin_ms':command+=['--baseline-axis-regularizers',*map(str,weights)]
        subprocess.run(command,check=True)
