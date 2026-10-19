"""Syntax/data/model/link validation and reproducibility manifest.

Does not rerun control experiments or certify their scientific assumptions.
"""
import argparse,ast,hashlib,json,platform,re,sys
from pathlib import Path
import numpy as np
import scipy,numba,torch,casadi,matplotlib

def main(write=False):
    root=Path(__file__).resolve().parent;counts={}
    for path in root.glob('*.py'):ast.parse(path.read_text());counts['python']=counts.get('python',0)+1
    for path in root.glob('*.json'):
        json.loads(path.read_text(),parse_constant=lambda s:(_ for _ in ()).throw(ValueError(s)))
        counts['json']=counts.get('json',0)+1
    for path in root.glob('*.npz'):
        with np.load(path) as data:
            for key in data.files:
                if np.issubdtype(data[key].dtype,np.number):assert np.isfinite(data[key]).all(),(path,key)
        counts['npz']=counts.get('npz',0)+1
    def finite(v):
        if isinstance(v,torch.Tensor):assert torch.isfinite(v).all()
        elif isinstance(v,dict):
            for x in v.values():finite(x)
        elif isinstance(v,(tuple,list)):
            for x in v:finite(x)
    for path in root.glob('*.pt'):finite(torch.load(path,map_location='cpu',weights_only=True));counts['models']=counts.get('models',0)+1
    for path in root.glob('*.md'):
        for target in re.findall(r'\]\(([^)]+)\)',path.read_text()):
            target=target.split('#')[0]
            if target and '://' not in target:assert (path.parent/target).exists(),(path,target)
    manifest={'date':'2026-10-19','python':sys.version,'platform':platform.platform(),
        'packages':{k:m.__version__ for k,m in [('numpy',np),('scipy',scipy),('numba',numba),('torch',torch),('casadi',casadi),('matplotlib',matplotlib)]},
        'cpu_quota':Path('/sys/fs/cgroup/cpu.max').read_text().strip(),'threads_per_process':1,
        'scope':'synthetic planar task, finite numerical approximation; statistical and modelling limits in README.md',
        'validation_counts':counts,'sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(root.iterdir()) if p.is_file() and p.name!='environment.json'}}
    target=root/'environment.json'
    if write:target.write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    else:
        stored=json.loads(target.read_text());assert stored['sha256']==manifest['sha256'],'artifact changes after manifest generation'
    print(json.dumps({'validated':counts,'manifest_written':write,'sha256_files':len(manifest['sha256'])}))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write-manifest',action='store_true');args=ap.parse_args();main(args.write_manifest)
