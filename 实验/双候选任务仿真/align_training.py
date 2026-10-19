"""Correct descriptor timestamps in the first causal training collection.

Robot/force histories and fits are unchanged. Only propagate the held estimate
from its 100 ms update time to each 20 ms actor query, using its own velocity.
New collections in run.py already perform this propagation.
"""
import hashlib,json
from pathlib import Path
import numpy as np

root=Path(__file__).resolve().parent;source=root/'training_contexts.npz';target=root/'training_aligned_contexts.npz'
if target.exists():raise FileExistsError(target)
with np.load(source) as data:
    x=data['x'];estimate=data['estimate'].copy();truth=data['truth']
age=(np.arange(len(x))%80%5)*.02
estimate[:,7]+=estimate[:,8]*age;estimate[:,9]+=estimate[:,11]*age
np.savez_compressed(target,x=x,estimate=estimate,truth=truth)
(root/'training_alignment.json').write_text(json.dumps({'source':source.name,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
    'target':target.name,'target_sha256':hashlib.sha256(target.read_bytes()).hexdigest(),'contexts':len(x),'physical_cases':48,
    'operation':'propagate held estimate by its own speed and surface velocity to actor timestamp; no new or future observations'},indent=2))
