"""Read separate raw force/pose captures; never pair them with image episodes."""
from pathlib import Path
from urllib.request import urlopen
import csv
import hashlib
import io
import json
import numpy as np

HERE=Path(__file__).resolve().parent
BASE='https://raw.githubusercontent.com/suung23/FR5-for-RUS/c765e23be3b2ccebdd971b521d9e1f8c2662d9f4/'


def main():
    cache=Path('/tmp/r15-real-scan/cache');manifest=[];results=[]
    for stem in ['B_t3p0_pose','S_stiff']:
        sources={}
        for suffix in ['meta.json','samples.csv']:
            path=f'force_hold_validation/runs/{stem}_{suffix}'
            dst=cache/path;dst.parent.mkdir(parents=True,exist_ok=True)
            if not dst.exists():dst.write_bytes(urlopen(BASE+path,timeout=60).read())
            data=dst.read_bytes();sources[suffix]=data
            manifest.append(dict(path=path,bytes=len(data),sha256=hashlib.sha256(data).hexdigest()))
        meta=json.loads(sources['meta.json'])
        rows=list(csv.DictReader(io.StringIO(sources['samples.csv'].decode())))
        arr=lambda key:np.array([float(r[key]) if r[key] else np.nan for r in rows])
        t=arr('t_s');f=arr('mag_n');fn=arr('fn_n')
        xyz=np.array([arr(k) for k in ['fx_n','fy_n','fz_n']]).T
        pose=np.array([arr(k) for k in ['tip_x_m','tip_y_m','tip_z_m']]).T
        pose=pose[np.isfinite(pose).all(axis=1)]
        dt=np.r_[np.diff(t),0.]
        results.append(dict(label=stem,metadata=meta,rows=len(rows),duration_s=float(t[-1]-t[0]),
            median_record_interval_s=float(np.median(np.diff(t))),
            p95_record_interval_s=float(np.quantile(np.diff(t),.95)),
            force_min_N=float(f.min()),force_max_N=float(f.max()),
            time_weighted_mean_force_N=float(dt@f/dt.sum()),
            max_abs_signed_magnitude_vs_norm_N=float(np.max(np.abs(np.abs(f)-np.linalg.norm(xyz,axis=1)))),
            max_abs_signed_magnitude_vs_normal_N=float(np.max(np.abs(f-fn))),
            valid_pose_rows=len(pose),distinct_xyz_values=len(np.unique(pose,axis=0)),
            max_translation_from_initial_m=float(np.max(np.linalg.norm(pose-pose[0],axis=1))),
            marked_events=[dict(time_s=float(r['t_s']),event=r['event']) for r in rows if r['event']]))
    out=dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        scope='Two fixed captures from a separate force-hold experiment; not synchronized with the image export cohort',
        interpretation=['~1ms capture intervals do not imply 1kHz independent pose updates',
            'mag_n is signed norm; compare its absolute value with Euclidean norm',
            'The S_stiff filename does not establish a qualified stiffness calibration',
            'Unknown syringe surface motion prevents inferring causal stiffness from force/robot-position correlation alone'],
        sources=manifest,captures=results)
    (HERE/'force_audit.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(results,indent=2))


if __name__=='__main__':main()
