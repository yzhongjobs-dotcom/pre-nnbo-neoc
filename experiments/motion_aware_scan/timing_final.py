"""Computation time of the final method (single process): python -X utf8 timing_final.py -> results_v2/timing_final.json"""
import json, os
import numpy as np
import common, cohort_eval as CE
if __name__ == '__main__':
    sc = CE.cohort('final', 4); ts = []; per = []
    for fam, sid, p, T in sc:
        pol = CE.build('Proposed', p); log, pol = common.run(pol, p, T, sid)
        t = np.array(pol.t_call[1:])*1e3; ts += list(t); per.append(float(np.percentile(t, 99)))
    ts = np.array(ts)
    out = dict(version='Proposed-IMMb (meas_bound=1.0, enc_y=True)', scenes=len(sc), calls=int(len(ts)), p50_ms=float(np.percentile(ts, 50)),
               p99_ms=float(np.percentile(ts, 99)), worst_scene_p99_ms=max(per), max_ms=float(ts.max()), calls_over_8ms=int((ts > 8).sum()),
               note='the call over 8 ms is the first window fit in a fresh process (loading compiled numba code)')
    json.dump(out, open(os.path.join(common.HERE, 'results_v2', 'timing_final.json'), 'w'), indent=1); print(out)
