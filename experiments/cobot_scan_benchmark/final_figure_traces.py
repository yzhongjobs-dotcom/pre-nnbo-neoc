"""Traces for the final experiment figures: fixed stiff baseline, Haninger, EtR-A, C1, C2 (seed 0)
on T3 scene 0 and T2 scenes 0 (soft push) and 2 (stiff push) of the seed-7200 cohort."""
import numpy as np
import complete_task_evidence as E
from fixed_extremes_eval import EXT

if __name__ == '__main__':
    out = {}
    for task, idx in (('T3_combined', [0]), ('T2_motion', [0, 2])):
        ss = E.scenes(task); sp = {s[0]: s for s in E.specs(task) if s[3] in (-1, 0)}
        for i in idx:
            p, T = ss[i]
            log, _, _ = E.episode(p, np.array(EXT['Fixed_stiff']), T, E.SEED + i, None)
            out[f'{task}_{i}_Fixed_stiff'] = log
            for m in ('EtR-A', 'EtR-A-Haninger', 'Candidate1_v2', 'Candidate2'):
                E.warm(sp[m]); pol, upd = E.build(sp[m]); E.bind(pol, p)
                log, _, _ = E.episode(p, pol, T, E.SEED + i, upd)
                out[f'{task}_{i}_{m}'] = log
            out[f'{task}_{i}_p'] = p
    np.savez_compressed('completed_tasks/final_figure_traces.npz', **out)
    print('saved', len(out))
