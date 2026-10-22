"""Field-style fixed-impedance baselines on the seed-7200 completed-task scenes.

As in Jin TMECH2023 (fixed small/large impedance at the method's own operating bounds),
the two fixed controllers use C1's own gain rule D = max(20, k/80), KPF = min(2/k, 1e-3),
M = 0.3, K = 0, ki = 30, evaluated at the softest (k = 1750 N/m, compliant) and stiffest
(k = 25000 N/m, stiff) local contact stiffness of the scene distribution."""
import json, csv, numpy as np
import complete_task_evidence as E

KS, KH = 1750., 25000.
EXT = {'Fixed_compliant': (0.3, max(20., KS/80), 0., 30., min(2/KS, 1e-3)),
       'Fixed_stiff': (0.3, max(20., KH/80), 0., 30., min(2/KH, 1e-3))}
TASKS = ('T1_static', 'T2_motion', 'T3_combined')

if __name__ == '__main__':
    rows = []
    for task in TASKS:
        ss = E.scenes(task)
        for name, g in EXT.items():
            for i, (p, T) in enumerate(ss):
                log, _, _ = E.episode(p, np.array(g), T, E.SEED + i, None)
                rows.append(dict(task=task, method=name, scene=i, **E.metrics(log, T, task)))
    # check against the archived Best_fixed implementation on scene 0 of T3
    arch = {(r['task'], r['method'], r['scene']): r for r in csv.DictReader(open('completed_tasks/per_scene.csv'))}
    out = {'rows': rows}
    json.dump(out, open('completed_tasks/fixed_extremes.json', 'w'), default=float, indent=0)
    for task in TASKS:
        for name in EXT:
            r = [x for x in rows if x['task'] == task and x['method'] == name]
            f = lambda k: np.mean([x[k] for x in r if x.get(k) is not None])
            print(task, name, 'rmse %.3f max %.3f peak %.3f worst %.2f nocontact %.3f' % (
                f('rmse_N'), f('max_abs_error_N'), f('peak_force_N'), max(x['peak_force_N'] or 1e9 for x in r), f('no_contact_s')))
