# Task: locate which channel causes abnormal force peaks, on a few representative cases

Simulation only. The adopted Method is `Huber-D4`: `B.build('Huber-D4', O.chosen()['Huber-D4'], sc)` with `import overlay_compare as O; O.setup(); B = O.B`. See `experiments/rossi_failure_study/gain_need_sweep.py` for how an episode is run and truth obtained (`P3.episode(sc, pol, seed, diagnostic=True)`, `o.truth_k`), and `estimator_swap.py` and `../motion_aware_scan/proposed.py` for where k_bar and k_good are used. Do not read project_records/HANDOFF.md.

Phase 1 only: find the cause on a few cases. No statistics, no method suggestions.

**Cases.** Scene A cell `('breath_rib','rib','neural_white',28.,6)` with `B.scene(cell, 674500, 0.)` and seed `O.sensor_seed(cell)`, peak 35.5 N. Then from `experiments/rossi_failure_study/results/memory_closed_loop_true_base.json` take the highest-peak episode in each of: name `Huber-true` variant `slow_thin`; name `Huber-true six` scene `combo`; name `Huber-true six` scene `cough`. Build these the way `memory_closed_loop.py` stage `true_base` does (base 774500, `scene_for`, same seed index).

**Runs, each case:**
- `own`: unchanged. Its peak must match the stored value; if not, stop and say so.
- `kbar`: k_bar = 1.2*truth_k every tick; k_good from the estimator.
- `kgood`: k_good = 1.2*truth_k every tick; k_bar from the estimator.
- `ff_off`: feed-forward increment zero; stiffness channels untouched.

Override on the policy instance, not in shared code, and confirm once that the command really uses the replaced value.

**Report, under one page in Chinese, in your final message:** a table of peak force and L per case and run; for each `own` case, around the peak: k/k_bar, k/k_good, what the force did just before (oscillation frequency if any, contact loss, a jump in k_good or k_bar). Say plainly which channel explains each case, whether the four cases agree, and what is unexplained.

Script: `experiments/rossi_failure_study/tail_channel_diagnosis.py`. Use `python -X utf8`. Edit no other repo file. No commit, no push.
