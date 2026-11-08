# Task: compare Huber-D4 with and without adaptive memory, gates on and off, evaluated per segment

Simulation only. Read `research/identification_module_scheme_20261108.md` sections 1, 3, 5 for context; do not read project_records/HANDOFF.md. Code: `experiments/rossi_failure_study/memory_filters.py`, `memory_closed_loop.py`, `tail_channel_diagnosis.py`, `estimator_swap.py` (`_stiffness_est`), `../motion_aware_scan/proposed.py`.

Earlier comparisons reported only whole-episode numbers, and one earlier reference was invalid because a wrapper changed it. Phase 1 here: a small set, to see what gates do and where memory helps. No method suggestions, no statistics beyond what is asked.

## Configurations (2x2)

Host estimator: Huber-D4 window, unchanged.
- Gates ON, no memory: the adopted Method with nothing attached.
- Gates ON, memory: `MF.attach(pol, 'forgetting', lam_max=.998)` as it stands.
- Gates OFF, no memory: every usable window estimate goes straight through.
- Gates OFF, memory: every window estimate goes through the forgetting filter, then straight through.

Gates OFF means, on the policy instance only: `k_good` takes the output each update (keep the existing centering factor); the precision test, the chi test and the b <= ln1.2 test no longer block anything; the k_check/k_good feed-forward jump rule is off. Keep unchanged: `k_bar = max(k + c*sigma, (1-rho_d)*k_bar)`, the first-acquisition logic, contact-loss supervision, servo identification. In the memory case use sigma = the filter's, as in `attach`, and feed all windows (`usable` always true). State exactly what you overrode.

Do not reuse `Window` from memory_filters as a reference. Its output is held and its sigma doubled.

## Episodes

Base 520000, seeds as `memory_closed_loop.py` stage `dev3` (SEED+q, q the index in `O.cells()`), probe mass present:
- Scene A `breath_rib`, nominal: all 36 cells.
- Scene A `breath_rib`, `slow_thin` variant of `scene_for`: all 36 cells.
- Scene B `breath_fast`: all 12 cells.

Check first: Gates ON, no memory must reproduce `results/memory_closed_loop_true_base.json` entries named `Huber-true dev` (L and peak) for all these episodes. If not, stop and say so.

## Segment evaluation

Use truth only for scoring. Per tick, from t>=1 and in contact, with k = logged |K(20.8)|, r = |d ln k/dt| with k smoothed over 0.2 s:
- `transition`: r > 0.2 /s, plus 0.1 s either side.
- `plateau`: r < 0.05 /s and not within 0.3 s of a transition.
- `event`: inside `sc['events']` plus 0.3 s after (scene B). Events take priority.
- other ticks: not scored by segment.

For each configuration, scene set and segment report: share of time; for e_g = ln(k_good/k) and e_bar = ln(k_bar/k): median, robust sd (IQR/1.349), fraction of ticks with k_good/k < 1.0 and < 0.85; force error RMSE and P99 of |F-8| N; number of episodes with peak > 15 N attributed to that segment (segment containing the peak). Also the raw window estimate e_raw for the same segments, so that estimator noise and consumed noise can be compared. Also whole-episode L mean/median and max peak.

Pool ticks across episodes but give each episode equal weight in the medians and fractions (average per-episode values).

## Report

Under 2 pages in Chinese, in your final message: the baseline check; the exact overrides; one table per scene set, configuration by segment; what the data show about (1) what gates contribute and (2) where memory helps or hurts, with and without gates. State only what was measured.

Script `experiments/rossi_failure_study/memory_segment_eval.py`, results `experiments/rossi_failure_study/results/memory_segment_eval.json`. Use `python -X utf8`, up to 8 workers. Edit no other repo file. No commit, no push.
