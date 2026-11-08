# Task: can each method track rib-crossing stiffness changes faster and more accurately, and does that improve constant force there?

Simulation only. Do not read project_records/HANDOFF.md. Read `research/identification_module_scheme_20261108.md` section 6; reuse `experiments/rossi_failure_study/memory_segment_eval.py` (gates-off host, baseline check) and `memory_filters2.py`, `memory_filters.py`, `memory_closed_loop.py` (`scene_for`, seeds).

The unit of evaluation is one stiffness transition, not the episode. Earlier reports pooled segments and emphasised plateau noise, which does not answer the question.

## Transitions

From true ell = ln|K(20.8)| smoothed over 0.1 s, in contact and t>=1: contiguous runs with |d ell/dt| > 0.2 /s, merging gaps under 0.1 s. Keep runs with |delta ell| >= 0.15. Classify rising or falling, and small (<0.3) or large (>=0.3). Report the count, duration and size distribution per scene set. Exclude transitions with contact loss before their start.

## Metrics per transition

Let ell_hat be ln of the consumed k_good; also do k_bar. Let b = median(ell_hat - ell) over the 0.5 s before the transition (the method's own local offset; if unavailable, use the episode median and count such cases).
1. Tracking speed: lag tau = argmin over tau in [0, 1.5] s of mean|ell_hat(t) - ell(t-tau) - b| on [t0, t1+0.5]; and settling time, the first time after t1 from which |ell_hat - ell - b| <= 0.1 holds for 0.2 s (censor at 2 s and report the censored share).
2. Transition accuracy on [t0, t1+0.5]: most negative and most positive (ell_hat - ell - b), and the raw (without b) most negative value ln(k_good/k), since under-estimation is the unstable side.
3. Plateau accuracy, per episode, from `memory_segment_eval.py` plateau ticks: median and robust sd of ln(k_good/k).
4. Force on [t0, t1+0.5] and recovery: peak |F-8| N, positive and negative extremes, RMS of force error in the 15-35 Hz band over this interval relative to the episode's plateau value (oscillation), and recovery time after t1 until |F-8| <= 0.5 N holds for 0.2 s (censor at 2 s). Also baseline it against the same metric computed on matched plateau intervals of the same length, so the transition's own effect is visible.

## Methods, closed loop

Update every 6 ticks with no longer hold, so no method carries the 0.192 s delay. Filters use gain variance R = 4 s^2 for overlap; use tempered likelihood (power 1/4) for IMM and EWA because the innovations are correlated, and say so. Use theory-fixed parameters from `memory_filters2.py`, no tuning.

References: `on_none` (adopted Method, nothing attached, must reproduce `Huber-true dev`), `off_none`, `oracle` (window estimate = 1.2*truth through the gates-off host; perfect tracking, the ceiling for force).
Memory on the gates-off host: Kalman tau=2, Kalman tau=0.2, IMM, Lepski, EWA, forgetting with floor, forgetting with distance floor, Huang ARD fixed precision.
Other estimators as they are in the adopted comparison: `D4`, `Rossi-ID`, `Kalman` and `MHE` via `B.build`.

## Episodes

Base 520000, seeds as stage `dev3`, probe mass present, scene A `breath_rib`: nominal all 36 cells; `slow_thin` all 36 cells. Scene `rib` alone is not needed. Keep the catastrophic episodes (peak > 36.5 N) out of the transition metrics but count them per method and list which cells, so they do not dominate medians.

## Report

Under 2 pages in Chinese, final message only. Baseline check. Transition counts. Four tables, one per question, each method by rising/falling (large transitions, nominal and slow_thin separately): medians and P90, with the small-transition result in one line. Then a plain answer, per question: which methods track faster or more accurately than `on_none` and `off_none`, by how much, and how close force gets to `oracle` at transitions. Say where differences are too small to matter or too few to tell. State only what was measured. No fresh seeds, no winner selection.

Files: `experiments/rossi_failure_study/rib_transition_eval.py`, results `results/rib_transition_eval.json`. Up to 8 workers, `python -X utf8`. Edit no other repo file. No commit, no push.
