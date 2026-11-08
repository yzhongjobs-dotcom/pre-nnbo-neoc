# Task: fix the H preprocessing defect and measure what it changes

Simulation only. Do not read project_records/HANDOFF.md. Read `project_records/agent_runs/2026-11-08_codex_low_estimate_forensics/output.md`, `experiments/rossi_failure_study/low_estimate_forensics.py`, `rib_transition_eval.py`, `memory_segment_eval.py`, and `contact_candidates.py` class `Preprocess`.

## The defect (verified in code)

`Preprocess._log`, HC branch: `y_out = y0 + ((y-y0) + tau_h*dy/TC)/mag_h`. mag_h and tau_h change online with the phase, the factor multiplies the displacement accumulated since the episode start, and rows are cached. So rows in one window use different factors, which injects steps into Y. A factor change of 1e-3 at y-y0 = 18 mm is 18 um, against a probe displacement of 4-10 um on hard tissue.

## The fix

Increment form with the current coefficients:
`y_out_j = y_out_{j-1} + (dy_j + tau_h*(dy_j - dy_{j-1})/TC)/mag_h`
A coefficient change no longer moves earlier rows or creates a jump. Implement on the estimator instance through a flag, leaving shared code behaviour unchanged when the flag is off. Check: at constant coefficients the two forms give the same increments to round-off; and the flag off reproduces stored baselines exactly. If the increment form needs a correction, derive it and say so. Check also whether `MEMORY` branch or any other path has the same defect, and whether Rossi-ID or other hosts go through it; report, do not change them.

## Phase 1: the 7 forensic cases

Runs: `own` (must reproduce stored peaks), `hfix`, `off_hfix` (gates off as in `memory_segment_eval.py`, plus the fix). Report peak and L, and whether the first low event still occurs and the recomputed k_hat of the analysed low windows. Also redo the low-tail table for k_hat < 0.5*c*k by sigma/k bin on A nominal and slow_thin (dev base).

## Phase 2: only if phase 1 shows a clear effect

Dev seeds (base 520000, stage `dev3` seeds): scene A nominal 36, scene A slow_thin 36, scene B 12. Configs: `on_none`, `on_hfix`, `off_none`, `off_hfix`, `off_hfix_forget` (forgetting with floor, 6-tick updates as in `rib_transition_eval.py`), `off_hfix_ewa` (same interface). Use `rib_transition_eval.py` metrics: large transitions rising and falling, four questions: lag and settling; low and high extremes with and without local offset; plateau bias and sd; force peak, oscillation ratio, recovery. Add whole-episode L mean/median, max peak, episodes above 15 N and above 36.5 N with the list of cells. Keep the catastrophic exclusion rule of that script but report the count.

## Report

Under 2 pages in Chinese, final message only. State exactly what changed in the code, checks, phase 1 table, low-tail table before and after, phase 2 tables, which failures remain with their cell ids. State only what was measured. No fresh seeds, no winner selection.

Files: `experiments/rossi_failure_study/h_fix_eval.py`, results `results/h_fix_eval.json`. Up to 8 workers, `python -X utf8`. Only add the flag in `contact_candidates.py`, default off; edit nothing else. No commit, no push.
