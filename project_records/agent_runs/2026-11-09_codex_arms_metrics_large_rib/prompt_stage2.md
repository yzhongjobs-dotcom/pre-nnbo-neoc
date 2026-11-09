# Continue: correct the transition detector and bound the windows. Re-score only.

Acceptance found a defect. `R.transitions` builds runs from abs(rate) > .2 and merges gaps < .1 s regardless of direction, so a rise and the following fall merge into one event with a small net delta. Example, case 36: crossings every ~.6 s with |d ln k| ~ 1, scored as (1.17-5.6 s, delta -.2). Case 39: 16 crossings, 1 event. Case 1 (scene A): (1.4-7.82 s, delta -.19). This affects all scenes.

I moved your earlier report to `project_records/agent_runs/2026-11-09_codex_arms_metrics_large_rib/report_stage1.md`. Do not write reports into `experiments/`.

Same constraints: no commit, no push, no research-doc edits, no retuning, no closed-loop or replay simulation. Only `score` and `tables` stages reuse saved ticks.

## Changes

1. `R.transitions`: detect rising runs (rate > .2) and falling runs (rate < -.2) separately. Merge two runs only if they have the same sign, the gap is < .1 s, contact holds across the gap, and no opposite-sign run lies between them. Boxcar, rate threshold, |delta| >= .15 and the t >= 1 s and contact rules stay as they are.

2. In `R.score`, the next change bounds each event's windows. next_t0 = t0 of the next detected transition with |delta| >= .15, kept or excluded.
   - `raw_under_magnitude`, `raw_over_magnitude`, `under_fraction`, `missing_fraction`: window [t0, min(t1+.5, next_t0)).
   - `recovery_raw`: `hold_time` gets a limit. A qualifying in-band run of .2 s must complete before next_t0. Otherwise censored, with `censored_by` = `next_change`, `episode_end` or `horizon` (full 2 s observed, not recovered). observed_horizon = min(2, limit - t1). A zero or negative horizon is still counted, as censored at 0.
   - `change_mae` unchanged. Lag and force fields unchanged.

3. `summarize`: recovered, unrecovered (`horizon`), truncated by episode end, truncated by next change. The four sum to n. KM median as now. Assert the sum.

4. Table 3: replace 回合结束截断 by two columns 下次变化截断 and 回合结束截断. Add one sentence: 下一次刚度变化开始前未恢复的事件按右删失处理，最大低估窗口也在下一次变化开始时结束。Table 3's event-count note must match.

5. Tests in `arms_reproduction_check.py`: a synthetic ln k triangle wave gives one event per flank with the right sign; a monotone ramp with a .05 s sub-threshold gap stays one event; recompute `change_mae`, under and recovery for 3 real A_large events independently.

6. Re-run `identification_completion.py --stage score` and `--stage tables`, and `rib_contrast_evidence.py`'s scored section.

## Report (final message only, at most 1 page)
Event counts by family and direction before and after. The new Table 3 lines for Huber-D4 three arms, D4 three arms, and Rossi-ID: A_large and existing large, k_good and raw. Test results. Anything that failed. Confirm the closed-loop JSON is untouched by hash.
