# Task: independent check of one research note and its code. Find errors; do not propose methods.

Read `research/identification_module_scheme_20261108.md`, then `experiments/rossi_failure_study/memory_filters.py` and `experiments/rossi_failure_study/memory_closed_loop.py`. Do not read project_records/HANDOFF.md. Report under 2 pages in Chinese, in your final message only: for each item correct / wrong / correct only under conditions, with the reason and the corrected value if wrong.

1. Section 2 formulas: tau = sqrt(N/q), sd^2 = N/(2 tau), lag rho*tau, tau* = (N/(4 rho^2))^(1/3) and the three RMSE values; the claim tau >= 1.9 s from half width 0.09 at 2 sd.
2. Huang ARD closed form alpha* = max(nu^2 - S, 0) for the scalar case of eq (14) in the paper at `literature/manual_inbox/Novel sparseness-inducing*.pdf` (use pdftotext). Check the gain 1 - r/nu^2, the 32% firing rate, and whether 0.33 r is the injected variance per step or per firing.
3. `HuangARD` and `Forgetting` code: do they implement what the note says? Is the eq (15)/(17) precision learning implemented correctly in scalar form?
4. `attach` in memory_filters.py: is the overlap inflation R*overlap consistent with hop=6 and nw=24? Does the Method still get raw-window behaviour where the note says it does (chi gate and jump check)? Any path where the filter output enters the gates differently from what the note states? Is the raw `usable` gate reasonable?
5. `boot` and the delta-L intervals in `memory_closed_loop.py`: does the corrected clustering do what the note says? Recompute two intervals from `experiments/rossi_failure_study/results/memory_closed_loop_test.json` (Huber-D4 minus Proposed, nominal, all scenes; Rossi-ID minus Proposed, slow_thin) and compare with the note.
6. Is "f=1.0 is at the cliff edge" supported by `experiments/rossi_failure_study/results/gain_need_sweep.json`? Check the numbers in the Section 1 table.

Constraints: no method suggestions. Do not edit any repo file. No commit, no push. Leave `00_基线快照_2026-10-14/` alone. No closed-loop batches.
