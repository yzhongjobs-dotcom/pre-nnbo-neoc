# Task: find why the stiffness estimate goes low before failures. Causal analysis on a few representative cases.

Simulation only. Do not read project_records/HANDOFF.md. Read `research/identification_module_scheme_20261108.md` sections 5-6, `experiments/rossi_failure_study/tail_channel_diagnosis.py` (cases, how to record the policy), `estimator_swap.py` (`_stiffness_est`, gates), `d_robust.py`, `contact_candidates.py`, `../motion_aware_scan/proposed.py`, `physics_plant3.py` (`replay`, `_tick`, depth equation). No method suggestions. Another Codex run may be using 8 workers; use at most 6.

Known direct mechanism in some failures: k low -> feed-forward or excitation amplified -> delayed loop oscillates. Unknown: why k goes low.

Competing explanations to test, per case, without assuming any:
- H1 Noise: sensor noise alone produces the low window, helped by low signal.
- H2 Window spanning a stiffness change.
- H3 Disturbance in the probe band: motion, or anything else entering z.
- H4 Probe mass or sampling path in the regression.
- H5 Gate or estimator logic accepts what it should not.
- H6 Self-starving probe: an over-estimated k_good shrinks the excitation (force amplitude A*k/k_good), which raises relative sigma and makes low outliers likelier.

## Cases

Adopted Method (nothing attached) on: the 4 cases of `tail_channel_diagnosis.py`; dev cells `('breath_rib','rib','neural_white',28.,1)` and `(...,4)` in `scene_for(..., 'slow_thin')` at base 520000 with seed SEED+q as in `memory_closed_loop.py` stage `dev3`; plus one hard non-failure control, overlay cell `('breath_rib','rib','neural_white',28.,2)` base 674500 (peak 10.2 N). Each `own` run must reproduce its stored peak.

## Evidence per case

**1. First destabilizing event.** Find the first time k/k_good or k/k_bar leaves the band [0.75, 1.3] before the force goes abnormal, and the window that caused it. Tabulate, for that window and the 10 before it: true k and true ln k span in the window; k_hat, sigma, chi; prior k_good; every gate decision (precise, consistent, bound, within, better, the acceptance path taken, k_ref); z = (k_hat - prior)/sigma; probe force amplitude A*k/k_un and relative sigma/k; Y'MY.

**2. Attribution by counterfactual replay.** For the same command sequence up to that window (`P3.replay`, or an equivalent that keeps depth), recompute the window estimate with one factor removed:
- sensor noise off;
- exogenous motion off (add u to the command so that depth is preserved; check depth matches);
- spatial stiffness frozen at its value at window end (rib or field off);
- tool mass off;
- any scripted event off, when the scene has one.
Report k_hat for each and which removal restores it. Check the plain reproduction first: replay with nothing removed must give the original window estimate. Say what this design cannot separate.

**3. Control response.** After the event: magnitudes of feedback, feed-forward and excitation command increments, k/k_good against the feed-forward admissible upper bound 1.23-1.29, oscillation frequency, what reached the force peak.

**4. Statistics of the low tail, on the full existing scene sets.** Run scene A nominal 36 and slow_thin 36 at dev base (no new variants). Over all windows of the adopted Method, as a function of relative sigma (sigma/k) in bins: empirical distribution of (k_hat - c*k)/sigma, where c is the median k_hat/k in that bin; fraction of k_hat < 0.5*c*k against the Gaussian prediction Phi(-0.5*c*k/sigma); the fraction of those low windows the gates accepted into k_good; the size of accepted downward jumps in sigma units. Also relative sigma as a function of k/k_good, to test H6, and as a function of k.

**5. Minimal interventions, on these 7 cases only**, to test causal direction, each changing one thing on the policy instance:
- `kgood_limit`: k_good may fall by at most the same per-update fraction rho_d as k_bar; rising unchanged.
- `exc_from_kbar`: excitation scaled by k_bar instead of k_good.
- `gate_prior_free`: precision gate uses k_ref = k_hat instead of the prior-based k_ref.
- `both`: `kgood_limit` and `exc_from_kbar`.
Report peak and L per case and whether the first destabilizing event still occurs.

## Report

Under 2 pages in Chinese, final message only. One table of causes per case (first event, what restored the estimate, which hypotheses are supported or rejected, with the numbers); the statistics of section 4 in one table; the intervention table. State where cases differ, and what the evidence cannot decide. State only what was measured.

Files: `experiments/rossi_failure_study/low_estimate_forensics.py`, results `results/low_estimate_forensics.json`. `python -X utf8`. Edit no other repo file. No commit, no push.
