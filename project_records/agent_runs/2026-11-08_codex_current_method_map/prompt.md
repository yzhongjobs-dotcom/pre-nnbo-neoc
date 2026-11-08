# Task: map the CURRENT default Method as actually implemented (read-only)

Repository: E:/github_project/scan-impedance-update. Read code only. Do NOT edit any file, run experiments, commit or push. Do not touch `00_基线快照_2026-10-14/`.

Context: the current default Method ("Huber-current", called Huber-D4 in `experiments/rossi_failure_study/baseline_interpretation.md` and `overlay_interpretation.md`) is the configuration used as "当前 Method" in `experiments/rossi_failure_study/baseline_compare.py` / `overlay_compare.py` (selected config in `experiments/rossi_failure_study/results/baseline_test_r2.json` or `overlay_plan.json`). The paper draft `research/paper_draft_zh.md` §三 describes an OLDER version (projection + IMM); the code has evolved. I need the exact current pipeline.

Trace, from `overlay_compare.py`/`baseline_compare.py` down through `method3.py`, `sweeps.py`, `estimators.py`, `contact_candidates.py`, `motion_aware_scan/common.py` and whatever they import, the per-control-cycle (8 ms) data flow of the current default Method. Report for each stage: what it computes, the exact equation (as implemented), inputs/outputs, parameter values in the default config, and file:line.

Stages to cover (merge/rename if the code differs):
1. Sensing pre-processing (ForceSensing: inertia compensation, delay alignment, filter inversion, how z and Y are built; is Y complex/phase-led? what is the "H" contact model: gain + phase lead, how is the phase obtained online).
2. Excitation: frequency, amplitude, how converted to position (which k used).
3. Window estimator: exact objective of Huber-current (k and d; D4 prior; Huber parameter; window length; rho/fc; solver/iterations), outputs (k_hat, sigma_k, chi, d estimate?).
4. Gates and trust: sigma/chi rules, acceptance (b ≤ ln1.2 still used?), FF_JUMP rule, prior gating (GATE_ON_PRIOR), conservative k_bar update (c, rho_d, clipping).
5. Feedback controller: admittance A1 form and parameters, how scaled by k_bar, integrator or not.
6. Feed-forward: what signal is predicted (eta from force/k̂ minus displacement? the estimated d? something else), which tracker/predictor (IMM? Kalman? parameters), prediction horizon, how it enters y_cmd, when it is disabled.
7. Contact loss/recovery supervisor.
8. Anything else that affects the command (limits, saturations, rate limits).
Also state: total computation per cycle if recorded somewhere (file/result), and list which parameters were selected on dev data.

Deliverable: write a concise report (max ~2.5 pages markdown, equations in LaTeX-ish plain text) to the -o output file. Include a one-paragraph "difference from paper_draft_zh.md §三" section. No process narration, no long code dumps; quote at most 3 lines of code per stage.
