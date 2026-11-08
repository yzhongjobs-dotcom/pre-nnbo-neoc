# Task: root-cause analysis of the identification layer on identical observations

Simulation only, corrected code (H increment form is the default). Do not read project_records/HANDOFF.md. Read `research/identification_module_scheme_20261108.md` sections 7-8, then reuse `experiments/rossi_failure_study/low_estimate_forensics.py` (recording, counterfactual replay), `rib_transition_eval.py` (transitions), `memory_segment_eval.py` (segment masks), `memory_filters2.py` (filters).

This is an estimator study, not a closed-loop comparison. All methods see the same recorded measurements.

## Data

Generate each episode once in closed loop with the adopted gated Huber-D4 on corrected code. Record per tick: true |K(20.8)| (col 21), exogenous motion u and breathing, force, encoder, window k_hat, sigma, chi, the Method's k_good and k_bar, gate decisions, probe force amplitude.

Episodes: scene A `breath_rib`, base 874500, seeds 910000 + index as in `corrected_comparison.py`, 6 episodes: nominal i = 2, 6 with (`neural_white`, 28 deg) and i = 6 with (`viscous`, 14 deg); `slow_thin` i = 1, 5 with (`neural_white`, 28 deg) and i = 5 with (`viscous`, 14 deg).

The simulation is small: keep the script simple and do not over-build it. Aim for one run of a few minutes.

## Estimators replayed causally on those data

- gated: the Method's recorded k_good (original Huber-D4).
- raw: the window k_hat with no gates and no memory.
- memory on raw window estimates, all windows, 6-tick updates, R = 4 s^2, parameters as in `rib_transition_eval.py`, no tuning: Kalman tau=2, Kalman tau=0.2, IMM, Lepski, EWA, forgetting with floor, Huang ARD fixed precision.

Replay of `gated` must equal the recorded value, and raw window estimates must recompute to the recorded ones. If not, stop and say so. State clearly that replay leaves out the effect of each estimator on the excitation and on the loop.

Truth: e = ln(k_est/k). The logged k is known to differ from the gain the sampled loop sees (`research/identification_module_scheme_20261108.md` section 4 of `stiffness_module_problem_20261108.md`). Report that gap per episode by measuring the noise-free sampled gain, and treat it as a named part of the systematic error, not as estimator error. Also count the 1.054 centering in `gated` as a named part.

## Step 1. Stable segments (plateau mask of `memory_segment_eval.py`)

Per method and per episode: bias (median e), random part (robust sd of e), RMSE. Compare the sd with the white-noise prediction for that memory. Separate what history can reduce from what it cannot:
- Allan deviation of the raw window estimate ln k_hat on plateaus, against averaging time from 0.2 s up to the plateau length (use overlapping estimator, non-overlapping for the shortest scale): does error keep falling as 1/sqrt(T), and where does it flatten?
- Bias as a function of breathing phase on plateaus: is there a breathing-locked systematic error?
- How much of the bias remains after noise is removed, via noise-off counterfactual on a few plateau windows.
Answer: does history really make up for the short window's lack of information, and how far?

## Step 2. Changing segments

Per transition (rule of `rib_transition_eval.py`, rising and falling; report counts; include the small ones separately since large ones are few): lag, settling time, mean and extreme error during transition with and without the pre-transition offset, largest under-estimate, and share of ticks with k_est/k < 1 and < 0.85. Make one figure: for each method, plateau sd on the x axis against transition lag and under-estimate on the y axis (the noise-versus-tracking trade-off).

## Step 3. Error causes, curve by curve

Figures (PNG, readable, labelled in English, consistent method colours), one per episode, each with panels sharing a time axis: (a) true k(t) and k_hat(t) for raw, gated and each memory method, transitions shaded; (b) error ln(k_est/k) for the same; (c) breathing u(t), force and probe force amplitude. A separate zoom figure for the two largest transitions and the two largest under-estimates. Write them to `research/figures/id_layer/`.

Attribute window error by counterfactual replay on the same commands, in a stated sequential order, for windows every 0.2 s on plateaus and every tick in transitions:
- noise: full minus noise-off;
- disturbance: then motion off with depth preserved, as in `low_estimate_forensics.py`;
- k change within the window: then stiffness frozen at the window-end value;
- remainder: systematic (includes sampling, probe mass, model).
Excitation: relate error size to probe force amplitude and relative sigma.
Memory lag: difference between each filter's error and the concurrent raw window error, split by plateau and transition.
Report the share of each cause by segment, per method where it applies, and state that the order affects the split, by giving at least one alternative order on a subsample.

## Report

In Chinese, final message only, under 3 pages. Replay checks. Tables for steps 1, 2 and cause shares for step 3. For each method, one line: accuracy on plateau, accuracy and lag in transition, and why, from the evidence. Say what the evidence cannot decide. No method suggestions, no winner selection, report what was measured.

Files: `experiments/rossi_failure_study/id_layer_root_cause.py`, results `results/id_layer_root_cause.json`, figures as above. Up to 8 workers, `python -X utf8`. No commit, no push. Edit no existing file.
