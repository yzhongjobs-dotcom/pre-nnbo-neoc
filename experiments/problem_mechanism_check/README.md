# Problem-stage mechanism check

Minimal linear analysis supporting the research-question argument in
[paper draft §一–二](../../research/paper_draft_zh.md). It is a failure
demonstration for the two closest routes, not a benchmark of the proposed method.

Plant: the cobot benchmark's normal axis (8 ms cycle, one-cycle command delay,
19 ms first-order servo, 30 Hz force filter, 0.08 N noise, 4 N target), linear
spring contact `f = k (y - s)`, admittance with K = 0 and the C1 rule-prior
scaling (M, D proportional to the design stiffness, K_PF = 2/k_design).

Run: `python -X utf8 gain_and_bias.py > results.txt` and `python -X utf8 loop_margin.py > margin_results.txt` (numpy only, deterministic).

Scope: the loop is the C1 rule-prior admittance, not DOSMAC/DO-MPAC; it bounds that loop only and is not evidence about the ESO-based methods.

| Check | Result |
|---|---|
| Stiffness-normalized loop vs ratio κ = k / k_design, nonlinear simulation (`gain_and_bias.py`) | settles for κ ≤ 1; oscillating at κ = 2 (13 N peak-to-peak), divergent from κ ≈ 4 |
| Same loop, exact linear stability (`loop_margin.py`, discrete closed-loop eigenvalues, no contact loss/saturation) | stable for κ ∈ (0, 1.19]; force bandwidth ≈ 33 / 8.3 / 2.2 Hz at κ = 1 / 0.3 / 0.033 |
| Static-surface windowed regression (intercept only) while 1 mm 0.25 Hz breathing + 5 mm push is regulated at 4 N | k̂/k median 0.003–0.18, i.e. severe underestimation, with or without 20.8 Hz excitation |
| Cubic-drift projection, no excitation | median 0.97–1.07, 128 ms spread [0.85, 1.46] |
| Cubic-drift projection with 0.02 mm 20.8 Hz excitation | median 1.00, 128 ms spread [0.93, 1.07] |

A first script version paired each force sample with the next position sample
and showed a 0.5 bias with excitation (60° phase at 20.8 Hz); the fixed version
pairs the position that produced the force. This is a reminder that encoder/force
time alignment matters for any excitation-based estimator.
