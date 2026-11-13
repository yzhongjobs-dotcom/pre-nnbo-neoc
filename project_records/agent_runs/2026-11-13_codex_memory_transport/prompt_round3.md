# Task (round 3): the trajectory correction inside the gated default Method

Reuse the runner, scenes and metrics of `experiments/rossi_failure_study/memory_transport_loop.py` / `trajectory_loop.py` / `trajectory_robust.py` (latest copies; if only in your temporary mirrors `D:/pip_env/temp/finite_memory_round2/` or earlier ones, start from there). Windows: `python -X utf8`. **Do not change any existing file, do not commit or push, do not edit research documents, leave `00_基线快照_2026-10-14/` alone, no dynamic arm or other plant option, no tuning, no extra arms.** New code in a new file `gated_delag_loop.py`. If the sandbox refuses repository writes, write to a temporary mirror of the repository paths and say where.

## Question

All trajectory-correction runs so far used the ungated host (`CC.own_mechanism(pol)`), where the estimate collapses in case 7 and fresh scene V4. The default Method keeps its acceptance gates (precision, consistency, feed-forward jump rule, CENTRE). Does the default gated Method already avoid those collapses, and does the stateless trajectory correction (`delag2`) still help on large rib crossings when it sits inside the gated Method?

## Arms

- `gated`: exactly as the existing arms but **without** calling `CC.own_mechanism(pol)` (the policy as built by `B.build('Huber-D4', O.chosen()['Huber-D4'], sc)`), estimator unwrapped. For the 30 original cases its L must equal the stored `Huber-D4` rows of `results/completion_closed_loop.json` (report max |dL|; if not 0, stop and report).
- `gated_delag2`: the same gated policy with the estimator wrapped by the unchanged `delag2` correction (returns `(k_c, sigma*k_c/k_w, chi)`).

Reuse stored rows for `huber` (ungated) and `delag2` (ungated + correction).

## Runs

The 30 original cases (0–23, 36–41) and the 12 fresh A_large scenes, both new arms once each, 8 workers.

## Predictions written before the run (report met / not met with numbers)

- G1. Gated Method avoids the collapses: `gated` L below 0.5 on case 7 and on fresh V4.
- G2. Inside the gated Method the correction still helps large rib crossings: mean L of `gated_delag2` below `gated` by at least 0.05 on original A_large and on fresh A_large; no A_large scene worse than `gated` by more than 0.05.
- G3. Cases 0–23: `gated_delag2` within ±0.01 of `gated`, case by case.
- G4. `gated_delag2` is the best of the four arms in mean L over all 42 scenes.

## Output

`results/gated_delag_loop.json`, one row per run (arm, scene identity, L, rmse, peak, over_15, loss_time, window acceptance rate from `M.metric`'s `accept`). No traces.

Final message (at most 1.5 pages):
1. Per group (A, A_slow_thin, B, A+B, A_large original, A_large fresh, all 42): mean/median L for `huber`, `delag2`, `gated`, `gated_delag2`.
2. Per case for A_slow_thin, A_large original and fresh: L for the four arms, and acceptance rate for the two gated arms.
3. Identification on A_large original fast cases (36, 37, 39, 40) and fresh A_large for the two gated arms, using the stiffness the controller uses for feed-forward (`k_good`: record it per tick from the policy after each step) : median |e| / mean e on plateau, rise, fall, reversal; sustained under-estimate median/worst.
4. G1–G4 verdicts with numbers and the harness check.

State plainly what failed.
