# Task (round 4): remove the force-ratio factor from `delag` and re-validate

Continue from round 3 (`trajectory_robust.py`, `trajectory_loop.py`; if they are only in your temporary mirror `D:/pip_env/temp/trajectory_robust_round3/`, start from there). Windows: `python -X utf8`. **Do not change any other existing file, do not commit or push, do not edit research documents, leave `00_基线快照_2026-10-14/` alone, no dynamic arm or other plant option, no tuning.** If the sandbox refuses repository writes, write to a temporary mirror and say where.

## Why

Round 3, block 2, scene V4 (fresh nominal A scene, viscous 14°, index 4): L went from 0.029 (`huber`) to 0.190 (`delag`) although the learned gain stayed 0. With `beta = 0` the only difference from `huber` is the factor `F6/F24` in `k_c = k_w*(F6/F24)/r`. That factor says the stiffness is proportional to the force whenever the trajectory carries no information, which is not physical (for a power-law contact the exponent is about 0.2–0.3, not 1), and it feeds force fluctuations between about 5 and 20 Hz straight into the gain. Hypothesis: this factor caused V4.

## Work

1. **Diagnose V4.** Rerun block-2 scene V4 with `delag` saving the trace. Report: the time the two arms' forces first differ by more than 0.2 N, the range and rms of `F6/F24` before and after that time, the dominant frequency of the force error after that time, peak force, and whether `beta` was 0 throughout.

2. **New arm `delag2`.** Identical to `delag` except

       k_c = clip(k_w / r, 200, 3e5),   r = clip(1 + beta*x/C_w, 0.5, 2.0)

   (no `F6/F24`). Gain learning unchanged. With `beta = 0` it must return `k_w` exactly; check this on one run by comparing the returned values with `k_w` bit for bit while `beta` is 0.

3. **Runs with `delag2`** (reuse stored `huber` and `delag` rows): the 30 original cases (0–23, 36–41); round-3 blocks 1, 2, 3 and 5 (same scenes, seeds and event shifts).

Predictions written before the run (report met / not met with numbers):
- S1. Block 2 (fresh A scenes), `delag2`: |ΔL| versus `huber` at most 0.01 in every scene, including V4.
- S2. `delag2` keeps the rib-crossing benefit: original A_large mean L within ±0.01 of `delag` (0.189); block 1 mean within ±0.02 of `delag` (0.260); no A_large or block-1 scene worse than `huber` by more than 0.02.
- S3. Original cases 0–23: |ΔL| versus `huber` at most 0.005 in every case.
- S4. Block 5 (late push/pull): mean ΔL versus `huber` at most +0.02; report the away-from-probe runs separately.

## Output

Add rows to `results/trajectory_robust.json` (and the original-case rows to `results/trajectory_loop.json` under arm `delag2`). Final message (at most 2 pages):
1. V4 diagnosis.
2. Original 30 cases: per group mean/median L for `huber`, `delag`, `delag2`; per case for A_large.
3. Block 1: per scene L `huber` / `delag` / `delag2`, median `beta` (delag2) after 3 s; identification error on rise, fall, reversal (median |e| / mean e) for the three arms; crossing force rms.
4. Block 2: per scene L three arms. Block 3: per (scene, seed) L three arms.
5. Block 5: per run L three arms, peak force, `beta` before the event, event-window force rms; split by direction.
6. S1–S4 verdicts with numbers.

State plainly what failed; no tuning, no extra analyses.
