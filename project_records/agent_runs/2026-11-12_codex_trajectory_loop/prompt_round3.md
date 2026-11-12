# Task (round 3): robustness of the stateless trajectory correction (`delag`)

Continue from `experiments/rossi_failure_study/trajectory_loop.py` (arms `huber` and `delag`, round 2). Windows: `python -X utf8`. **Do not change any other existing file, do not commit or push, do not edit research documents, leave `00_基线快照_2026-10-14/` alone, no dynamic arm or other plant option, no change to the `delag` constants.** Put the new code in a new file `experiments/rossi_failure_study/trajectory_robust.py` that imports the arms from `trajectory_loop.py`. If the sandbox refuses to write inside the repository, write to a temporary mirror and say where.

## Scientific question

Round 2 gave one closed loop per case: on the six A_large cases mean L fell from 0.253 to 0.189 and the other 24 cases were unchanged. Is that improvement stable, and what happens when a push/pull arrives after the gain has been learned?

Predictions written before the run (report met / not met with numbers):
- R1. Fresh A_large scenes: mean L of `delag` below `huber`; no scene worse by more than 0.02; at least 9 of 12 improved or equal (within 0.005).
- R2. New sensor seeds on the six original A_large scenes: in every (scene, seed) `delag` is not worse than `huber` by more than 0.02, and the mean difference is negative.
- R3. ±0.5% perturbation: the A_large mean of `delag` stays below the mean of `huber` for every perturbation; for cases 36 and 37 the spread of `delag` L across perturbations is smaller than the `huber`−`delag` gap of round 2.
- R4 (expected weakness). Push/pull after the gain is learned: `delag` is worse than `huber` around the event. Bound expected: mean L increase at most 0.05, no peak force above 15 N. Report the pull-away direction separately.
- R5. Fresh A scenes: |ΔL| at most 0.01 in every scene and `beta` stays 0.

## Blocks (all with arms `huber` and `delag`, same harness as round 2, 8 workers)

1. **Fresh A_large scenes.** `CC.cases()` uses overlay cells with index `i` in (0, 2, 6) for A_large. Build the same kind of case for `cond = 'breath_rib'`, variant `'large_rib'`, `i` in (1, 3, 4, 5, 7, 8), both settings in `CC.SETTINGS`: 12 scenes. Sensor seed: `CC.SENSOR_BASE + 1000 + running number`. Check and report that the cells exist and that the scene builder accepted them.
2. **Fresh A scenes.** Same, variant `'nominal'`: 12 scenes, seeds `CC.SENSOR_BASE + 2000 + running number`.
3. **New sensor seeds.** The six original A_large cases (36–41) with seeds `case['seed'] + 10000` and `+ 20000`.
4. **Perturbation.** Cases 36–41, the value returned to the controller multiplied by `(1 + eps)`, `eps` in (+0.005, −0.005), applied after each arm's own processing (for `delag`, after the correction). Original seeds.
5. **Late push/pull on large ribs.** `cond = 'breath_rib_fast'`, variant `'large_rib'`, `i` in (0, 2, 6), both settings: 6 scenes, seeds `CC.SENSOR_BASE + 3000 + running number`. In `scenes3.py` the push/pull starts at `0.2*T` (2 s), before the gain is learned. After building the scene dict, move the event to start at t = 6.0 s **in your script only**: rebuild the event wave from `event_info` with the same formula as `scenes3.py` (`amp*(_smooth((t-start)/dur) - _smooth((t-start-dur-hold)/dur))`, `dur = 1.875*amp/peak_limit`, `hold = 0.15`), find its sign from the stored `u`, subtract it from `sc['u']` and add it back shifted to 6.0 s; update `events` and `event_info`. Verify and report: the reconstructed wave matches the original contribution to `u` (max abs difference of `u_original - wave - (breathing + cardiac part)` cannot be checked directly, so instead check that `u_original - wave` is smooth across the original event interval: report its max second difference there against the same measure elsewhere), and that nothing else in the scene depends on the event time (search `scenes3.py`; if muscle activation `act` or anything else is tied to this event kind, say so and handle it consistently). Run each scene twice: with the original sign, and with the event wave's sign flipped. 12 runs per arm. Label each run "toward probe" or "away from probe" by whether the contact force rises or falls at the event onset in the `huber` run.

## Output

`results/trajectory_robust.json` with one row per run (block, scene identity, seed, arm, eps or sign, L, rmse, peak, over_15, loss_time). For blocks 1 and 5 also save traces as in round 2 to `results/trajectory_robust/`.

Final message (at most 2 pages):
1. Block 1: per scene L `huber` / `delag`, median `beta` after 3 s; means and medians; identification error on rise, fall, reversal (median |e| / mean e, pooled over fast scenes — call a scene fast when its scan speed is above 15 mm/s) for both arms; force rms on crossings.
2. Block 2: per scene L both arms, max `beta`.
3. Block 3: per (scene, seed) L both arms.
4. Block 4: per case L for eps = 0 (round 2), +0.005, −0.005, both arms.
5. Block 5: per run L both arms, peak force both arms, `beta` just before the event, and during the event window (start to end + 0.6 s): mean and minimum of `ln(k_ret/k_eff)` and force rms, both arms; split by direction.
6. R1–R5 verdicts with numbers.

State plainly what failed; no tuning, no extra analyses.
