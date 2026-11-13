# Task (round 2): finite memory with every past window transported directly to now

Continue from `experiments/rossi_failure_study/memory_transport_loop.py` (round 1; if it is only in the temporary mirror `D:/pip_env/temp/memory_transport_mirror_99ou0lkl/`, start from there) and reuse its runner, scenes, metrics and stored rows. Windows: `python -X utf8`. **Do not change any other existing file, do not commit or push, do not edit research documents, leave `00_基线快照_2026-10-14/` alone, no dynamic arm or other plant option, no tuning, no extra arms.** If the sandbox refuses repository writes, write to a temporary mirror of the repository paths and say where.

## Why

Round 1's recursive memory (`kf`) removed the collapse (case 7, fresh V4) but lost the rib-crossing benefit. Trace of case 37: from 1.6 to 2.0 s the state was carried up by the per-tick transport to 18 times the truth while the window stayed near 1; the per-tick measurement weight was only about 2%, and the transport uncertainty was added in quadrature per tick although the gain error is common to all ticks. The over-estimate then shrank the probe and the window became uninformative. Round 2 removes the recursion: each stored window output is transported to the present by its own bounded factor, so nothing accumulates.

## New arms

Wrap the Huber-D4 window estimator as before (ungated host). Keep the `delag2` histories and gain learner unchanged (`y6`, `y24`, `F24`, `C_w`, `rel = sigma/k_w`, `beta`). Store per tick `ln k_w`, `C_w`, `y24`, `rel` (NaN when `k_w` is not finite).

At call `j`, with `beta` as learned up to the previous tick, use the stored ticks `i = j-H … j` that have finite `k_w`, `rel_i <= 3`, and lie in the current unbroken contact run (no reset between `i` and `j`):

    d_i = ln( clip(1 + beta*(y6[j] - y24[i]) / C_w[i], 0.5, 2.0) )
    v_i = rel_i^2 + (0.5*d_i)^2
    w_i = 1/v_i
    lam = sum_i w_i*(ln k_w[i] - d_i) / sum_i w_i
    n_used = number of ticks used;  span = j - (oldest tick used);  n_eff = (span + 24)/24
    P = (1/sum_i w_i) * n_used/n_eff
    return ( clip(exp(lam), 200, 3e5), k*sqrt(P), chi )

If no tick qualifies, return the original triple. Then update `beta` as in `delag2`. With `H = 0` this is `delag2`; check that on one run (returned values equal to the stored `delag2` trace within 1e-12).

Arms:
- `tmem24`: `H = 24`.
- `tmem48`: `H = 48` (a planned second length, not a tuning sweep).
- `mem24`: `H = 24` with `beta` forced to 0 (finite memory without transport).

## Runs

Same 42 scenes as round 1 (original cases 0–23, 36–41; the 12 fresh A_large scenes), each arm once, 8 workers. Reuse stored rows for `huber`, `delag2`, `kf`.

## Predictions written before the run (report met / not met with numbers)

- M1. `tmem24` keeps the rib-crossing benefit: original A_large mean L at most 0.195; fresh A_large mean at most 0.265; no A_large scene (original or fresh) worse than `delag2` by more than 0.05.
- M2. `tmem24` removes the collapse: case 7 and fresh V4 both below 0.5.
- M3. Original cases 0–23 except 7: `tmem24` within ±0.02 of `huber`, case by case.
- M4. Transport matters: over original and fresh A_large together, mean L of `mem24` is at least 0.03 above `tmem24`.
- M5. No lock-up: no run of the three new arms has `k_ret/k_eff > 3` for more than 0.5 s in total after t = 0.5 s.
- M6. Identification: pooled original plateau median |e| at least 15% below `huber` for `tmem24` and at least 25% for `tmem48`; on the original fast A_large cases (36, 37, 39, 40) the mean error on fall ticks is at most +0.15 for `tmem24`.

## Output

Add rows to `results/memory_transport_loop.json`. Save traces for cases 7, 36, 37 and fresh N4, N5, V4, V5 (new arms) to `results/memory_transport/` with `t`, `fax`, `k_eff`, `t_est`, `k_w`, `k_ret`, `P`, `beta`.

Final message (at most 2 pages), same layout as round 1 with columns `huber`, `delag2`, `kf`, `tmem24`, `tmem48`, `mem24`:
1. Per group mean/median L.
2. Per case L for A_slow_thin, A_large original and fresh.
3. Identification per group: median |e| / mean e on plateau, rise, fall, reversal; sustained under-estimate median/worst.
4. Ripple and crossing force rms per group.
5. M1–M6 verdicts with numbers, the M5 statistic of the worst run per arm, and the `H = 0` check.

State plainly what failed.
