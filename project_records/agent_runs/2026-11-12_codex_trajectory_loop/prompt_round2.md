# Task (round 2): closed-loop test of a stateless, bounded trajectory correction

Continue from `experiments/rossi_failure_study/trajectory_loop.py` and `results/trajectory_loop.json` (round 1, already in the repository). Windows: `python -X utf8`. **Do not change any other existing file, do not commit or push, do not edit research documents, leave `00_基线快照_2026-10-14/` alone, no dynamic arm or other plant option.** You may extend `trajectory_loop.py` (add the arm and a `--round2` entry that reuses the stored `huber` and `traj` rows instead of rerunning them). If the sandbox refuses to write inside the repository, write to a temporary mirror and say where.

## Why round 2

Round 1's filter kept its own level state. In closed loop the first, short window gave a level several times too high, the over-estimate shrank the probe motion (excitation is divided by the stiffness), the window lost its information and the filter never recovered (cases 1, 11, 14, 15, 20 from t = 0.09 s). The correction itself worked on the rib crossings. Round 2 keeps the level entirely on the memoryless window and uses the trajectory only for the change across the window's lag, bounded, with a gain that starts at zero and is learned from the data.

## New arm `delag`

Wrap the estimator as in round 1. At each call with finite `k_w` (and `sigma`):

Histories kept by the wrapper (reset them all if `n` does not advance by one):
- `y[j] = o.log[n-1, 2]`, `F[j] = o.log[n-1, 5]`.
- `y6`, `F6`: means of the last 6 values; `y24`, `F24`: means of the last 24 values (available prefix at start).
- `C_w[j] = F24/k_w`, `rel[j] = sigma/k_w`; store `C_w`, `y24`, `rel` per tick (NaN when `k_w` is not finite).

Correction, using the gain learned up to the previous tick:

    x = y6 - y24
    r = clip(1 + beta*x/C_w, 0.5, 2.0)
    k_c = clip(k_w * (F6/F24) / r, 200, 3e5)
    return (k_c, sigma*k_c/k_w, chi)

Skip the correction (return the original triple) when `F24 < 1` or fewer than 24 ticks of history exist.

Gain update, after the correction, when ticks `j` and `j-24` both have finite `C_w` and `rel <= 1` and every tick between them had `F >= 1`:

    dx = y24[j] - y24[j-24]
    dC = C_w[j] - C_w[j-24]
    v  = C_w[j]^2*rel[j]^2 + C_w[j-24]^2*rel[j-24]^2      # noise variance of dC
    lam = exp(-Tc/4.0)
    Sxx = lam*Sxx + dx*dx
    Sxy = lam*Sxy + dx*dC
    Sxn = lam*Sxn + dx*dx*v
    beta = clip((Sxy - sqrt(24*Sxn)) / (Sxx + 2e-5), 0, 1)

`Sxx = Sxy = Sxn = 0`, `beta = 0` at start. (The subtraction is one standard deviation of the noise in `Sxy`; windows 24 ticks apart overlap 24-fold in the sum.)

No tuning of these constants.

## Predictions written before the run (report met / not met with numbers)

- Q1. A_large (36–41): mean L at most 0.22; no case worse than `huber` by more than 0.02.
- Q2. A, A_slow_thin, B, A+B (0–23): no case worse than `huber` by more than 0.05; each group mean within ±0.02 of `huber`, except that A_slow_thin may improve through case 7 (report it either way).
- Q3. A_large fast cases (36, 37, 39, 40): mean error of the returned stiffness on rise ticks within ±0.06, on fall ticks below +0.12.
- Q4. Learned gain: on A_large fast cases `beta` is above 0.05 for most ticks after t = 3 s; on B cases its median is below 0.1.

## Runs and output

Arm `delag` on cases 0–23 and 36–41, 8 workers. Reuse the stored round-1 rows for `huber` and `traj`. Save traces as in round 1 to `results/trajectory_loop/case_{index:02d}_delag.npz`, with `k_w`, `k_ret`, and `beta` per estimator tick. Add the rows to `results/trajectory_loop.json`.

Final message (at most 2 pages), same table layout as round 1 with `delag` beside `huber` and `traj`:
1. Per case: L for `huber`, `traj`, `delag`, stored `Huber-D4_imm`, `oracle_keff_off`; peak force for `huber` and `delag`.
2. Per group: mean/median L per arm; force rms on crossings and on plateau.
3. Per group: identification table (median |e| / mean e per segment; sustained under-estimate median/worst), A_large split fast/slow.
4. Per case: median and final `beta` after t = 3 s.
5. Q1–Q4 verdicts with numbers.

State plainly what failed; do not tune or add analyses to rescue a prediction.
