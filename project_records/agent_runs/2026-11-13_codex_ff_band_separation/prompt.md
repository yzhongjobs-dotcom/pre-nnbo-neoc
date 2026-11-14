# Task: closed-loop test of feeding the motion tracker an average of its measurement

Work in `experiments/rossi_failure_study/`. Reuse the runner, scenes, metrics and stored rows of `trajectory_loop.py`, `trajectory_robust.py` (arms `huber`, `delag2`; the 30 original cases, 12 fresh A_large scenes, 12 late push/pull runs). Windows: `python -X utf8`. **Do not change any existing file, do not commit or push, do not edit research documents, leave `00_基线快照_2026-10-14/` alone, no dynamic arm or other plant option, no tuning, no extra arms.** New code in a new file `ff_band_separation_loop.py`. If the sandbox refuses repository writes, write to a temporary mirror of the repository paths and say where.

Theory and predictions: `research/ff_band_separation_20261113.md` (read it first). In short: in the linearised loop model, feeding the motion tracker the mean of the last 6 motion measurements, **without** lengthening its prediction horizon, widens the stiffness error the feed-forward loop tolerates; lengthening the horizon by the 2.5 ticks of delay cancels the effect.

## The change (controller side, estimator untouched)

The motion measurement reaches the tracker through `pol.tracker.update(...)`, which `contact_compare.ModelPolicy._wrap` redirects to the value returned by `ModelPolicy._next_z()` (overridden in `contact_candidates.CandidatePolicy`); the prediction is `self.tracker.predict(L)` in `motion_aware_scan/variants.py::MotionAwareIMM.drift_pred` with `L = (ndel+0.5)*Tc`. Read these to place the change correctly.

Variant `(n, extra)`:
- the tracker receives the mean of the last `n` values that `_next_z()` returned (fewer at the start of a motion run; clear this history whenever the tracker is reset or the motion state is cleared, i.e. whenever the policy restarts the accumulated motion `Z`);
- the prediction uses `predict(L + extra*Tc)`;
- nothing else changes: `pol.Z`, `pol.Zd`, the tracker's parameters and measurement variance, the estimator, the feedback path.

Do this by wrapping on the policy instance after it is built (and re-applying after every `reset`, since the tracker object is recreated). Required check: with `(n, extra) = (1, 0)` the run must reproduce the stored `huber` L exactly on cases 0, 7 and 36 (report max |dL|; it must be 0, otherwise stop and report).

## Block S: stiffness-error sweep (the theory's direct claim)

Inject a scaled true stiffness into the controller: reuse the mechanism of the stored reference `oracle_keff_off` (see `identification_completion.py` / `corrected_comparison.py` for how the true `k_eff` is injected with the Method's gates off) and multiply the injected value by `f`. Cases 0, 1, 3, 36; `f` in (0.5, 0.7, 0.85, 1.0, 1.2, 1.5, 2.0); tracker variants `(1, 0)` and `(6, 0)`. 56 runs. Check that `f = 1`, `(1, 0)` reproduces the stored `oracle_keff_off` L on these cases (report max |dL|).

## Block M: the Method

Ungated Huber-D4 host (`huber`) with tracker variants `(6, 0)`, `(6, 1)`, `(6, 2.5)` on the 30 original cases and the 12 fresh A_large scenes; `delag2` with `(6, 0)` on the same 42 scenes; `huber` with `(6, 0)` on the 12 late push/pull runs. Reuse stored rows for plain `huber` and `delag2`.

## Output

`results/ff_band_separation_loop.json`, one row per run (block, scene identity, arm, n, extra, f, L, rmse, peak, over_15, loss_time). No traces.

Final message (at most 2 pages):
1. Block S: table of L by case, `f` and tracker variant; peak force alongside.
2. Block M, per group (A, A_slow_thin, B, A+B, A_large original, A_large fresh, all 42): mean/median L for `huber`, `huber (6,0)`, `huber (6,1)`, `huber (6,2.5)`, `delag2`, `delag2 (6,0)`. Per case for A_slow_thin, A_large original and fresh. Per group also the mean rmse (N).
3. Late push/pull: per run L and event-window force rms for `huber` and `huber (6,0)`, split toward / away.
4. Verdicts on predictions B1–B5 of the research note with numbers; the two reproduction checks.

State plainly what failed.
