# Motion-aware contact estimation and stiffness-normalised admittance (simulation study)

Implements and evaluates the method of [paper draft §三](../../research/paper_draft_zh.md) against
literature comparators. Plant: `../cobot_scan_benchmark/plant.py` (direct-command mode added:
policies return an absolute position command each 8 ms; same delay, servo lag, speed limit and
force filter for every method; the controller reads the measurement at the start of the period).

| File | Role |
|---|---|
| `common.py` | direct-command base class, window regression (numba core), metrics |
| `design_A1.py` | offline unit-stiffness admittance design; frozen stability and circle-criterion test |
| `A1_feasible.json` | 535 circle-certified designs (κ ∈ [0.5, 2.14], k_i = 0) |
| `proposed.py` | proposed controller, true-state reference controller, fixed-impedance ends |
| `comparators.py` | Rossi-type, Wijayarathne-type (privileged stiffness), DO-MPAC, Rossi+DO-MPAC (first draft by Codex, reviewed; Rossi excitation set to the paper's fixed displacement amplitude) |
| `tune_equal.py` | equal-budget tuning on development seeds (41 random + 20 local per method) |
| `tuned/*.json` | selected parameters and all tuning evaluations |
| `run_experiments.py` | S1 estimation, S2 stiffness span, S3 main task, ablation, timing |
| `scaling_check.py` | re-optimisation per stiffness vs scaled design |
| `make_tables.py` | builds `results/tables.md` |
| `results/`, `results_replication/` | main test cohort (seeds 9500–9711, 9800+) and fresh replication cohort (10500–10711) |

Reproduce (local CPU, about 3 minutes on 16 cores after numba compilation):

```
python -X utf8 design_A1.py
python -X utf8 tune_equal.py && TUNE_STAGE=2 python -X utf8 tune_equal.py
python -X utf8 run_experiments.py && python -X utf8 scaling_check.py
REPLICATE=1 python -X utf8 run_experiments.py
python -X utf8 make_tables.py
```

Development history that matters for interpreting results: the first test run used a faulty
reference controller and a pure-Python estimator; both were fixed (the proposed method's
mathematics unchanged) and the same test seeds were re-run, then confirmed on a fresh cohort.
An information-driven excitation and a sample-based drift feed-forward were tried on development
seeds only and rejected. Results are simulation evidence, not robot or imaging validation.

## Post-audit study (current results)

After an independent setup audit (19 findings; fixes listed in paper draft §5.7) the study was
re-run with `scenes_v2.py` families (S2 static span, S3a–c moderate motion, S3d stress test),
`tune_v2.py` (equal 61-evaluation budget, objective mean log(1+J)), `space_patch.py` (fixed
comparators), `imm.py`/`variants.py` (Kalman/IMM motion trackers), `select_final.py`
(pre-declared selection on seeds 80000+: IMM chosen after the 1 ms alignment fix; the misaligned round had chosen single-model Kalman, archived in results_v2/pre_alignment_fix/) and `run_v2.py`.
Primary results: `results_v2/confirm/` (seeds 90000+, not used for tuning or selection; evaluated once before the alignment fix) and `results_v2/tables.md`
(`make_tables_v2.py`). `results_v2/test` and `results_v2/replication` are an earlier round (before the 1 ms
alignment fix) with the IMM tracker as final method (it destabilised in 2/12 T3 scenes); `results/` and
`results_replication/` are the pre-audit round. Diagnostics: `diag_phases.py`, `test_imm.py`.
Self-check checks: `ff_loop.py` (frozen-gain analysis of the force-closed motion feed-forward loop -> results_v2/ff_loop.json), `consistency_check.py` (rho_d capped to the A1 screening range, supervisor removed / reduced to the FF gate -> results_v2/consistency.json), `mechanism_checks.py` (comparators with the same supervisor, predictor swap with all else fixed, normalisation removed -> results_v2/mechanism_checks.json), `tracker_dilemma.py` (single-bandwidth two-ends test and IMM mode probabilities, quiet vs event segments -> results_v2/tracker_dilemma.json), `tracker_refine.py` (displacement-domain measurement-noise floor, pre-declared adoption rule; not adopted -> results_v2/tracker_refine.json).

```
FINAL=Proposed-IMM SPLIT=confirm python -X utf8 run_v2.py && FINAL=Proposed-IMM python -X utf8 scaling_check.py && python -X utf8 make_tables_v2.py
```
