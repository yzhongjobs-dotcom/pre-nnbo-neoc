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
(pre-declared selection on seeds 80000+: single-model Kalman chosen) and `run_v2.py`.
Primary results: `results_v2/confirm/` (seeds 90000+, untouched) and `results_v2/tables.md`
(`make_tables_v2.py`). `results_v2/test` and `results_v2/replication` are the earlier round with
the IMM tracker as final method (it destabilised in 2/12 T3 scenes); `results/` and
`results_replication/` are the pre-audit round. Diagnostics: `diag_phases.py`, `test_imm.py`.

```
FINAL=Proposed-KF SPLIT=confirm python -X utf8 run_v2.py && python -X utf8 scaling_check.py && python -X utf8 make_tables_v2.py
```
