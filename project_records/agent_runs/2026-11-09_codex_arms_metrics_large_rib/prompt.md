# Task: three method arms, core tracking metrics, one large rib-contrast scene

Repository: E:/github_project/scan-impedance-update. Work in `experiments/rossi_failure_study/`.
Read `identification_completion.py`, `corrected_comparison.py`, `rib_transition_eval.py`, `memory_closed_loop.py:scene_for`, `memory_segment_eval.py:gates_off` first. Use `python -X utf8`.

Constraints: no commit or push. Do not edit research documents, HANDOFF or decision logs. Leave `00_基线快照_2026-10-14/` alone. No retuning, no new estimators, no change to moduli, `k_rib`, thresholds or the frozen configurations (`O.chosen()`). No extra experiments beyond those listed. If something does not work, report it; do not tune around it.

## 1. Method arms

Baselines stay as they are (`Rossi-ID`, `Kalman`, `MHE`, `Xue-const`, `Simple-1`, `Simple-2`, full controllers).

Our methods get three arms, on `Huber-D4` and `D4`:

| id | meaning |
|---|---|
| `Huber-D4`, `D4` | gated, unchanged |
| `Huber-D4_nogate`, `D4_nogate` | `CC.own_mechanism(pol)`, the same definition the baselines use |
| `Huber-D4_imm`, `D4_imm` | `CC.own_mechanism(pol)` plus `R.make_filter(pol, 'imm')` |

`P5` keeps its gated row only.

Reason for the clean definition: `off_none` and the memory arms use `SE.gates_off` alone, which leaves CENTRE on and leaves `rel_max`, `chi2_max` acting in the phase detector. Evidence in `results/identification_replay.json`: `off_none` k_good minus raw pre-transition bias is about ln 1.054. So these existing arms are not gate-free.

In replay, the other memory forms (`kal_long`, `kal_short`, `lepski`, `ewa`, `forget_time`, `ard_fixed`) stay as replay-only rows on Huber-D4 and also use the clean definition. `off_none` is replaced by `Huber-D4_nogate`; `imm` by `Huber-D4_imm`. In `results/corrected_comparison.json` add a note that the existing `off_none` and memory closed-loop rows carry CENTRE and the phase-detector sigma gate and are superseded. Do not delete them.

## 2. One new scene, `A_large`

Append to `CC.cases()`, after the existing 36, so their seeds do not change:
`('A_large', 'breath_rib', (0, 2, 6), 'large_rib')`, both `SETTINGS`, 6 cases. Update the count assertion to 42. Append `'A_large'` to `SCENES`.

In `scene_for`, variant `large_rib`: `B.S3.cell_scene3(c, t, 'nominal', i, base=base, overrides=dict(slide='short', force=8., dynamics=dy, phase=ph, cover=.018))`. Nothing else changes: same rib cells, field and motion as scene A, because only the probe's orientation relative to the ribs and the cover differ. `slide='short'` means the probe's long axis lies parallel to the ribs and it moves across them, so the footprint does not span two ribs.

Existing 36 scenes must stay byte-identical. Check this with one scene hash per group and by reproducing one stored closed-loop row per method group.

## 3. Evidence script `rib_contrast_evidence.py`, output `results/rib_contrast_evidence.json`

Port `D:/pip_env/temp/claude/8fbcb2674095e6be/e4c4f610-c229-451e-9777-043bc8fe9776/scratchpad/rib_contrast_check.py` and `rib_contrast_check2.py`. Keep the Fig. 4C values as written there, labelled "read by eye from Al-Zogbi 2021 Fig. 4C". Produce:

1. Fig. 4C adjacent ratios and overall max/min per curve.
2. Static, ribs only (`field_sigma=0`, `surface_amp=0`): per-crossing ratio of the 20 N contact depth, and of 8 N tangent, for `slide` long and short at cover 7, 18, 32 mm.
3. Per-crossing ratio of logged k_eff (`log[:,26]`, .1 s boxcar, 3% prominence) for A, A_slow_thin, A+B, A_large.
4. Sizes of scored no-event transitions per group.
5. Soft-region check: same scene with ribs removed, |G| at 8 N, long versus short. Report the ratio. Do not adjust anything if it differs.

## 4. Metrics (Table 3)

In `rib_transition_eval.score`, add per channel `change_mae`: mean |ln(k_hat/k_eff)| over valid ticks in [t0, t1] only.

Table 3 columns, for each of raw and k_good, rising and falling:
方法 | 方向 | 事件数 | 变化段误差中位 | 最大低估中位 | 最大低估最差 | 恢复时间中位 | 已恢复 | 未恢复 | 回合结束截断

- 变化段误差: median across events of `change_mae`.
- 最大低估: existing `raw_under_magnitude`, window [t0, t1+.5 s], never negative. Not debiased.
- 恢复时间: existing `recovery_raw`, within ln 1.2 held .2 s, from the end of the change, censored at 2 s, Kaplan–Meier median.
- Remove lag and the under-fraction and over-estimate columns from the tables. Keep lag and the others in the JSON. One sentence under the table: lag is kept in the data but is not a main metric, because a curve shift has no direct response meaning.

Families, separate subtables, all no-event: existing scenes |Δln k| >= .3; existing scenes .15 <= |Δln k| < .3; `A_large` |Δln k| >= .3. `A_large` is excluded from the two existing-scene families. Event-coupled stays in the JSON only.

Rows: Huber-D4 three arms, D4 three arms, P5, Rossi-ID, Kalman, MHE, Simple-1, Simple-2. Xue stays out of Tables 2 and 3.

## 5. Tables 1 and 2

Add the same arms, and `A_large` as a new group column in Table 1 and a new group in Table 2. Table 1 keeps Xue-const, full controllers and the diagnostic Oracle rows. Label arms in Chinese: 门控, 无门控, 无门控＋IMM. Regenerate `comparison_tables.md` from saved data. Update its header sentence: A_large = same scene as A with the probe parallel to the ribs and an 18 mm cover. Memory other than IMM is not in the tables; say so in one line.

## 6. Runs

- Closed loop: the 4 new arms on all 42 cases; every method already in Table 1, including `oracle_keff_off`, `oracle_keff_idealcal`, `Huber-D4_idealcal`, on the 6 `A_large` cases. Reuse stored rows for the old 36 cases. Add to `results/completion_closed_loop.json`.
- Replay: sources Huber-D4 and Rossi-ID on 42 cases, per-tick data as now in `results/id_replay_ticks/`. Then score and tables.
- Keep the replay self-check: the source method's replayed commands must match the log, with the same error bound as now.

## 7. Report (`output.md`, at most 2 pages)

Conclusion first. Then: files changed; the reproduction checks; Table 1 rows for Huber-D4 three arms and D4 three arms, with A_large; the evidence numbers in section 3; how many A_large transitions were scored; any A_large case that loses contact or diverges with every method including Oracle. Report failures as failures.
