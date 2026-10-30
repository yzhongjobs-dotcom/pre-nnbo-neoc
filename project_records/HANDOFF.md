# Handoff: rapid impedance updates under environmental change

Updated 2026-10-30. This file records the current state and what the project does next. The latest user request overrides it. Earlier checkpoint narratives were removed; recover them with `git show b7a3f91:project_records/HANDOFF.md`.

## 1. Start here

- Read this file, then the relevant parts of [long-term requirements](long_term_requirements.md). When entering a research stage, read only the matching file in [research_workflow](research_workflow/). Check [failure review](research_failures_and_corrections.md) before repeating a judgment that failed before.
- Standing user preferences:
  - Report concrete numbers, tables and conclusions in the chat, not only file links.
  - State the quota cost before any large run; do not auto-resume large work.
  - No further simple-controller tuning or scenario expansion.
  - The best-tuned fixed impedance is not a main baseline; fixed impedance uses the method's own lower/upper bounds.
  - The simulation section follows top-paper form: 2–4 scenarios, a few representative comparators, one ablation, computation time.
  - Real-robot experiments come later.
  - New long-term constraints are written only after the user confirms them.

**Current user scope (2026-10-30):** First assess and refine **research status → research question → proposed method**, using both rules and reference layers of `D:/桌面/科研经验总结/SKILL.md`. Do not advance to experiments before these are reviewed. The user supplied C10/C11/S04; all three are now archived and have accepted ARIS original-passage cards. This checkpoint adds scientific assessment and mechanism analysis, not controller changes or experimental evidence.

**Current result (2026-10-30, post-audit research round):** Simulation setup audited (Codex, 19 findings) and fixed; comparators strengthened; scenarios split into static span, moderate motion and a stress test; equal-budget tuning with a robust objective. Research improvement: motion feed-forward now uses a Kalman motion tracker (transfer from manoeuvring-target tracking / radiotherapy latency compensation) fed by force changes converted with the motion-invariant stiffness; two-model IMM was rejected after it destabilised at stiffness edges (motion feed-forward is a disturbance-observer loop whose bandwidth is delay-limited). Final method selected on an independent selection set by a pre-declared rule and reported on an untouched confirmation cohort ([tables](../experiments/motion_aware_scan/results_v2/tables.md)): significantly better than all comparators under tissue motion (incl. privileged-stiffness proxy), better than non-privileged comparators on static span, behind the privileged proxy there. [Paper draft](../research/paper_draft_zh.md) §一–六 updated. Remaining: quantitative bandwidth–delay condition for the motion loop, coverage on fast stiff ramps, real-robot experiments.

## 2. Question and current method

- **Parent question** ([research_question.md](../research/research_question.md)): when the environment changes during continuous robotic scanning of soft tissue, use new contact information to rapidly update impedance and improve actual contact force and motion.
- **Refined question (user-confirmed 2026-10-22):** in constant-force ultrasound scanning of moving, heterogeneous tissue with only contact force and probe position available, tissue motion and local stiffness change are confounded in the contact force. The method must separate them online, using stiffness for impedance adjustment and drift for reference compensation.
- **B0** ([method](../research/method_design.md), [appendix](../research/technical_appendix.txt), [proofs](../research/theoretical_proofs.txt)) remains the revisable starting point. The current C1 is a reduced implementation, not complete B0.
- **C1, current focus.** It runs on the cobot benchmark: [controller](../experiments/cobot_scan_benchmark/contact_update.py), [results](../experiments/cobot_scan_benchmark/contact_update_results/README.md); equations and the local identifiability condition are in [methods §2.7](../research/optimal_impedance_question_and_methods.md). It has three parts:
  - a windowed regression that profiles a cubic force-drift polynomial out of the stiffness estimate;
  - a 128 ms OCP over D, K_PF and u_ff, run with three Gauss-Newton iterations every 8 ms;
  - stiffness-based contact recovery with command-state compensation.

  Admittance K is fixed at 0. "Drift" is an equivalent contact-force drift that also absorbs viscoelastic relaxation; it is not a recovered physical surface position.
- **C2 is set aside at the user's request** ("先不管c2"). It is a residual MLP distilled from a true-stiffness teacher (DAgger, CMA-ES). Keep its artifacts; do not present it as a coequal main claim.
- **Assessment from the 2026-10-26 user conversation ([transcript](../conversation_transcript_2026-10-26.txt) [059]–[065]):**
  - C1's method design already has the technical organization and completeness of related top-journal methods.
  - What is not established is its advantage and contribution: it still lacks fair closest-prior comparisons, the key gain ablation, and real-robot evidence.
  - The three simulation tasks are sufficient.

## 3. Current evidence (seed-7200, 40 scenes per task)

Force RMSE in N, for T1 stiffness change / T2 push–retract soft / T2 stiff / T3 combined scan:

| Method | T1 | T2 soft | T2 stiff | T3 |
|---|---|---|---|---|
| C1 repaired | 0.087 | 0.068 | 0.291 | 0.142 |
| C1 v2 (before repair) | 0.060 | 0.057 | 1.372 | 0.147 |
| C2 | 0.083 | 0.143 | 0.493 | 0.210 |
| EtR-A | 0.107 | 0.147 | 0.551 | 0.230 |
| Haninger (damping-rule adaptation) | 0.122 | 0.178 | 0.433 | 0.273 |
| Jin MS | 0.117 | 0.265 | 1.268 | 0.285 |
| Xue nominal (adaptation) | 0.276 | 0.383 | 3.414 | 0.511 |

- **Effect of the repair on the 20 stiff T2 scenes:**
  - low-contact scenes (force below 0.5 N) fall from 6 to 0;
  - the worst sampled peak falls from 22.23 to 6.78 N;
  - two-event recovery rises from 80% to 100%.

  Restoring the old estimator gives 0.531 N, so the estimator change contributes.
- **Costs:** T1 RMSE worsens and its worst peak rises from 5.33 to 7.58 N. The T3 change is not statistically resolved.
- **Fairness limit:** the C1 repair used T2 development data, while the baselines keep their earlier T3-selected settings. This is not an equally retuned ranking.
- **Ablation (C1 v2, T3):** gains only gives 0.489 N against 0.147 N for the full method. This establishes a benefit of enabling reference compensation over that gains-only instance; it does not isolate the gain-optimization increment with the same reference channel.
- **"Near-optimal" is unsupported:** local OCP reference gaps reach a maximum of 198.96%.
- **Timing:** software P99 is 0.501 / 0.628 / 0.368 ms (T1/T2/T3). Physical communication is excluded, so this is not real-time certification.
- **Sources:** the final Chinese sections are [paper overview](../research/paper_draft_zh.md) and [simulation section](../research/paper_experiments_final_zh.md); closest-prior differences are in [methods §2.8](../research/optimal_impedance_question_and_methods.md).

## 4. Next work, in order

1. N01–N04 done. Keep claims aligned with their originals (pages in field_state §6).
2. Integrate verified points from the independent critique of draft §一–四; keep facts traceable to originals/evidence cards.
3. Design-stage verdict: all baselines met (methods B3 completed afterwards). Simulation plan is in [methods §10](../research/optimal_impedance_question_and_methods.md#10-仿真验证方案2026-10-30待执行) (S1 estimation mechanism, S2 stiffness span, S3 main task, one ablation, timing; fresh test seeds; ~10–12 scenes per condition). After user go-ahead, implement the §三 method changes in the cobot benchmark and run the §四 validation (three existing tasks, Wijayarathne-type / TASE2018-type / DO-MPAC-type / static-surface-MHE comparators, oracle reference, one ablation section, timing), stating quota before large runs.
4. Unproven items to keep explicit: variable-gain and recovery-switching stability, exact scaling under viscoelastic/HC nonlinearity and saturation, physical motion recovery (not claimed), real-robot evidence.

## 5. Where things are

| Need | Location |
|---|---|
| Field state and related-work scope | [field_state.md](../research/field_state.md), [related_work_scope_and_coverage.md](../research/related_work_scope_and_coverage.md) |
| Methods, C1/C2 specification, closest-prior comparison | [optimal_impedance_question_and_methods.md](../research/optimal_impedance_question_and_methods.md) |
| How top papers form methods and run experiments | [reverse analysis](../research/independent_reverse_analysis_of_top_journal_methods.md), [simulation setup analysis](../research/simulation_setup_reverse_analysis.md), [transferable mechanisms](../research/transferable_mechanisms.md) |
| Literature | [literature index](../literature/literature_index.md): 106 archived unique PDFs, 109 usage records; supplied2018 inbox PDF not yet archived |
| Experiments | [experiment index](../experiments/experiment_index.md), [cobot benchmark](../experiments/cobot_scan_benchmark/README.md) |

## 6. Rules and recovery

- Preserve `00_基线快照_2026-10-14/` byte for byte. File check: `python -X utf8 literature/verify_materials.py`.
- Early B0 and 2D experiment directories, plus superseded drafts, were removed from the working tree on 2026-10-26. Recover them from Git `b7a3f91`. Their lessons are in the experiment index and the failure review.
- **ARIS checkpoint2026-10-30:** same run/root/source policy.32 query attempts (27 complete, five historical failures),221 candidate/local/alias records,27 accepted selected-passage cards;106 unique PDFs/109 usage records. Counts are not core citations or coverage certificates.
- C10 is torque-executed explicit force/moment control with moving-phantom evidence; C11 is pre-calibrated geometry/model and image target-force regulation; S04 TAFAC is known-reference velocity/acceleration/control-cycle compensation, online2025-12-06. All are at their appropriate capability layer. C11 formula direction and S04 first-order delay-proof limits are explicit.
- [Fresh independent review](../literature/aris_coverage_review.json) request`13f969fbc87348f8adb511c7fbc5ddea` returned **CONTINUE**, accepted via an actual completed native-child transcript (no Hook claim). It identifies the2018 source gap and does not claim all221 records independently cleared. Controller is **QUERY_PLANNING** with that gap; no formal download batch is currently pending. The exact mirrored candidate map saysSUFFICIENT because that was the reviewed candidate; the adverse review did not accept the landscape.
- C10/C11/S04 download requests are closed. [Current manual checklist](../literature/manual_fulltext_requests.md) requests only the2018 original. No need to redownload completed files or fill JSON. No scientific-core ARIS pipeline was started.
- Assessment derivations explain current C1 only: projection information/error relationship; linearized gain/reference/servo transfer; ideal frozen-gain conditionD>M*k_i. They do not certify actual delayed variable-gain/switching execution. B0 proof text was not independently certified in this checkpoint.
- A prior premature gain diagnostic was stopped before this checkpoint; no new experimental results were added. Its diagnostic files are outside the project working tree, not part of the present evidence. Do not resume it automatically.

- After each completed checkpoint, run `python -X utf8 scripts/research_sync.py --message "research: ..."`. Never force-push.
