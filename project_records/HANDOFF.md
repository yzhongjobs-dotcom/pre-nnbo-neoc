# Handoff: rapid impedance updates under environmental change

Updated 2026-10-29. This file records the current state and what the project does next. The latest user request overrides it. Earlier checkpoint narratives were removed; recover them with `git show b7a3f91:project_records/HANDOFF.md`.

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

**Latest user scope and stop (2026-10-29):** Continue ARIS supplementation, self-correct drift, admit only scientifically relevant research status, then stop for acceptance. Available-material delivery is ready: ten valid new originals, three reused primary representatives, 20 accepted partial-text cards; main prose updated. Do not automatically resume experiments, ideas or further searches. Two original/timing issues remain (invalid DOB-MPC HTML and missing ACT-US), with explicit coverage limits; see the manual checklist.

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
- **Ablation (C1 v2, T3):** gains only gives 0.489 N against 0.147 N for the full method. Most of the benefit comes from drift estimation plus reference compensation, not yet from impedance optimization.
- **"Near-optimal" is unsupported:** local OCP reference gaps reach a maximum of 198.96%.
- **Timing:** software P99 is 0.501 / 0.628 / 0.368 ms (T1/T2/T3). Physical communication is excluded, so this is not real-time certification.
- **Sources:** the final Chinese sections are [paper overview](../research/paper_draft_zh.md) and [simulation section](../research/paper_experiments_final_zh.md); closest-prior differences are in [methods §2.8](../research/optimal_impedance_question_and_methods.md).

## 4. Next work, in order

1. **Fair closest-prior comparison** (main set in methods §2.8):
   - Xue (full adaptation);
   - Wang TIE2026 G02;
   - DOSMAC (lumped DOB plus sliding-mode admittance; use the standard form if the full text is unavailable, and say so);
   - Cao-style vision compensation, using noisy, delayed true surface motion.

   Fixed high/low impedance only shows the classical tradeoff. Give every method the same tuning data.
2. **Rule-gain ablation:** the same estimator and reference optimization, but gains set by a rule. This tests whether gain optimization helps beyond reference compensation. If it does not, decide with the user whether K should enter the optimization or the contribution should be reframed.
3. **Identifiability derivation:** stiffness versus equivalent drift under constant-force regulation, with and without dither, and why a lumped disturbance observer cannot schedule gains.
4. **Stability/passivity condition for the variable gains:** account for K = 0, force integral and delayed servo. Kronander H04 does not apply directly; see the energy-tank treatment in Fu H01.
5. **"Near-optimal":** support it with a converged reference in a declared action class, or drop it. Wording changes are the user's call.
6. **Real-robot readiness** (transcript [065]):
   - calibrate the measurement chain (tool gravity, force zero, timestamps, filter and execution delays);
   - fit the estimator to measured noise (robust loss, a position measurement-error model if needed);
   - gate updates by reliability: hold stiffness when excitation is insufficient, fade drift compensation when unreliable, use hysteresis for contact switching;
   - check the actual excitation signal-to-noise ratio;
   - then rerun the three tasks with measured noise, and only then run physical experiments (heterogeneous phantom plus motion stage).

   Platform availability is unconfirmed; ask the user.
7. **Literature gaps:**
   - DOSMAC and Zhang MHE-VIC now have the recorded original-passage evidence; remaining access and coverage issues are in section6;
   - the DOB plus model-predictive admittance ultrasound lead (January2027 issue; earliest online date pending);
   - joint-estimation citation chains.

   Available-material supplementation is delivered with20 selected-passage Evidence Cards; current remaining originals and coverage limits are in section6. Wait for user acceptance; no comparison or simulation is started.
8. **Update** the paper overview and the simulation section after items 1–5.

**Open user decisions:**
- the question wording ("near-optimal"; stiffness versus D/K_PF);
- target tier (RA-L/TMECH or T-RO);
- the eventual fate of C2.

**Not to do:**
- C2 optimization;
- new scenarios;
- extra simple-controller tuning;
- reframing the question to protect C1.

## 5. Where things are

| Need | Location |
|---|---|
| Field state and related-work scope | [field_state.md](../research/field_state.md), [related_work_scope_and_coverage.md](../research/related_work_scope_and_coverage.md) |
| Methods, C1/C2 specification, closest-prior comparison | [optimal_impedance_question_and_methods.md](../research/optimal_impedance_question_and_methods.md) |
| How top papers form methods and run experiments | [reverse analysis](../research/independent_reverse_analysis_of_top_journal_methods.md), [simulation setup analysis](../research/simulation_setup_reverse_analysis.md), [transferable mechanisms](../research/transferable_mechanisms.md) |
| Literature | [literature index](../literature/literature_index.md): 101 unique PDFs, 104 usage records |
| Experiments | [experiment index](../experiments/experiment_index.md), [cobot benchmark](../experiments/cobot_scan_benchmark/README.md) |

## 6. Rules and recovery

- Preserve `00_基线快照_2026-10-14/` byte for byte. File check: `python -X utf8 literature/verify_materials.py`.
- Early B0 and 2D experiment directories, plus superseded drafts, were removed from the working tree on 2026-10-26. Recover them from Git `b7a3f91`. Their lessons are in the experiment index and the failure review.
- **ARIS delivery for user acceptance (2026-10-29):**
  - run `c1-related-work-20261026`, root `D:/pip_env/temp/scan_impedance_literature_20261026/aris_run`; approved policy/cutoff unchanged, canonical SHA256 `b7d45f09459c9ed5d5d9b36278fe1cce280d8a76da068e148e25fc00814f91c7`.
  - current stage **HUMAN_SEARCH_REQUIRED / human_fulltext_batch**. Original eight-paper batch completed. Current formal request is ACT-US only; the [consolidated checklist](../literature/manual_fulltext_requests.md) also identifies the invalid DOB-MPC file. Do not claim completion without these direct-neighbor/timing gaps.
  - 11 user files: ten valid PDFs retained, one Elsevier security-verification HTML not admitted as paper Evidence. Library is 101 unique PDFs /104 usage records. Existing PDF/text bytes and B0 snapshot unchanged.
  - 25 query attempts (20 completed, five preserved historical adapter failures), 158 candidate/local/alias records; 103 OUT_OF_SCOPE, two DUPLICATE, two HOLD_IDENTITY, 51 IN_SCOPE. Do not use these as core-paper counts. Latest backbone refinement returned two below-threshold diagnostic abstracts; it does not certify high-citation coverage.
  - 20 accepted selected-passage cards: earlier E32/G22/E57/H01/G02/G13/G21; new C01–C08/S01–S02; reused E27/G23/E41. [Field state](../research/field_state.md) incorporates task/control capabilities and conditions; reviews/resources do not automatically enter the main comparison.
  - DOSMAC S01 now has original equations: lumped stiffness/position compensation in zero-K admittance; hardware combined variation comes from fixed curved surfaces, simulation uses true nominal position. MHE-VIC S02 assumes a premeasured/static normal surface; MHE estimates Ke/De/robot disturbance and rules update gains. E27/G23 learning interfaces and E41 OS/MS optimization are source-grounded, not algorithm-name overlap judgments.
  - independent fresh-context [review](../literature/aris_coverage_review.json) returned **CONTINUE**, accepted via the actual native-child transcript. Registry omissions for learned/predictive representatives were corrected afterward; no second sufficient-coverage review claimed. The [formal map](../literature/aris_active_field_map.md) remains the last17-card PARTIAL version, while3 newer cards await allowed map revision after the pending batch. Historical search and unavailable decisive evidence remain explicit limits.
  - [Search/decisions](../literature/aris_literature_search.json), [scope and self-corrections](../research/related_work_scope_and_coverage.md), [current request](../literature/aris_human_fulltext_request.json), [manifest template](../literature/aris_fulltext_manifest_template.json). Do not replay old download requests or start from scratch.
  - **Stop for user acceptance now.** If the user supplies the remaining originals and asks to resume: verify originals/date, register and read the current ACT-US batch, inspect DOB-MPC, then revise the same map/account and review the remaining material coverage gaps. No automatic C2 work, tuning, simulations or method selection.

- After each completed checkpoint, run `python -X utf8 scripts/research_sync.py --message "research: ..."`. Never force-push.
