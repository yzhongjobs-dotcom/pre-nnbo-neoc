# Handoff: rapid impedance updates under environmental change

Updated 2026-10-30. This file records the current state and what the project does next. The latest user request overrides it. Earlier checkpoint narratives were removed; recover them with `git show b7a3f91:project_records/HANDOFF.md`.

## 1. Start here

- Read this file, then the relevant parts of [long-term requirements](long_term_requirements.md). When entering a research stage, read only the matching file in [research_workflow](research_workflow/). Check [failure review](research_failures_and_corrections.md) before repeating a judgment that failed before. Key decisions and their reasons are in [decision log](decision_log.md); raw Codex/ARIS prompts and outputs in [agent runs](agent_runs/).
- Standing user preferences:
  - Report concrete numbers, tables and conclusions in the chat, not only file links.
  - State the quota cost before any large run; do not auto-resume large work.
  - No further simple-controller tuning or scenario expansion.
  - The best-tuned fixed impedance is not a main baseline; fixed impedance uses the method's own lower/upper bounds.
  - The simulation section follows top-paper form: 2–4 scenarios, a few representative comparators, one ablation, computation time.
  - Real-robot experiments come later.
  - New long-term constraints are written only after the user confirms them.

**Current user scope (2026-10-30):** First assess and refine **research status → research question → proposed method**, using both rules and reference layers of `D:/桌面/科研经验总结/SKILL.md`. Do not advance to experiments before these are reviewed. The user supplied C10/C11/S04; all three are now archived and have accepted ARIS original-passage cards. This checkpoint adds scientific assessment and mechanism analysis, not controller changes or experimental evidence.

**Related-work update in progress (2026-10-31, user request):** the user questioned that the closest works are old (Rossi 2016, Moreira 2014, Yuen 2010) and asked for 2022–2026 top-venue coverage via the formal ARIS flow, shared with Codex (`gpt-6.1-sol`, medium), done in rounds because of quota; literature is for related work, not all for simulation comparison. Formal run `recent-related-work-20261031` (RUN_ROOT `D:/pip_env/temp/scan_impedance_literature_20261031/aris_run`, stage PAPER_READING): 31 Scholar queries, 347 candidates, 42 admitted, 13 withdrawn at user request, 24 full texts + 5 arXiv. Round 1 read 13 papers (summary in [field state §7](../research/field_state.md); notes in `agent_runs/2026-10-31_read1/`). Key: Wei TMECH2025 assumes a static environment and does not close the control loop; Lei 2025 / Qian ICRA2026 observe motion from US images; Beber ICRA2024 uses an offline viscoelastic map. **Next:** (1) Codex reads remaining route/background papers and converts Claude notes into Evidence Cards; (2) register Wei/Lei/Zakeri/Chen user sources when the Controller allows, read Zakeri and Chen; (3) check cited leads (Li et al. RA-L moving environment + variable admittance; Jiang Actuators2024 IAAC; Tan TASE tumor tracking; Jiang RA-L2021/2022); (4) Field Map + independent coverage review; (5) rewrite draft §1 against skill stage-1 standard. The user will then continue question-by-question review.

**Acceptance round (2026-10-31):** Skill-based horizontal (anchors Pol, K16, Wij, L10, C03) and vertical review; added dose-response sweep (`sweep.py`; best median at every push speed and ratio >=10x, runaway tails at ramps <=0.15 s and ratio >=30x), three figures (`make_figures.py`, research/figures), rewritten title/abstract/contributions; independent Codex re-review: no data or figure errors, two framing defects fixed (coupling arises whenever the controlled probe is also the motion sensor and the gain has error; claims are a tolerance analysis plus empirical risk reduction, 3/50 -> 0/50, not a stability guarantee). Remaining: real robot.

**Current result (2026-10-31, tail-failure round):** Large untouched cohorts (50/family) exposed rare runaways (S3c 5/50 on the diag cohort) that the 6–12-scene sets had hidden. Oracle tests located them in the measurement-gain channel (biased accepted windows, chi2 in [1,3)). Two analysis-derived changes, no new tuned parameters: a property-1 bias bound gating the measurement gain (b <= ln 1.2) and encoder-endpoint displacement (robust to servo-lag/delay mismatch). Final method = tuned_v2/Proposed-IMMb.json. Confirm cohort (seeds 130000+): S3c 0/50 runaways; significantly better than Rossi, DO-MPAC, their combination and a new Yuen-type route-3 comparator (privileged motion sensor) in all five families, and than the true-stiffness Wij proxy in the four motion families (S2 tie). Robust to lag/delay/noise; limits: unknown force filter, faster servo, ridge >40x, S3d long single losses. ARIS survey narrowed novelty (F10, Cortesão/Dominici, Liang IMM exist); 8 full texts requested in literature/manual_fulltext_requests.md. Draft §一–六 rewritten; tables in experiments/motion_aware_scan/results_v2/tables_cohort.md; record in methods §16.

**Previous result (2026-10-30, skill self-check round):** The confirmed method (motion-invariant projection estimator + excitation; force-derived equivalent-motion measurement converted with the unbiased stiffness, IMM tracking and delay prediction as feed-forward; conservative-stiffness normalised admittance; contact supervisor) is unchanged. A usage-B self-check (Claude + independent Codex review) repaired three closure breaks, recorded in [methods §15](../research/optimal_impedance_question_and_methods.md): (1) literature/problem now include predictive physiological-motion compensation (M01–M03, full texts requested in literature/manual_fulltext_requests.md) and difficulty 4 (stiffness is the gain of the motion measurement); (2) new structural property: with exact measurement stiffness the feed-forward loop is open by construction, tolerated k/k_ref ≈ [0.75, 1.2] (common range, Kalman gain recomputed with the controller R) (`ff_loop.py`), which derives the two-stiffness split and FF shut-off at contact loss; rho_d capped to the screened range changes confirm RMSE by −3%..+11% (`consistency_check.py`); (3) matched checks (`mechanism_checks.py`): tracking+prediction is the main gain, IMM vs matched Kalman is condition-dependent (no general switching claim), normalisation mainly helps soft tissue, comparators with the same supervisor do not overturn any result. Second round (user: every component must answer a real problem and grow naturally from status/question): status rewritten around 'where the two kinds of information come from' (three routes + table), question shows failures first and grounds four difficulties; IMM grounded by `tracker_dilemma.py` (equivalent motion includes (y-s) dln k, switches between quiet and abrupt; single-bandwidth two-ends dilemma on soft pushes; IMM mode probability follows push events); measurement-noise-floor refinement tested and not adopted (`tracker_refine.py`). Third round: whole chain re-walked against the skill with M01–M03 full texts read (route 3 premises corrected, Cao TIM2025 added; problem states mechanism cause and dilemma forms; process history moved out of the draft; frozen-stiffness qualification on the FF-loop structure; S3 table shows longest single loss; A1 law and tuning J defined). Draft §一–六 updated. Remaining: real robot; read M01–M03; FF-loop stability under mode switching and time-varying stiffness; residual peaks on stiff pushes. Results: [tables](../experiments/motion_aware_scan/results_v2/tables.md).

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
