# Rapid impedance updates under environmental change

Research goal: during continuous robotic scanning of soft tissue, use new contact information to generate near-optimal impedance promptly, with actual force and motion-control benefits. The target paper requires a substantive scientific contribution and evidence sufficient to compete at leading journals; a complete method and evidence at that level have not yet been established.

Read the [handoff](project_records/HANDOFF.md) first for the current state and next work. Confirmed goals and stage standards are in the [long-term requirements](project_records/long_term_requirements.md); stage-by-stage guidance is in [research_workflow](project_records/research_workflow/). The current method is C1 on the [cobot scan benchmark](experiments/cobot_scan_benchmark/README.md): it repairs stiff push/retract contact loss on three simulated tasks, but fair closest-prior comparisons, the gain ablation, a stability condition and real-robot evidence are still missing.

The research proceeds iteratively through field analysis, question refinement,
method design, experimental verification, failure review and long-term improvement.
Use the entries below to connect each claim to its method, evidence and limitations.

| Entry | Purpose |
|---|---|
| [Long-term requirements](project_records/long_term_requirements.md) | Confirmed goals, workflow, experiment order and reporting requirements |
| [Research question](research/research_question.md) | Core task, B0 starting point and unverified premises |
| [Field state](research/field_state.md) | Existing capabilities, conditions and full-text access boundaries |
| [Two candidate methods](research/optimal_impedance_question_and_methods.md) | Correct candidate definitions and distinction between specific experiments and complete concepts |
| [Method design](research/method_design.md) | Complete B0 model, assumptions, decisions and execution relationships |
| [Transferable mechanisms](research/transferable_mechanisms.md) | Reviewed mature mechanisms and adaptation conditions; no newly selected final method |
| [Literature index](literature/literature_index.md) | Local originals, text, sources and reading scope from two external repositories and existing project literature |
| [Experiment index](experiments/experiment_index.md) | Reusable models, solvers, comparators, results and reproduction entries |
| [Research failures and corrections](project_records/research_failures_and_corrections.md) | Causes, consequences and corrections of actual project deviations |
| [Independent analysis of leading-journal methods](research/independent_reverse_analysis_of_top_journal_methods.md) | How sample papers organize questions, methods and evidence; complements the failure review |

The B0 [method](research/method_design.md), [appendix](research/technical_appendix.txt) and [conditional proofs](research/theoretical_proofs.txt) remain a revisable research starting point. Preserve the frozen snapshot `00_基线快照_2026-10-14/` unchanged. Git retains complete older material; do not create another historical-copy directory.

Repository: [yzhongjobs-dotcom/pre-nnbo-neoc](https://github.com/yzhongjobs-dotcom/pre-nnbo-neoc).
