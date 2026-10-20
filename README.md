# Rapid impedance updates under environmental change

Research goal: during continuous robotic scanning of soft tissue, use new contact information to produce timely near-optimal impedance with actual force and motion control benefits. The intended paper must offer a substantive scientific contribution and evidence strong enough to compete at a leading journal. A complete method and evidence package meeting that goal has not yet been established.

Start with the [handoff](project_records/HANDOFF.md). The user's confirmed understanding of method formation and contribution scale is recorded in the [long-term requirements](project_records/long_term_requirements.md), with paper-specific evidence in the [reverse analysis](research/independent_reverse_analysis_of_top_journal_methods.md). The latest [dual-candidate task simulations](experiments/dual_candidate_task_simulation/README.md) implement both candidates under the original question, with [24 scenarios, ablations, distribution shifts and release-delay experiments](experiments/dual_candidate_task_simulation/results.md), followed by an explicit task adaptation of Jin's multistep predictive variable impedance method. Force and tangential errors still trail that comparator. A journal-level contribution has not yet been established for either candidate. Candidate 1 uses local quadratic value feedback, not the B0 fixed-collocation root.

| Entry | Purpose |
|---|---|
| [Long-term requirements](project_records/long_term_requirements.md) | Confirmed goals, workflow, experiment order and reporting requirements |
| [Research question](research/research_question.md) | Core task, B0 starting point and unverified premises |
| [Field state](research/field_state.md) | Prior capabilities, conditions and full-text access limits |
| [Two candidate methods](research/optimal_impedance_question_and_methods.md) | Correct candidate definitions and the distinction between implementations and complete concepts |
| [Transferable mechanisms](research/transferable_mechanisms.md) | Examined mature mechanisms and adaptation conditions, not selected new designs |
| [Literature index](literature/literature_index.md) | Local papers/texts, sources and reading scope for the two external repositories and original project literature |
| [Experiment index](experiments/experiment_index.md) | Reusable models, solvers, comparators, results and reproduction entry points |
| [Failures and corrections](project_records/research_failures_and_corrections.md) | Actual project mistakes, consequences and corrections |
| [Independent reverse analysis](research/independent_reverse_analysis_of_top_journal_methods.md) | How sample papers organize questions, methods and evidence; complementary to the failure review |
| [Decisions and changes](project_records/decisions_and_changes.md) | Necessary history, recovery versions and consolidation rationale |

The B0 [method](research/method_design.md), [appendix](research/technical_appendix.txt) and [conditional proofs](research/theoretical_proofs.txt) remain a revisable research starting point. The frozen snapshot `00_基线快照_2026-10-14/` remains unchanged. Git retains complete older materials; do not create a second historical archive.

Operational documentation and working paths use English. Scientific prose and experiment explanations remain Chinese for the user's reading. The frozen snapshot, original papers/data, historical provenance and pinned Git URLs retain their original form. See the [migration record](project_records/english_migration.json) for old-to-current paths. Claude Code loads the shared [agent instructions](AGENTS.md) through [CLAUDE.md](CLAUDE.md); both agents use the same working files.

Audit the migrated checkpoint with `python -X utf8 scripts/verify_english_migration.py`. This checks paths, source changes, links and preserved evidence without overwriting historical experiment manifests. Old experiment hash checks describe their original Git versions; do not regenerate those records merely to make renamed sources pass. The migration audit is specific to this checkpoint, so future intentional scientific edits need their own recorded checks.

Repository: [yzhongjobs-dotcom/pre-nnbo-neoc](https://github.com/yzhongjobs-dotcom/pre-nnbo-neoc). After completing authorized changes, check and record them under [AGENTS.md](AGENTS.md), then run:

```text
python -X utf8 scripts/research_sync.py --message "research: specific change and purpose"
```
