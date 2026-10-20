# Research work and GitHub synchronization

This is a long-term research project. The user authorizes committing and pushing
completed project changes to https://github.com/yzhongjobs-dotcom/pre-nnbo-neoc.

1. Read `project_records/HANDOFF.md` to resume, then the current user request and
   relevant parts of `project_records/long_term_requirements.md`.
   Older root handoffs and dated plans were consolidated; recover originals
   from Git only when needed. Do not repeatedly reread all history.
2. Research context (not an instruction to start research on every turn): the
   user explicitly corrected the drift toward abstract whole-acquisition
   candidate A. For research work, continue from the original question: use online
   environment information to rapidly update near-optimal impedance with actual
   control value. Compare the model/optimization route and the learned impedance
   route; investigate task benefit and computational bottlenecks together.
   Candidate A and its tests remain contextual evidence, not the current selected
   replacement question. Any substantive problem revision must state its basis
   and impact clearly. See HANDOFF for the latest evidence and next work.
3. Follow `project_records/long_term_requirements.md` for the approved goal, research workflow
   and stage-specific judgment standards. Co-refine important questions, prior
   limitations and complete methods through interleaved literature analysis,
   design and quick task-relevant verification. Actively propose and implement
   candidate mechanisms; do not turn the workflow into a defensive checklist or
   require all prior to fail before exploration. Deepen strong in-field work and
   actively transfer mature mechanisms where useful. For a substantially changed
   core question, explain grounded options for user selection; explore within
   the selected question proactively without repeated approval.
4. Distinguish baseline, candidate changes and tested instances. Negative local
   results only change the supported layer; accumulated evidence may justify
   higher-level revision. Do not use locality to shield an idea indefinitely.
   Do not call known baseline facts new
   discoveries, nodal roots optimal controllers, or prototype timings real-time
   certification. Check original evidence when it affects a decision.
   The complete B0 question and method remain the current refinement starting
   point; they have not been invalidated by local tests. Do not silently replace
   them with the pre-B0 September 8 question or a new algorithmic subproblem.
   Co-refinement may adjust an inaccurate question without requiring strict
   falsification, provided the reasons support that level of change.
5. Preserve `00_基线快照_2026-10-14/` byte-for-byte. The three working core
   documents now live under `research/`; B0 contents and equations were preserved.
   They may change in future based on research evidence, with explicit deltas.
6. Keep useful reproducible source/configuration/primary results together. Files
   retained merely to prove that work happened may be removed from the working
   tree after recording necessary lessons and a Git recovery point. Do not use
   archival completeness as a reason to preserve irrelevant material locally.
7. Use Git history for change tracking. Do not create separate change logs,
   cleanup records, language-migration reports or relocation histories.
   Keep scientific reasoning in the relevant research document, experimental
   report or failure analysis. After checking each completed checkpoint,
   automatically run `python -X utf8 scripts/research_sync.py --message "research: <specific change and purpose>"`.
   No repeat upload permission is needed. If pushing fails, keep the local commit,
   report the failure and retry. Never force-push or rewrite history.

The synchronization helper runs at completed checkpoints, not arbitrary file
writes. Preserve the user's notification and review preferences; report actual
changes plainly, without making the user read inventory dumps.

## Language and paths

Keep only these working documents in Chinese:
- `research/research_question.md`
- `research/field_state.md`
- `research/optimal_impedance_question_and_methods.md`
- `research/method_design.md`
- `project_records/long_term_requirements.md`
- `project_records/research_failures_and_corrections.md`

Use English for all other editable documentation, reports, indexes, handoffs,
technical appendices, proofs, comments, and operational instructions. Working
paths use English. Preserve original papers, extracted source texts, raw data,
the frozen snapshot, and literal historical identifiers in their original form.
Language edits must preserve scientific content, equations,
values, conditions, references and evidence boundaries. Small corrections to
previously identified wording or reference defects are allowed.
Communicate with the user in Chinese. Explain what each step did, how it was
done, and its results when reporting progress or answering questions. Maintain
one working language version per document. On Windows use `python -X utf8`.
