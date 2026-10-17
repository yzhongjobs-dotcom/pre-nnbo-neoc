# Research work and GitHub synchronization

This is a long-term research project. The user authorizes committing and pushing
completed project changes to https://github.com/yzhongjobs-dotcom/pre-nnbo-neoc.

1. Read `项目记录/HANDOFF.md` to resume, then the current user request and
   relevant parts of `项目记录/长期要求.md` and `项目记录/研究决策与变更.md`.
   Older root handoffs and dated plans were consolidated; recover originals
   from Git only when needed. Do not repeatedly reread all history.
2. Current scope: the user has started independent theory/experiment review and
   authorized a quick evaluation of candidate 2 (learned adaptive impedance).
   This authorization does not select candidate 2 as the main method or change
   the core research question. See the current HANDOFF and candidate report.
3. Follow the user-approved scientific principles in `项目记录/长期要求.md`:
   co-refine the question, prior limitations and method from field understanding,
   domain taste and required capabilities. Preserve the field-state research
   module while updating its evidence. Do not shrink the question to preserve an
   algorithm, evade prior because modules overlap, or require NEOC to be uniquely
   indispensable. Failure analysis and mechanism transfer are useful routes,
   not mandatory formation paths for every contribution. For a substantially
   changed core question, present grounded options for user selection; explore
   within the selected question proactively without repeated approval.
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
   documents now live under `研究/`; B0 contents and equations were preserved.
   They may change in future based on research evidence, with explicit deltas.
6. Keep useful reproducible source/configuration/primary results together. Files
   retained merely to prove that work happened may be removed from the working
   tree after recording necessary lessons and a Git recovery point. Do not use
   archival completeness as a reason to preserve irrelevant material locally.
7. For each completed logical change, update `项目记录/研究决策与变更.md` with
   its reason, comparison version, evidence, scope and affected assumptions or
   proofs; perform appropriate checks, then automatically run:
   `python -X utf8 scripts/research_sync.py --message "research: <specific change and purpose>"`
   Run before reporting completion and at the end of a turn that changes files.
   No repeat upload permission is needed. If pushing fails, keep the local
   commit, report the failure and retry. Never force-push or rewrite history.

The synchronization helper runs at completed checkpoints, not arbitrary file
writes. Preserve the user's notification and review preferences; report actual
changes plainly, without making the user read inventory dumps.
