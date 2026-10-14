# Research work and GitHub synchronization

The user authorizes committing and pushing this project's baseline and every
subsequent completed research change to:
https://github.com/zyanzhi90-dot/pre-nnbo-neoc

1. Start from `00_基线与变更记录.md` and the user's latest research objective.
   Read `01_项目推进复盘与后续研究依据_2026-10-14.md` for the completed
   historical audit, selective reuse and the corrected next research workflow.
   Preserve the reference snapshot under `00_基线快照_2026-10-14/`.
2. Develop the idea from the field state, strongest-prior failure mechanisms and
   domain taste. Existing modules may be reused. Do not shrink the question to
   preserve a named algorithm or treat module overlap as a reason to evade prior.
3. After each completed logical research change, update the change record with
   its reason, evidence, baseline/version, affected assumptions/proofs and scope.
   This also applies to diagnostic or negative results; do not call them a
   completed new method or a new discovery when already present in the baseline.
4. Perform checks appropriate to the change, then automatically run:
   `python -X utf8 scripts/research_sync.py --message "research: <specific change and purpose>"`
   Run this before reporting the change complete and at the end of a turn that
   changed project files. No repeat permission is needed for this authorized sync.
5. The script commits project changes and pushes main. If pushing fails, retain
   the local commit, report the failure, and retry on the next checkpoint.
   Never force-push, reset user work, rewrite history, or hide a failed upload.
6. Ignore local caches and credentials. Keep research reports, source scripts,
   data, logs and proof changes together in the corresponding research commit.

Automatic synchronization is performed at completed research checkpoints, not
halfway through a file write or while presenting an unfinished result as complete.
