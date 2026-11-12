# Git and checkpoint synchronization (read when editing the repository)

- The user authorizes committing and pushing completed changes to `https://github.com/yzhongjobs-dotcom/pre-nnbo-neoc`; do not request permission again for routine completed checkpoints.
- Before editing in a local clone, run `git pull --rebase origin main` and resolve conflicts rather than discarding or overwriting other agents' work. Avoid repeated pulls within one uninterrupted editing operation unless the branch has changed.
- After validating a completed local checkpoint, run:
  `python -X utf8 scripts/research_sync.py --message "research: <specific change and purpose>"`
  The helper handles concurrent commits and baseline-tag recovery. This is for completed checkpoints, not arbitrary file writes.
- If using a direct GitHub API connector to create a commit on the remote branch, verify the current branch head and use a non-forced, concurrency-safe update; the remote commit is already pushed, so do not run a duplicate local sync for that same commit.
- Preserve `00_基线快照_2026-10-14/` exactly. Keep necessary reproducible code, configuration and primary results; record key scientific decisions in `project_records/decision_log.md` and delegated Codex/ARIS prompt plus final output in `project_records/agent_runs/`. Use Git history for change tracking instead of extra change logs.
- Never force-push or rewrite published history. On failure, keep valid local work, explain the blocker, and retry only after checking what changed upstream.
