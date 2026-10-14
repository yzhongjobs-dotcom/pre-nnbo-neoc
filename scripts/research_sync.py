"""Commit and push an authorized completed research checkpoint."""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

EXPECTED_REMOTE = "https://github.com/zyanzhi90-dot/pre-nnbo-neoc.git"


def git(root: Path, *args: str, timeout: int = 120) -> str:
    env = dict(os.environ, GIT_TERMINAL_PROMPT="0", GCM_INTERACTIVE="never")
    run = subprocess.run(["git", "-C", str(root), *args], capture_output=True,
                         text=True, encoding="utf-8", errors="replace",
                         env=env, timeout=timeout)
    if run.returncode:
        raise RuntimeError(f"git {args[0]} failed: {run.stderr.strip() or run.stdout.strip()}")
    return run.stdout.strip()


@contextmanager
def sync_lock(root: Path):
    gitdir = Path(git(root, "rev-parse", "--absolute-git-dir"))
    lock = gitdir / "research-sync.lock"
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        raise RuntimeError("Another sync is running, or a previous interrupted sync left its lock. Inspect before retrying.")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump({"pid": os.getpid(), "started_utc": datetime.now(timezone.utc).isoformat()}, handle)
        yield gitdir
    finally:
        lock.unlink(missing_ok=True)


def verify_snapshot(root: Path):
    manifest = root / "00_基线快照_2026-10-14" / "manifest.json"
    data = json.loads(manifest.read_text(encoding="utf-8"))
    for entry in data["files"]:
        snapshot = (root / entry["snapshot"]).resolve()
        if not snapshot.is_relative_to(manifest.parent.resolve()):
            raise RuntimeError("Baseline snapshot path escaped its directory")
        if hashlib.sha256(snapshot.read_bytes()).hexdigest() != entry["sha256"]:
            raise RuntimeError(f"Reference baseline snapshot changed: {entry['snapshot']}")
    # Future working method edits are allowed; frozen reference edits are not.
    changed = git(root, "diff", "baseline-b0-2026-10-14", "--", "00_基线快照_2026-10-14")
    if changed:
        raise RuntimeError("Reference baseline or its manifest changed; preserve the B0 snapshot")


def sync(root: Path, message: str, expected_remote: str = EXPECTED_REMOTE,
         check_snapshot: bool = True) -> dict:
    root = root.resolve()
    if Path(git(root, "rev-parse", "--show-toplevel")).resolve() != root:
        raise RuntimeError("Sync must run against the project repository root")
    if git(root, "remote", "get-url", "origin") != expected_remote:
        raise RuntimeError("Origin does not match the authorized repository")
    if git(root, "branch", "--show-current") != "main":
        raise RuntimeError("Automatic sync is restricted to main; inspect the current branch")
    with sync_lock(root) as gitdir:
        for marker in ["MERGE_HEAD", "CHERRY_PICK_HEAD", "REVERT_HEAD", "rebase-merge", "rebase-apply"]:
            if (gitdir / marker).exists():
                raise RuntimeError("An unfinished Git operation requires inspection")
        if (gitdir / "index.lock").exists():
            raise RuntimeError("Git index is locked by another operation")
        if git(root, "ls-files", "-u"):
            raise RuntimeError("Unresolved merge conflicts require inspection")
        if check_snapshot:
            verify_snapshot(root)
        before = git(root, "rev-parse", "HEAD")
        git(root, "add", "--all", "--", ".")
        staged = git(root, "diff", "--cached", "--name-only", "-z")
        if staged:
            git(root, "commit", "-m", message)
        commit = git(root, "rev-parse", "HEAD")
        # Retry previously unpushed commits even when no new files changed.
        git(root, "push", "origin", "HEAD:refs/heads/main", timeout=300)
        remote = git(root, "ls-remote", "origin", "refs/heads/main").split()[0]
        if remote != commit:
            raise RuntimeError("Remote main does not match the local checkpoint")
        return {"created_commit": commit != before, "commit": commit,
                "changed_files": len([p for p in staged.split("\0") if p]),
                "pushed": True, "remote_matches": True}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message", required=True)
    args = parser.parse_args()
    if not args.message.strip():
        parser.error("A meaningful checkpoint message is required")
    try:
        result = sync(Path(__file__).resolve().parents[1], args.message)
        print(json.dumps(result, ensure_ascii=False))
        return 0
    except (RuntimeError, OSError, ValueError, subprocess.TimeoutExpired) as exc:
        print(json.dumps({"pushed": False, "error": str(exc),
                          "note": "Existing local commits are retained; do not force-push."},
                         ensure_ascii=False), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
