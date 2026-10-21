"""Commit and push an authorized completed research checkpoint."""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time
from urllib.parse import urlsplit

EXPECTED_REMOTE = "https://github.com/yzhongjobs-dotcom/pre-nnbo-neoc.git"
BASELINE_TAG = "baseline-b0-2026-10-14"
# Transient transport failures worth retrying (not authentication or path errors).
NETWORK_ERRORS = ("Could not resolve host", "Connection reset", "Connection refused",
                  "Connection timed out", "Operation timed out", "RPC failed",
                  "early EOF", "unexpected disconnect", "HTTP 429", "HTTP 502",
                  "HTTP 503", "HTTP 504", "TLS connection", "Failed to connect")
REJECTED = ("non-fast-forward", "fetch first", "[rejected]", "Updates were rejected")


def git(root: Path, *args: str, timeout: int = 120) -> str:
    env = dict(os.environ, GIT_TERMINAL_PROMPT="0", GCM_INTERACTIVE="never")
    run = subprocess.run(["git", "-C", str(root), *args], capture_output=True,
                         text=True, encoding="utf-8", errors="replace",
                         env=env, timeout=timeout)
    if run.returncode:
        raise RuntimeError(f"git {args[0]} failed: {run.stderr.strip() or run.stdout.strip()}")
    return run.stdout.strip()


def _repo_identity(url: str):
    """(host, owner/repo) for URL remotes, None for local paths.

    Accepts the same repository written as https with or without `.git`,
    with a trailing slash, embedded credentials, scp-style ssh, or the
    `/git/<owner>/<repo>` path used by sandboxed git proxies (cloud agents).
    """
    url = url.strip()
    scp = re.match(r"^[\w.-]+@([\w.-]+):(.+)$", url)
    if scp:
        host, path = scp.group(1), scp.group(2)
    elif "://" in url:
        parts = urlsplit(url)
        host, path = (parts.hostname or ""), parts.path
    else:
        return None
    segs = [seg for seg in path.strip("/").split("/") if seg]
    if len(segs) < 2:
        return None
    repo = segs[-1][:-4] if segs[-1].lower().endswith(".git") else segs[-1]
    return host.lower(), f"{segs[-2]}/{repo}".lower()


def remote_matches(actual: str, expected: str) -> bool:
    a, e = _repo_identity(actual), _repo_identity(expected)
    if a is None or e is None:
        if a is not None or e is not None:
            return False
        return Path(actual).resolve() == Path(expected).resolve()
    if a[1] != e[1]:
        return False
    proxy = a[0] in ("127.0.0.1", "localhost") and "/git/" in actual
    return a[0] == e[0] or proxy


def git_network(root: Path, *args: str, timeout: int = 300, attempts: int = 5) -> str:
    """Run a fetch/push, retrying only transient network failures (2, 4, 8, 16 s)."""
    for attempt in range(attempts):
        try:
            return git(root, *args, timeout=timeout)
        except RuntimeError as exc:
            if attempt == attempts - 1 or not any(m in str(exc) for m in NETWORK_ERRORS):
                raise
            time.sleep(2 ** (attempt + 1))
    raise AssertionError("unreachable")


def ensure_baseline_tag(root: Path):
    """Shallow or tagless clones (e.g. cloud agents) fetch the baseline tag once."""
    try:
        git(root, "rev-parse", "--verify", "--quiet", f"refs/tags/{BASELINE_TAG}^{{commit}}")
    except RuntimeError:
        git_network(root, "fetch", "--no-tags", "origin",
                    f"+refs/tags/{BASELINE_TAG}:refs/tags/{BASELINE_TAG}")


def integrate_remote(root: Path):
    """Bring in commits other agents pushed, replaying local checkpoints on top.

    Never discards or force-overwrites remote work: a conflicting rebase is
    aborted and reported, leaving the local commits intact for inspection.
    """
    git_network(root, "fetch", "--no-tags", "origin", "+refs/heads/main:refs/remotes/origin/main")
    head, upstream = git(root, "rev-parse", "HEAD"), git(root, "rev-parse", "refs/remotes/origin/main")
    if head == upstream:
        return
    try:
        git(root, "merge-base", "--is-ancestor", upstream, "HEAD")
        return  # only local commits are ahead; nothing to integrate
    except RuntimeError:
        pass
    try:
        git(root, "rebase", "refs/remotes/origin/main")
    except RuntimeError as exc:
        try:
            git(root, "rebase", "--abort")
        except RuntimeError:
            pass
        raise RuntimeError("Remote main has conflicting changes; local commits kept, "
                           f"resolve manually before syncing ({exc})")


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
        # The frozen manifest was produced on Windows; accept its separators
        # without rewriting the byte-preserved reference manifest on POSIX.
        snapshot = (root / entry["snapshot"].replace("\\", "/")).resolve()
        if not snapshot.is_relative_to(manifest.parent.resolve()):
            raise RuntimeError("Baseline snapshot path escaped its directory")
        if hashlib.sha256(snapshot.read_bytes()).hexdigest() != entry["sha256"]:
            raise RuntimeError(f"Reference baseline snapshot changed: {entry['snapshot']}")
    # Future working method edits are allowed; frozen reference edits are not.
    ensure_baseline_tag(root)
    changed = git(root, "diff", BASELINE_TAG, "--", "00_基线快照_2026-10-14")
    if changed:
        raise RuntimeError("Reference baseline or its manifest changed; preserve the B0 snapshot")


def sync(root: Path, message: str, expected_remote: str = EXPECTED_REMOTE,
         check_snapshot: bool = True) -> dict:
    root = root.resolve()
    if Path(git(root, "rev-parse", "--show-toplevel")).resolve() != root:
        raise RuntimeError("Sync must run against the project repository root")
    if not remote_matches(git(root, "remote", "get-url", "origin"), expected_remote):
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
        # Retry previously unpushed commits even when no new files changed.
        # Another agent (Codex or Claude Code) may have pushed meanwhile:
        # integrate and retry a rejected push a few times, never force-push.
        for attempt in range(3):
            integrate_remote(root)
            if check_snapshot:
                verify_snapshot(root)
            try:
                git_network(root, "push", "origin", "HEAD:refs/heads/main", timeout=300)
                break
            except RuntimeError as exc:
                if attempt == 2 or not any(m in str(exc) for m in REJECTED):
                    raise
        commit = git(root, "rev-parse", "HEAD")
        remote = git_network(root, "ls-remote", "origin", "refs/heads/main").split()[0]
        if remote != commit:
            raise RuntimeError("Remote main does not match the local checkpoint")
        return {"created_commit": bool(staged), "commit": commit,
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
