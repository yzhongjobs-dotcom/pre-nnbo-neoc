"""Exercise synchronization against a disposable local bare Git remote."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from research_sync import git, remote_matches, sync


class SyncLifecycleTest(unittest.TestCase):
    def test_real_commit_push_failure_retry_and_snapshot_guard(self):
        runtime = Path(__file__).resolve().parents[1] / ".sync-runtime"
        runtime.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="sync-test-", dir=runtime) as name:
            base = Path(name).resolve()
            # Restrict automatic temporary-directory cleanup to this workspace.
            self.assertTrue(base.is_relative_to(runtime.resolve()))
            repo, remote = base / "work", base / "remote.git"
            repo.mkdir(); remote.mkdir()
            git(remote, "init", "--bare")
            git(repo, "init", "-b", "main")
            git(repo, "config", "user.name", "Sync Test")
            git(repo, "config", "user.email", "sync-test@example.invalid")
            git(repo, "config", "core.autocrlf", "false")
            git(repo, "remote", "add", "origin", str(remote))
            frozen = repo / "00_基线快照_2026-10-14" / "reference.txt"
            frozen.parent.mkdir()
            frozen.write_bytes(b"baseline\n")
            manifest = {"files": [{"snapshot": str(frozen.relative_to(repo)),
                       "sha256": hashlib.sha256(frozen.read_bytes()).hexdigest()}]}
            (frozen.parent / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
            git(repo, "add", ".")
            git(repo, "commit", "-m", "baseline")
            git(repo, "tag", "baseline-b0-2026-10-14")
            git(repo, "push", "origin", "main")

            note = repo / "research.md"
            note.write_text("Completed research change\n", encoding="utf-8")
            done = sync(repo, "research: test completed change", str(remote))
            self.assertTrue(done["created_commit"] and done["remote_matches"])
            count = git(repo, "rev-list", "--count", "HEAD")
            repeat = sync(repo, "research: no duplicate", str(remote))
            self.assertFalse(repeat["created_commit"])
            self.assertEqual(count, git(repo, "rev-list", "--count", "HEAD"))

            frozen.write_bytes(b"changed baseline\n")
            with self.assertRaisesRegex(RuntimeError, "snapshot changed"):
                sync(repo, "must reject snapshot modification", str(remote))
            frozen.write_bytes(b"baseline\n")

            git(repo, "checkout", "-b", "other")
            with self.assertRaisesRegex(RuntimeError, "restricted to main"):
                sync(repo, "must reject other branch", str(remote))
            git(repo, "checkout", "main")

            git(repo, "config", "remote.origin.pushurl", str(base / "missing.git"))
            note.write_text("Second completed change\n", encoding="utf-8")
            old = git(repo, "rev-parse", "HEAD")
            with self.assertRaisesRegex(RuntimeError, "git push failed"):
                sync(repo, "research: keep commit on push failure", str(remote))
            retained = git(repo, "rev-parse", "HEAD")
            self.assertNotEqual(old, retained)
            self.assertEqual("", git(repo, "status", "--porcelain"))
            git(repo, "config", "--unset", "remote.origin.pushurl")
            retried = sync(repo, "research: retry pending push", str(remote))
            self.assertFalse(retried["created_commit"])
            self.assertEqual(retained, retried["commit"])
            self.assertTrue(retried["remote_matches"])


def _make_origin(base: Path):
    """Bare remote with a tagged baseline snapshot, as in the real repository."""
    remote, seed = base / "remote.git", base / "seed"
    remote.mkdir(); seed.mkdir()
    git(remote, "init", "--bare")
    git(remote, "symbolic-ref", "HEAD", "refs/heads/main")
    git(seed, "init", "-b", "main")
    _identity(seed, "Seed")
    git(seed, "remote", "add", "origin", str(remote))
    frozen = seed / "00_基线快照_2026-10-14" / "reference.txt"
    frozen.parent.mkdir()
    frozen.write_bytes(b"baseline\n")
    manifest = {"files": [{"snapshot": str(frozen.relative_to(seed)),
               "sha256": hashlib.sha256(frozen.read_bytes()).hexdigest()}]}
    (frozen.parent / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    (seed / "shared.md").write_text("line\n", encoding="utf-8")
    git(seed, "add", ".")
    git(seed, "commit", "-m", "baseline")
    git(seed, "tag", "baseline-b0-2026-10-14")
    git(seed, "push", "origin", "main", "--tags")
    return remote


def _identity(repo: Path, name: str):
    git(repo, "config", "user.name", name)
    git(repo, "config", "user.email", f"{name.lower()}@example.invalid")
    git(repo, "config", "core.autocrlf", "false")


def _clone(remote: Path, dest: Path, name: str, shallow: bool = False):
    args = ["clone", "--no-tags", "--depth", "1"] if shallow else ["clone"]
    git(dest.parent, *args, remote.as_uri() if shallow else str(remote), dest.name)
    _identity(dest, name)
    return dest


class TwoAgentTest(unittest.TestCase):
    """Codex (full local clone) and Claude Code (shallow cloud clone without
    tags) both sync to main in turn; neither may lose the other's work."""

    def test_alternating_agents_conflict_and_shallow_clone(self):
        runtime = Path(__file__).resolve().parents[1] / ".sync-runtime"
        runtime.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="sync-two-", dir=runtime) as name:
            base = Path(name).resolve()
            self.assertTrue(base.is_relative_to(runtime.resolve()))
            remote = _make_origin(base)
            codex = _clone(remote, base / "codex", "Codex")
            cc = _clone(remote, base / "cc", "Claude", shallow=True)
            self.assertEqual("", git(cc, "tag"))  # cloud clone starts without the tag
            url = remote.as_uri()

            (codex / "a.md").write_text("codex 1\n", encoding="utf-8")
            sync(codex, "research: codex change", str(remote))
            # Claude has not seen Codex's commit: sync must rebase, not reject.
            (cc / "b.md").write_text("claude 1\n", encoding="utf-8")
            done = sync(cc, "research: claude change", url)
            self.assertTrue(done["remote_matches"])
            self.assertEqual("baseline-b0-2026-10-14", git(cc, "tag"))
            log = git(cc, "log", "--format=%s", "refs/remotes/origin/main")
            self.assertIn("research: codex change", log)
            self.assertIn("research: claude change", log)

            # And back the other way.
            (codex / "c.md").write_text("codex 2\n", encoding="utf-8")
            again = sync(codex, "research: codex after claude", str(remote))
            self.assertTrue((codex / "b.md").exists() and again["remote_matches"])

            # Conflicting edits to the same line: stop, keep local commit, no force.
            (codex / "shared.md").write_text("codex version\n", encoding="utf-8")
            sync(codex, "research: codex edits shared", str(remote))
            (cc / "shared.md").write_text("claude version\n", encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "conflicting changes"):
                sync(cc, "research: claude edits shared", url)
            self.assertEqual("claude version\n", (cc / "shared.md").read_text(encoding="utf-8"))
            self.assertEqual("research: claude edits shared", git(cc, "log", "-1", "--format=%s"))
            self.assertFalse((cc / ".git" / "rebase-merge").exists())
            remote_head = git(codex, "ls-remote", str(remote), "refs/heads/main").split()[0]
            self.assertEqual(git(codex, "rev-parse", "HEAD"), remote_head)


class RemoteMatchTest(unittest.TestCase):
    def test_equivalent_spellings_only(self):
        expected = "https://github.com/yzhongjobs-dotcom/pre-nnbo-neoc.git"
        for ok in ["https://github.com/yzhongjobs-dotcom/pre-nnbo-neoc",
                   "https://github.com/YZhongJobs-dotcom/pre-nnbo-neoc/",
                   "git@github.com:yzhongjobs-dotcom/pre-nnbo-neoc.git",
                   "https://user:token@github.com/yzhongjobs-dotcom/pre-nnbo-neoc.git",
                   "http://local_proxy@127.0.0.1:43211/git/yzhongjobs-dotcom/pre-nnbo-neoc"]:
            self.assertTrue(remote_matches(ok, expected), ok)
        for bad in ["https://github.com/someone-else/pre-nnbo-neoc",
                    "https://gitlab.com/yzhongjobs-dotcom/pre-nnbo-neoc",
                    "https://github.com/yzhongjobs-dotcom/other-repo",
                    "/tmp/pre-nnbo-neoc"]:
            self.assertFalse(remote_matches(bad, expected), bad)


if __name__ == "__main__":
    unittest.main()
