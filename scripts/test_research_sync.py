"""Exercise synchronization against a disposable local bare Git remote."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from research_sync import git, sync


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


if __name__ == "__main__":
    unittest.main()
