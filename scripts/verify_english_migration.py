"""Audit the 2026-10-20 migration without rerunning scientific experiments.

This checks the migration checkpoint, not arbitrary future scientific edits.
Original hashes describe historical evidence and are not silently regenerated.
"""
from pathlib import Path
import ast
import hashlib
import json
import re
import subprocess
from urllib.parse import unquote

from migrate_english_paths import COMPARISON, FROZEN, NAMES, ROOT, mapped_path, rewrite_references
from research_sync import verify_snapshot

TRANSLATED = {
    "experiments/experiment_index.md", "literature/literature_index.md",
    "README.md", "project_records/HANDOFF.md",
    "project_records/long_term_requirements.md",
    "project_records/decisions_and_changes.md",
    "project_records/research_failures_and_corrections.md",
    "project_records/references/method_design_workflow_and_standards.txt",
    "project_records/cleanup_record.json",
    "literature/local_materials_manifest.json",
}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def originals(paths):
    """Stream historical blobs with one Git process rather than hundreds."""
    process = subprocess.Popen(["git", "cat-file", "--batch"], cwd=ROOT,
                               stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    try:
        for path in paths:
            process.stdin.write(f"{COMPARISON}:{path}\n".encode("utf-8"))
            process.stdin.flush()
            header = process.stdout.readline().split()
            assert len(header) == 3 and header[1] == b"blob", path
            size = int(header[2])
            data = process.stdout.read(size)
            assert len(data) == size and process.stdout.read(1) == b"\n", path
            yield path, data
    finally:
        process.stdin.close()
        process.stdout.close()
        process.wait()


def main():
    target = ROOT / "project_records/english_migration.json"
    manifest = json.loads(target.read_text(encoding="utf-8"))
    counts = dict(mapped_files=0, renamed_files=0, unchanged_evidence_files=0,
                  scientific_text_files=0, python_syntax=0, json_files=0,
                  local_links=0, python_path_rewrites=0)
    failures = []
    historical_paths = subprocess.check_output(
        ["git", "ls-tree", "-rz", "--name-only", COMPARISON], cwd=ROOT).decode("utf-8").split("\0")
    before_evidence = hashlib.sha256()
    after_evidence = hashlib.sha256()
    for old, historical in originals([p for p in historical_paths if p]):
        new = mapped_path(old)
        path = ROOT / new
        assert path.is_file(), new
        data = path.read_bytes()
        counts["mapped_files"] += 1
        counts["renamed_files"] += old != new
        if not new.startswith(FROZEN + "/"):
            assert not re.search(r"[\u3400-\u9fff]", new), new
        preserved = (new.startswith((FROZEN + "/", "literature/papers/", "literature/extracted_text/"))
                     or (new.startswith("experiments/") and path.suffix not in {".py", ".md", ".txt"}))
        if preserved:
            assert data == historical, f"Changed evidence: {new}"
            prefix = old.encode("utf-8") + b"\0"
            before_evidence.update(prefix + hashlib.sha256(historical).digest())
            after_evidence.update(prefix + hashlib.sha256(data).digest())
            counts["unchanged_evidence_files"] += 1
        elif new not in TRANSLATED and new != "AGENTS.md":
            expected = digest(historical)
            if path.suffix in {".md", ".txt", ".py"} or new.startswith("literature/"):
                expected = digest(rewrite_references(historical.decode("utf-8")).encode("utf-8"))
            if new == "literature/verify_materials.py":
                old_text = historical.decode("utf-8")
                expected = digest(rewrite_references(old_text).replace(
                    "Material path outside 文献:", "Material path outside literature:").encode("utf-8"))
            assert digest(data) == expected, f"Unexpected content change: {new}"
            if path.suffix in {".md", ".txt"}:
                counts["scientific_text_files"] += 1
        # Existing scientific code may change path strings, but not algorithms.
        if path.suffix == ".py":
            expected_text = rewrite_references(historical.decode("utf-8"))
            if new == "literature/verify_materials.py":
                expected_text = expected_text.replace("Material path outside 文献:", "Material path outside literature:")
            assert ast.dump(ast.parse(data.decode("utf-8-sig"))) == ast.dump(ast.parse(expected_text)), new
            counts["python_path_rewrites"] += 1
    paths = [p for p in ROOT.rglob("*") if p.is_file() and ".git" not in p.parts
             and "__pycache__" not in p.parts]
    for path in paths:
        if path.suffix == ".py":
            ast.parse(path.read_text(encoding="utf-8-sig"), filename=str(path))
            counts["python_syntax"] += 1
        if path.suffix == ".json":
            json.loads(path.read_text(encoding="utf-8-sig"))
            counts["json_files"] += 1
        if path.suffix != ".md" or FROZEN in path.parts:
            continue
        text = path.read_text(encoding="utf-8-sig")
        for reference in re.findall(r"\]\(([^()\s]*(?:\([^()]*\)[^()\s]*)*)\)", text):
            reference = unquote(reference.split("#")[0]).strip("<>")
            if not reference or re.match(r"[a-zA-Z][\w+.-]*:", reference):
                continue
            resolved = (path.parent / reference).resolve()
            if not resolved.exists():
                failures.append([path.relative_to(ROOT).as_posix(), reference])
            counts["local_links"] += 1
    assert not failures, json.dumps(failures, ensure_ascii=False)
    # The operational translations must not retain Chinese prose. Frozen names
    # and legacy paths in the historical cleanup record are explicit exceptions.
    for name in TRANSLATED - {"project_records/cleanup_record.json", "literature/local_materials_manifest.json"}:
        text = (ROOT / name).read_text(encoding="utf-8").replace(FROZEN, "FROZEN")
        text = re.sub(r"https?://[^\s)]+", "", text)
        assert not re.search(r"[\u3400-\u9fff]", text), f"Untranslated operational text: {name}"
    verify_snapshot(ROOT)
    counts["frozen_snapshot_verified"] = True
    counts["broken_local_links"] = failures
    counts["scientific_experiments_run"] = False
    counts["original_evidence_aggregate_sha256"] = before_evidence.hexdigest()
    counts["current_evidence_aggregate_sha256"] = after_evidence.hexdigest()
    manifest.pop("files", None)
    manifest["component_map"] = NAMES
    manifest["original_files"] = counts["mapped_files"]
    manifest["renamed_files"] = counts["renamed_files"]
    manifest["note"] = "Git at comparison_commit retains original paths and bytes. Apply component_map to path components, except inside the frozen snapshot. Aggregate hashes avoid a redundant per-file historical inventory."
    manifest["validation"] = counts
    manifest["translation_files"] = sorted(TRANSLATED)
    manifest["validation_scope"] = (
        "Static dependency paths, source AST, links, JSON, scientific text/path equivalence and evidence bytes. "
        "Historical environment/source hashes remain historical. Import smoke checks and existing sync-helper tests are recorded separately.")
    target.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(counts, ensure_ascii=False))


if __name__ == "__main__":
    main()
