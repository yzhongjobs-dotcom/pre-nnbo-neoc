"""Reproduce the approved English path migration from its recorded checkpoint.

Run --apply only on the original, clean comparison commit. Historical URLs,
frozen B0 files, paper contents, and experiment data remain original evidence.
The manifest records old names and hashes; it is not a new experiment record.
"""
from pathlib import Path, PurePosixPath
import argparse
import hashlib
import json
import re
import subprocess
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
COMPARISON = "067ca1048fc492c1a61d39e88c7f839c80c29961"
FROZEN = "00_基线快照_2026-10-14"
NAMES = {
    "实验": "experiments", "文献": "literature", "研究": "research",
    "项目记录": "project_records", "参考": "references",
    "原文": "papers", "提取文本": "extracted_text",
    "逆向分析样本": "reverse_analysis_samples", "领域文献": "field_papers",
    "领域补充": "field_supplements", "方法借鉴": "method_transfer",
    "B0_问题依据核查": "b0_problem_evidence_audit",
    "MHE引导阻抗快验": "mhe_guided_impedance_pilot",
    "NNBO原文复现": "nnbo_source_reproduction",
    "两候选分线验证": "separate_candidate_validation",
    "双候选任务仿真": "dual_candidate_task_simulation",
    "最优阻抗更新": "optimal_impedance_update",
    "动态接触强对照": "dynamic_contact_strong_comparators",
    "参数与观测边界": "parameter_observation_boundaries",
    "同次扫描阻抗反馈": "within_scan_impedance_feedback",
    "持续缓变更新": "continuous_slow_variation_updates",
    "接触历史机制核查": "contact_history_mechanism_audit",
    "接触约束阻抗更新": "contact_constrained_impedance_update",
    "有限时域接触": "finite_horizon_contact",
    "空间接触估计核查": "spatial_contact_estimation_audit",
    "约束切换与值根": "constraint_switching_value_roots",
    "观测与局部根": "observation_local_roots",
    "接触闭环": "contact_closed_loop",
    "01_复用动作计算灵敏度": "01_reused_action_sensitivity",
    "02_平滑二次critic对照": "02_smooth_quadratic_critic_comparator",
    "03_小邻域与密集更新_探索性": "03_small_neighborhood_dense_updates_exploratory",
    "研究问题.md": "research_question.md", "研究现状.md": "field_state.md",
    "方法设计.md": "method_design.md", "技术附录.txt": "technical_appendix.txt",
    "理论证明.txt": "theoretical_proofs.txt",
    "最优阻抗更新_问题与方法.md": "optimal_impedance_question_and_methods.md",
    "两候选分线方法设计与迭代.md": "separate_candidate_design_and_iteration.md",
    "两候选改进_定向文献与设计.md": "candidate_improvements_targeted_literature_and_design.md",
    "可借鉴机制.md": "transferable_mechanisms.md",
    "顶刊论文问题与方法设计_独立逆向分析.md": "independent_reverse_analysis_of_top_journal_methods.md",
    "长期要求.md": "long_term_requirements.md",
    "研究决策与变更.md": "decisions_and_changes.md",
    "科研失败与纠偏.md": "research_failures_and_corrections.md",
    "清理记录.json": "cleanup_record.json",
    "Method_Design_科研流程与标准.txt": "method_design_workflow_and_standards.txt",
    "顶会顶刊_Method_经验标尺_Reference_Distribution_v1_1_Audited.txt": "top_venue_method_reference_distribution_v1_1_audited.txt",
    "实验索引.md": "experiment_index.md", "文献索引.md": "literature_index.md",
    "候选1报告.md": "candidate1_report.md", "候选2报告.md": "candidate2_report.md",
    "候选说明.md": "candidate_description.md", "测试报告.md": "test_report.md",
    "结果.md": "results.md", "独立结果.md": "independent_results.md",
    "设计与协议.md": "design_and_protocol.md",
    "平滑对照的解析条件.md": "analytic_conditions_for_smooth_comparator.md",
    "本地资料清单.json": "local_materials_manifest.json",
    "资料核验.py": "verify_materials.py", "顶刊样本复核.py": "verify_top_venue_samples.py",
    "阻抗文献重建_来源与筛选.json": "impedance_literature_sources_and_screening.json",
    "顶刊样本数学核对.json": "top_venue_sample_math_checks.json",
    "顶刊样本来源与阅读范围.json": "top_venue_sample_sources_and_reading_scope.json",
}


def mapped_path(path):
    if path == FROZEN or path.startswith(FROZEN + "/"):
        return path
    return "/".join(NAMES.get(part, part) for part in PurePosixPath(path).parts)


def rewrite_references(text):
    """Replace path components without translating prose or historical URLs."""
    protected = []
    def protect(match):
        protected.append(match.group())
        return f"\x00PROTECTED{len(protected)-1}\x00"
    # Frozen paths and external URLs must retain their original spelling.
    text = re.sub(r"https?://[^\s<>\"')]+|" + re.escape(FROZEN) + r"[^\s`\"')<>]*", protect, text)
    # Some local Markdown targets use URL-encoded Chinese names. External URLs
    # were protected above and must continue to address their pinned versions.
    text = re.sub(r"(?<=\]\()([^\s)]+)(?=\))",
                  lambda m: mapped_path(unquote(m.group())) if "%" in m.group() else m.group(), text)
    for old, new in sorted(NAMES.items(), key=lambda pair: len(pair[0]), reverse=True):
        if "." in old:
            text = text.replace(old, new)
        else:
            # At a path boundary, or as a quoted component in Path expressions.
            text = re.sub(r"(?<![\w\u3400-\u9fff])" + re.escape(old) + r"(?=[/\\])", new, text)
            for quote in ("'", '"', '`'):
                text = text.replace(quote + old + quote, quote + new + quote)
    for i, value in enumerate(protected):
        text = text.replace(f"\x00PROTECTED{i}\x00", value)
    return text


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def migrate():
    current = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    if current != COMPARISON:
        raise RuntimeError("Migration is restricted to the recorded comparison commit")
    if subprocess.run(["git", "diff", "--quiet", "HEAD", "--"], cwd=ROOT).returncode:
        raise RuntimeError("Commit or preserve existing tracked changes before migration")
    tracked = subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT).decode("utf-8").split("\0")
    tracked = [path for path in tracked if path]
    destinations = [mapped_path(path) for path in tracked]
    if len(set(path.casefold() for path in destinations)) != len(destinations):
        raise RuntimeError("Case-insensitive destination collision")
    records = []
    for old, new in zip(tracked, destinations):
        source, target = ROOT / old, ROOT / new
        if not source.is_file():
            raise RuntimeError(f"Missing migration source: {old}")
        before = source.read_bytes()
        after = before
        # Historical JSON/CSV/arrays/models are evidence, not live path indexes.
        live_json = old.startswith("文献/") and old.endswith(".json")
        editable = source.suffix in {".py", ".md", ".txt"} or live_json
        protected = old.startswith((FROZEN + "/", "文献/原文/", "文献/提取文本/"))
        if editable and not protected:
            after = rewrite_references(before.decode("utf-8")).encode("utf-8")
        if old != new:
            if target.exists():
                raise RuntimeError(f"Destination exists: {new}")
            if not target.resolve().is_relative_to(ROOT.resolve()):
                raise RuntimeError("Migration destination escaped the workspace")
            target.parent.mkdir(parents=True, exist_ok=True)
            source.rename(target)
        if after != before:
            target.write_bytes(after)
        records.append(dict(old_path=old, path=new, before_sha256=sha256(before),
                            after_path_rewrite_sha256=sha256(after)))
    # Remove only empty original directories; never delete ignored/untracked files.
    for directory in sorted(ROOT.rglob("*"), key=lambda path: len(path.parts), reverse=True):
        if directory.is_dir() and ".git" not in directory.parts and not directory.name == "__pycache__":
            if not any(directory.iterdir()):
                directory.rmdir()
    manifest = dict(comparison_commit=COMPARISON, date="2026-10-20", component_map=NAMES,
                    original_files=len(records), renamed_files=sum(a != b for a,b in zip(tracked,destinations)),
                    scope="English paths and operational documents; scientific prose and experiment explanations remain Chinese.",
                    historical_evidence="Original experiment data, paper contents, frozen snapshot and historical URLs are preserved. Old provenance paths resolve through this map.",
                    note="Git at comparison_commit retains original paths and bytes. Apply component_map to path components, except inside the frozen snapshot. Validation records aggregate evidence hashes instead of duplicating a per-file historical inventory.")
    (ROOT / "project_records/english_migration.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(dict(files=len(records), renamed=sum(a != b for a,b in zip(tracked,destinations)))))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    if args.apply:
        migrate()
    else:
        print("Use --apply at the original comparison commit; see project_records/english_migration.json for the completed migration.")
