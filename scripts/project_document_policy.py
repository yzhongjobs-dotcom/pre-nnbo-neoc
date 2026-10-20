"""Current document-language policy and explicitly recorded reference repairs.

The comparison commit preserves original Chinese wording. Only working paths,
the corrections below, and delimited current-context notes may differ in the
restored documents. This is an integrity check, not a scientific proof audit.
"""
import re

CHINESE_RESTORED = {
    "README.md",
    "project_records/HANDOFF.md",
    "project_records/long_term_requirements.md",
    "project_records/decisions_and_changes.md",
    "project_records/research_failures_and_corrections.md",
    "project_records/references/method_design_workflow_and_standards.txt",
    "experiments/experiment_index.md",
    "literature/literature_index.md",
}

NOTE_START = "<!-- project-update:start -->"
NOTE_END = "<!-- project-update:end -->"
CLEANUP_SCOPE_NOTE = (
    "本记录描述2026-10-19清理时的事实、路径、数量与保护范围，不代表当前目录仍使用旧名称。"
    "旧路径和原始版本保留用于历史回查；当前位置见english_migration.json。后续用户确认与当前任务以HANDOFF和长期要求为准。"
)

CORE_REFERENCE_ALIASES = {
    "R15_RESEARCH_PROBLEM_SPATIAL_PREVIEW_CONTACT_MODEL_APPENDIX_2026-09-12.txt": "technical_appendix.txt",
    "R15_RESEARCH_PROBLEM_SPATIAL_PREVIEW_CONTACT_MODEL_2026-09-12.md": "method_design.md",
    "R15_ANEOC_CERTIFIED_DEPLOYMENT_PROOFS_2026-09-13.txt": "theoretical_proofs.txt",
}
CORE_FILES = {"research/method_design.md", "research/technical_appendix.txt", "research/theoretical_proofs.txt"}

TEXT_CORRECTIONS = {
    "research/field_state.md": {
        "进一步的标准分析留到新对话。":
        "后续经用户确认的方法形成认识与阶段要求已落实长期要求；最新候选进度见HANDOFF，本段不构成待审核或暂停指令。",
    },
    "project_records/references/top_venue_method_reference_distribution_v1_1_audited.txt": {
        "当前定位：": "历史定位（已由后续确认的长期要求调整，不作为当前执行指令）：",
        "可以作为 R15 正式 Method 审计的主经验标尺，但最终 novelty verdict 仍必须依赖 R15 自身 strongest-prior full-text research，而不能只凭本文件。":
        "当时曾定位为 R15 Method 审计的主经验标尺；当前只保留为历史样本与经验参考。科研判断依照用户最新确认的长期要求、具体原文和实际证据，不将本文固定篇数、先选全部对照或统一门槛恢复为强制流程。",
    },
}


def documented_reference_repairs(path, text):
    if path in CORE_FILES:
        for old, new in CORE_REFERENCE_ALIASES.items():
            text = text.replace(old, new)
    for old, new in TEXT_CORRECTIONS.get(path, {}).items():
        text = text.replace(old, new)
    return text


def strip_project_notes(text):
    pattern = r"(?ms)^" + re.escape(NOTE_START) + r"\r?\n.*?^" + re.escape(NOTE_END) + r"\r?\n"
    return re.sub(pattern, "", text)


def insert_project_note(text, note):
    assert NOTE_START not in text and NOTE_END not in text
    split = text.find("\n") + 1
    assert split > 0
    block = NOTE_START + "\n\n" + note.strip() + "\n\n" + NOTE_END + "\n"
    return text[:split] + block + text[split:]
