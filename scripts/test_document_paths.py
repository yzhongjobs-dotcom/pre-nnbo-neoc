"""Regression checks for prose corruption and missed terminal directory names."""
import unittest

from migrate_english_paths import rewrite_references
from project_document_policy import insert_project_note, strip_project_notes


class DocumentPathTests(unittest.TestCase):
    def test_chinese_scientific_slash_prose_is_not_a_directory(self):
        text = "直接反馈、参考/补偿；参考/前馈联合动作；原文/非零目标力适配尚待核准。"
        self.assertEqual(rewrite_references(text), text)

    def test_terminal_directory_in_execution_command(self):
        self.assertEqual(rewrite_references("cd 实验/接触约束阻抗更新\n"),
                         "cd experiments/contact_constrained_impedance_update\n")

    def test_nested_parent_and_short_child_names(self):
        self.assertEqual(rewrite_references("python 文献/资料核验.py\n"),
                         "python literature/verify_materials.py\n")
        self.assertEqual(rewrite_references("[PDF](../../文献/原文/NNBO.pdf)"),
                         "[PDF](../../literature/papers/NNBO.pdf)")

    def test_pinned_url_and_snapshot_keep_their_names(self):
        text = "https://github.com/example/repo/blob/abc/研究/方法设计.md\n`00_基线快照_2026-10-14/研究/方法设计.md`"
        self.assertEqual(rewrite_references(text), text)

    def test_real_local_link_and_encoded_local_link(self):
        self.assertEqual(rewrite_references("[方法](../研究/方法设计.md)"),
                         "[方法](../research/method_design.md)")
        self.assertEqual(rewrite_references("[闭环](%E6%8E%A5%E8%A7%A6%E9%97%AD%E7%8E%AF/README.md)"),
                         "[闭环](contact_closed_loop/README.md)")

    def test_removed_branch_provenance_remains_historical(self):
        text = "原位置 `实验/闭环策略滚动阻抗/`。"
        self.assertEqual(rewrite_references(text), text)

    def test_context_notes_do_not_change_original_body(self):
        text = "# 原标题\r\n\r\n原有中文要求和公式。\r\n"
        self.assertEqual(strip_project_notes(insert_project_note(text, "当前说明。")), text)


if __name__ == "__main__":
    unittest.main()
