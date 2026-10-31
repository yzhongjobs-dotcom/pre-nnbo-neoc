# 外部 agent 运行记录

保存 Codex CLI（含 ARIS 调研）等外部 agent 每次运行的**原始任务书与最终输出**，供以后复盘：当时问了什么、它答了什么、哪些被采纳。采纳与否及理由写在[决策记录](../decision_log.md)，不在这里写。

- 每次运行一个目录：`YYYY-MM-DD_HHMM_<agent>_<主题>/`，内含 `prompt.txt`（任务书原文）与 `output.md`（最终输出原文）。运行日志（数百 KB 的过程流）不保存。
- 输出不改写；其中的判断是该 agent 的意见，不等于项目结论，以决策记录和研究文档为准。
- 若输出另写成了项目文件（如文献调研备忘录），`output.md` 只是指针，正文在该文件：`2026-10-31_0536_codex_lit` 的正文是 [literature/aris_concurrent_survey_2026-10-31.md](../../literature/aris_concurrent_survey_2026-10-31.md)。
- 2026-10-30 至 10-31 的记录由会话临时目录补存；更早的运行未保留输出，结论见当时的 Git 提交信息。
