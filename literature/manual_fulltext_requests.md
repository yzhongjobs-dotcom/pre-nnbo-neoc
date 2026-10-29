# 本轮人工原文状态与剩余两项

更新：2026-10-29。本轮已核查你提供的11个PDF命名文件：**10份有效原文已保存并完成所列选段核读；1份是安全验证HTML。** 原8篇ARIS下载批次已完成，不需要再次下载。当前按你的要求交付现状补充并停下等待验收；以下是仍影响完整性判断的两项原文问题。

## 尚未闭合的两项

| 编号／建议文件名 | 文献 | 需要解决什么 |
|---|---|---|
| S03.pdf | **Autonomous robotic ultrasound scanning on complex surfaces: Automatic trajectory generation and disturbance observer-based model predictive admittance control**；[10.1016/j.conengprac.2026.107240](https://doi.org/10.1016/j.conengprac.2026.107240) | 现有同名文件62900字节，以`<!DOCTYPE`开头，内容是Elsevier Security verification，未取得论文正文。须用有效PDF／完整作者稿替换；既有索引为2027年1月期次，最早在线日期也待核准。该文可能改变观察器与预测导纳结合的任务近邻判断。 |
| C09.pdf | **ACT-US: End-to-End Visuomotor Policy Learning for Respiratory-Adaptive Robotic Ultrasound**；IJCNN，元数据2026；[10.1109/IJCNN69854.2026.11728409](https://doi.org/10.1109/IJCNN69854.2026.11728409) | 独立覆盖审查指出的呼吸适应超声策略近邻。身份核准，多次元数据补核仍无实际摘要、本地原文或已声明arXiv身份。需要核清动作、视觉补偿、学习和导纳的接口及最早公开时间，不能依据检索片段认定能力。 |

仍可放入本地 `E:\github_project\scan-impedance-update\literature\manual_inbox\`。只需这两项有更新时告知，无需手工填写JSON。正式Controller当前只请求C09；S03是既有关键缺口，二者在本表集中展示，避免逐篇往返。其他文件保持原状，不删除或覆盖用户提供的原始材料。

## 已处理的10份有效原文

| 编号 | 文献及处理位置 |
|---|---|
| C01 | Bi等机器人超声学习综述：交叉校准范围和术语，不堆入主比较表 |
| C02 | Roveda/Piga无力传感环境估计：在线辨识到阻抗规则 |
| C03 | Wu轻量超声恒力系统：主动／被动硬件任务路线 |
| C04 | Huang自适应鲁棒力控制：直接调力机制 |
| C05 | Zeng能量感知力导纳：支持性执行基础，简述装配条件 |
| C06 | Duan约束学习超声策略：学习目标力增量，与增益学习区别 |
| C07 | Chalasani连续触诊：联合几何／刚度估计与采样条件 |
| C08 | Zhu软组织切割CBF：执行机制资源，公开日期未核准前不列为截止日前核心能力 |
| S01 | DOSMAC：已由摘要线索深化为原文控制式、条件及实验 |
| S02 | MHE-VIC：已由摘要线索深化为原文观察器、增益规则及加工收益 |

此外复用了库内E27、G23、E41原文，补齐学习阻抗和预测增益优化的正式证据。每份实际核读范围、文件哈希和用途见[统一清单](local_materials_manifest.json)；研究现状正文见[field_state](../research/field_state.md)。

[ARIS research-lit/SKILL.md](<D:/桌面/科研Agent Harness设计/ARIS/skills/skills-codex/research-lit/SKILL.md>)明确规定："Automatic full-text acquisition is arXiv-only." 其余已选原文按一次人工批次处理。本清单说明证据缺口，不是再次请求检索授权。未取得原文时不宣布该路线不存在，也不宣布研究现状已查全。
