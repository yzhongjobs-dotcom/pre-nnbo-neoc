# 人工原文获取：新主线的决定性近邻（一批四篇）

更新：2026-10-30。TASE2018（Santos、Cortesão）已收到并核读，进入[研究现状§6](../research/field_state.md)和[初稿第一节](../research/paper_draft_zh.md)。Wijayarathne与Hammond ICRA2020已读arXiv版（arXiv:2004.09729）。

新主线是“恒力调节下组织运动使接触增益不可观；用对运动不变、带误差界的估计重标定阻抗并做漂移前馈”。下面四篇可能已经部分做到这件事，会直接影响第一节的缺口和第四节的新颖性判断，需要原文：

| 建议文件名 | 文献 | 为什么必须看原文 |
|---|---|---|
| N01.pdf | **Implicit force control for an industrial robot based on stiffness estimation and compensation during motion**；R. Rossi、L. Fossali、A. Novazzi、L. Bascetta、P. Rocco；ICRA 2016，pp.1138–1145 | 位置接口工业机器人，在运动中在线估计刚度并补偿，正是本文“仅位置接口”条件下的最直接近邻；需核其估计是否对表面运动/几何变化不变 |
| N02.pdf | **Real-time adaptive control for haptic telemanipulation with Kalman active observers**；R. Cortesão、J. Park、O. Khatib；IEEE T-RO，2006，22(5) | TASE2018的刚度重标定与在线刚度估计来自此文（其参考文献[8]）；需核其刚度估计方式与扰动补偿是否已处理运动 |
| N03.pdf | **Model predictive control architectures with force feedback for robotic-assisted beating heart surgery**；M. Dominici、R. Cortesão；ICRA 2014，pp.2276–2282 | 主动观察器生成前馈补偿心脏运动（Wijayarathne文中[27]）；需核是否同时处理刚度变化与增益重标定 |
| N04.pdf | **Wijayarathne与Hammond ICRA2020 正式版**：Identification of Compliant Contact Parameters and Admittance Force Modulation on a Non-stationary Compliant Surface | 已读arXiv版；请下载IEEE正式版以便核对页码并归档 |

可选（若方便）：Yuan、Lee，AIM 2024，DOI 10.1109/AIM55361.2024.10637095（机器人超声中利用呼吸准周期性的自适应力控，任务近邻）。

请放到 **`E:\github_project\scan-impedance-update\literature\manual_inbox\`**，原文件名也可以。
