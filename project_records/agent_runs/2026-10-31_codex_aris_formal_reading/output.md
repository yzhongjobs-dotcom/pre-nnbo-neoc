# ARIS formal run report 2

run_id: `recent-related-work-20261031`  
RUN_ROOT: `D:/pip_env/temp/scan_impedance_literature_20261031/aris_run`  
policy sha256: `b467b81f955d1690249d3e803f7c4fca75db9115f80f704d254a53e91cdc8665`  
query-plan sha256: `1cac16a4e916603dd43992ed118794bab96b4d9ad9612b0bbd0032a3746c36c4`

## A. Controller 最终状态

Controller 当前为 `HUMAN_SEARCH_REQUIRED`，`waiting_for=human_fulltext_batch`。本次唯一的 `finish-reading` 已完成；停止原因是：42 篇 `ADMIT_FOR_READING` 中，5 篇有稳定 arXiv URL，已由 Controller 读取并提交 Evidence；其余 37 篇无 arXiv 全文，按政策统一 defer 到一个人工全文批次。Controller 请求对象中的 `failed_paper_id` 记录为 `jkUy-uxymaMJ`，但本次不是该论文的 provider 失败，而是非 arXiv 批量人工下载关口。

## B. 实际检索与覆盖

- provider：SerpApi Google Scholar（`serpapi_google_scholar`，Google Scholar engine）。31 个 schema-v2 查询全部完成，包含政策全部查询轴、六条前向引文链、定向细化和 venue 分层；无 provider 故障、无 `HUMAN_SEARCH_REQUIRED` 的人工检索关口。
- 返回/去重：每个查询按计划返回 20 个结果位，共 620 个 SerpApi 结果位；Controller 去重后为 347 个唯一候选。由于候选后续 metadata/enrichment 会在 corpus JSONL 中追加记录，报告以 Controller 当前 347 个 paper ID 为去重口径。
- admission：`ADMIT_FOR_READING` 42 篇，`ADMIT_DISCOVERY_ONLY` 5 篇；另有 `HOLD_IDENTITY` 6 篇、`EXCLUDE_IRRELEVANT` 294 篇。in-scope 候选按真实题名+摘要筛查；明确不相关项按规则使用 `TITLE_ONLY` 快速排除，身份未核实项未进入阅读。
- coverage reviewer：尚未产生独立 native coverage reviewer 结论；Controller 在全文人工关口前尚未创建 Field Map/coverage review request。因此不能把本轮检索完成误报为覆盖闭合。未完成覆盖是：37 篇非 arXiv 的人工全文、其后逐篇 Evidence、Field Map 以及独立 coverage reviewer 审查。

表中缩写：信息来源 `F=力/力矩/接触力，I=超声图像/置信度，V=外部视觉/触觉，K=运动学/几何/动力学模型，E=专家示范/学习策略，S=生理或任务状态`；“已知”只表示政策中的既有记录/锚点，不等于本轮已形成正式 Evidence。

## C. admitted 文献

|编号|完整题名|作者|年份|venue|DOI/稳定身份|主题类|刚度/运动信息来源|执行接口|研究现状中的作用|全文状态|已知文献|
|---|---|---|---:|---|---|---|---|---|---|---|---|
|1|Technical challenges in autonomous robotic ultrasound examinations: perception, planning, and control|Dong Guo; Yongde Zhang; Xuequan Huang; Chuang He|2026|Journal of Robotic Surgery|10.1007/s11701-026-03509-z|综述|I/F/K|综述|背景/发现层|发现层不读|否|
|2|A Survey of Autonomous Robotic Ultrasound Scanning Systems|Khushboo Munir; Abdullah F. Al-Battal; Ammar Alsheghri; Harald Becher; Michelle Noga; Kumaradevan Punithakumar|2025|IEEE Access|10.1109/access.2025.3574464|综述|I/F/K|综述|背景|待人工|否|
|3|Toward Clinical Applications of Intelligent Robotic Ultrasound Systems|Taiyu Han; Guochen Ning; Hanying Liang; Zihan Li; Zhongliang Jiang; Fang Chen; Yan Kang; Jianwen Luo; Hongen Liao|2026|IEEE Reviews in Biomedical Engineering|10.1109/rbme.2025.3610605|综述|I/F/K|综述|背景|待人工|否|
|4|Robot-Assisted Medical Imaging: A Review|Septimiu E. Salcudean; Hamid Moradi; David G. Black; Nassir Navab|2022|Proceedings of the IEEE|10.1109/jproc.2022.3162840|综述|I/F/K|综述|背景|待人工|否|
|5|A Review and Perspective of Techniques for Autonomous Robotic Ultrasound Acquisitions|Yanding Qin; Lele Dang; Fan Ren; Zhuomao Li; Lijun Duan; Hongpeng Wang; Jianda Han|2026|Sensors|10.3390/s26072081|综述|I/F/K|综述|背景|待人工|否|
|6|A shared-control ultrasound-guided robotic system for liver puncture under respiratory motion|Wenlong Zhang; Weihao Lu; Qijie Lu; Wentuo Yang; Le Xie; Yuchen Yang|2026|Scientific Reports|10.1038/s41598-026-56017-3|呼吸/任务|F/I/S|位姿-力共享控制|最近邻|待人工|否|
|7|A review of robot-assisted ultrasound examination: Systems and technology|Haiyan Du; Xinran Zhang; Yongde Zhang; Fujun Zhang; Letao Lin; Tao Huang|2024|The International Journal of Medical Robotics and Computer Assisted Surgery|10.1002/rcs.2660|综述|I/F/K|综述|背景|待人工|否|
|8|Task Space Compliant Control and Six-Dimensional Force Regulation Toward Automated Robotic Ultrasound Imaging|Junchen Wang; Chunheng Lu; Yifei Lv; Siqin Yang; Mingbo Zhang; Yu Shen|2024|IEEE Transactions on Automation Science and Engineering|10.1109/tase.2023.3282974|直接超声控制|F/K|六维导纳/力控|路线代表|待人工|否|
|9|A Lightweight Robotic Ultrasound System With Hybrid Force Control and Modular Force-Range Switching|Zehao Wu; Weijian Zhang; Xianli Wang; I-Ming Chen; Qingsong Xu|2026|IEEE Transactions on Robotics|10.1109/tro.2026.3734315|直接超声控制|F/K|混合力控|最近邻|待人工|否|
|10|Design of a Cost-Effective Ultrasound Force Sensor and Force Control System for Robotic Extra-Body Ultrasound Imaging|Yixuan Zheng; Hongyuan Ning; Eason Rangarajan; Aban Merali; Adam Geale; Lukas Lindenroth; Zhouyang Xu; Weizhao Wang; Philipp Kruse; Steven Morris; Liang Ye; Xinyi Fu; Kawal Rhode; Richard James Housden|2025|Sensors|10.3390/s25020468|直接超声控制|F/K|力传感器+混合位置力控|路线代表|待人工|否|
|11|Autonomous Robotic Ultrasound Scanning for Liver Interventions: an Integrated Redundancy Resolution and Motion Control Framework|Peng Chen; Hiroshi Yokoi; Bo Zhang|2026|Intelligent Systems & Robotic Mechanics|10.62836/isrm.v3i1.0001|任务/运动|I/F/K|冗余解析+运动控制|最近邻|待人工|否|
|12|ACT-US: End-to-End Visuomotor Policy Learning for Respiratory-Adaptive Robotic Ultrasound|Xingyu Zhang; Ke Lu; Hongjuan Pei; Qingyuan Liu; Yihua Shao; Xiang Luo; Jian Xue|2026|2026 International Joint Conference on Neural Networks (IJCNN)|10.1109/ijcnn69854.2026.11728409|学习/呼吸|I/S/E|端到端视觉运动策略|发现层|发现层不读|否|
|13|Learning Variable Impedance Control for Contact Sensitive Tasks|Miroslav Bogdanovic; Majid Khadiv; Ludovic Righetti|2020|IEEE Robotics and Automation Letters|10.1109/lra.2020.3011379|学习阻抗|F/K/E|学习变阻抗|路线代表|待人工|否|
|14|On the Explicit Robust Force Control via Disturbance Observer|Emre Sariyildiz; Kouhei Ohnishi|2015|IEEE Transactions on Industrial Electronics|10.1109/tie.2014.2361611|鲁棒力控|F/K|扰动观测器力控|背景机制|待人工|否|
|15|Minimally Model-Based Trajectory Tracking and Variable Impedance Control of Flexible-Joint Robots|Emmanouil Spyrakos-Papastavridis; Jian S. Dai|2021|IEEE Transactions on Industrial Electronics|10.1109/tie.2020.2994886|变阻抗|K/F|最小模型轨迹+变阻抗|背景机制|待人工|否|
|16|Tactile-Guided Robotic Ultrasound: Mapping Preplanned Scan Paths for Intercostal Imaging|Yifan Zhang; Dianye Huang; Nassir Navab; Zhongliang Jiang|2025|2025 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS)|10.1109/iros60139.2025.11246299|触觉/运动补偿|F/I/K|触觉引导路径映射|最近邻|待人工|是（R08）|
|17|Dual-Rate Force-Image Control with Model-Based Orientation Limits for Robotic Ultrasound|Tyler Foster; Qiang Zhang; A B M Tahidul Haque; Anh Thu Nguyen|2026|arXiv preprint|arXiv:2610.05839v1|模型/多速率|F/I/K|双速率力控+姿态速率上限|路线代表|arXiv已读|否|
|18|AuSoScan: Automatic Scoliosis Assessment by Ultrasound Scanning With Soft Contact Control|Anqing Duan; Wanli Liuchen; Peng Zhou; Dezhen Song; Chenguang Yang; Yongping Zheng; David Navarro-Alarcon|2025|IEEE/ASME Transactions on Mechatronics|10.1109/tmech.2025.3583041|任务/软接触|I/F/K|软接触扫描控制|最近邻|待人工|否|
|19|Ultrasound-Guided Assistive Robots for Scoliosis Assessment With Optimization-Based Control and Variable Impedance|Anqing Duan; Maria Victorova; Jingyuan Zhao; Yuxiang Sun; Yongping Zheng; David Navarro-Alarcon|2022|IEEE Robotics and Automation Letters|10.1109/lra.2022.3186504|模型/变阻抗|I/F/K|优化控制+变阻抗|路线代表|待人工|否|
|20|Inverse-Reinforcement-Learning-Based Robotic Ultrasound Active Compliance Control in Uncertain Environments|Guochen Ning; Hanying Liang; Xinran Zhang; Hongen Liao|2024|IEEE Transactions on Industrial Electronics|10.1109/tie.2023.3250767|学习合规|F/I/E|逆强化学习主动合规|路线代表|待人工|否|
|21|Compliance Control of a Robotic Breast Ultrasound System Based on Variable Admittance Control|Yongde Zhang; Dong Guo; Xuequan Huang; Chuang He|2026|The International Journal of Medical Robotics and Computer Assisted Surgery|10.1002/rcs.70194|直接超声控制|F/I|变导纳|最近邻|待人工|否|
|22|Manifold Generative Annealing for Contact Constraint-aware Impedance Control with Reinforcement Learning Prior|Zhilin He; Guoquan Qi; Changliu Liu|2026|arXiv preprint|arXiv:2610.21140v1|模型-学习混合|F/K/E|短时域运动-阻抗优化|路线代表|arXiv已读|否|
|23|Explicit force-based control strategies for compliant robot-assisted ultrasound imaging|Adrian Piedra; R Brooke Jeffrey; Oussama Khatib|2026|The International Journal of Robotics Research|10.1177/02783649261441635|直接超声控制|F/K|显式力控|路线代表|待人工|否|
|24|ForceU-VLA: A Force-Aware Vision-Language-Action Model for Embodied Ultrasound Scanning|Xingzheng Wu; Cheng Zhang; Guihao Yan; Xifeng Hu; Zhi Liu; Qing Cai|2026|arXiv preprint|arXiv:2608.15009v1|学习阻抗/动作|F/I/V/E|多模态 VLA 动作|路线代表|arXiv已读|否|
|25|Dual-Agent-Based Robotic Ultrasound Path Planning and Interaction Control in External-Vision-Independent Environments|Xinye Wang; Zhiyuan He; Peng Chen; Zhe Wang; Tao Sun|2025|IEEE Transactions on Automation Science and Engineering|10.1109/tase.2025.3614996|规划/交互|F/I/K|双智能体路径+交互控制|最近邻|待人工|否|
|26|Breaking the Latency Barrier: Synergistic Perception and Control for High-Frequency 3D Ultrasound Servoing|Yizhao Qian; Yujie Zhu; Jiayuan Luo; Li Liu; Yixuan Yuan; Hongen Liao; Guochen Ning|2026|2026 IEEE International Conference on Robotics and Automation (ICRA)|10.1109/icra57385.2026.11697234|延迟/运动|I/K|高频超声伺服|最近邻|待人工|否|
|27|Contact-Aware Impedance Controller for Robot-Assisted Ultrasound Imaging|MD Miraj Arefin; M Efe Tiryaki|2026|arXiv preprint|arXiv:2609.34400v1|接触/阻抗|I/F/K|图像接触感知阻抗|最近邻|arXiv已读|否|
|28|Impedance Control on Arbitrary Surfaces for Ultrasound Scanning Using Discrete Differential Geometry|Michael Dyck; Arne Sachtler; Julian Klodmann; Alin Albu-Schaffer|2022|IEEE Robotics and Automation Letters|10.1109/lra.2022.3184800|几何/阻抗|I/K|离散微分几何阻抗|路线代表|待人工|否|
|29|Optimization-Based Variable Impedance Control of Robotic Manipulator for Medical Contact Tasks|Junling Fu; Ilaria Burzo; Elisa Iovene; Jianzhuang Zhao; Giancarlo Ferrigno; Elena De Momi|2024|IEEE Transactions on Instrumentation and Measurement|10.1109/tim.2024.3372209|模型/优化|F/K|优化变阻抗|模型路线代表|待人工|否|
|30|Precise Repositioning of Robotic Ultrasound: Improving Registration-Based Motion Compensation Using Ultrasound Confidence Optimization|Zhongliang Jiang; Nehil Danis; Yuan Bi; Mingchuan Zhou; Markus Kroenke; Thomas Wendler; Nassir Navab|2022|IEEE Transactions on Instrumentation and Measurement|10.1109/tim.2022.3200360|图像/运动补偿|I/K|置信度优化重定位|最近邻|待人工|是（R01）|
|31|Trajectory-aware feedforward admittance control for contact force tracking with time-varying reference trajectory in a teleoperated robotic ultrasound system|Sen Li; Yufan Pi; Weihua Li; Giuseppe Carbone; Kai Wu|2026|Control Engineering Practice|10.1016/j.conengprac.2025.106674|导纳/轨迹|F/K|轨迹前馈导纳|发现层|发现层不读|否|
|32|Dynamic Arm Impedance Regulation: Experiences From Medical Ultrasound Scanning Experts|Baoshan Niu; Dapeng Yang; Tao Geng; Le Zhang; Boyu Ma; Qi Huang; Hong Liu|2024|IEEE Transactions on Instrumentation and Measurement|10.1109/tim.2024.3428622|人类策略/阻抗|E/F/I|动态臂阻抗调节|路线代表|待人工|否|
|33|Observer-based sliding mode admittance control for ultrasound robot force-tracking in complex interaction environments|Sen Li; Qi Chen; Weihua Li; Mingfeng Wang; Sungsung Pang; Kai Wu|2026|Control Engineering Practice|10.1016/j.conengprac.2026.106764|鲁棒导纳|F/K|滑模+观测器导纳|发现层|发现层不读|否|
|34|Safe Contact Force Generation for Robotic Thyroid Ultrasound Imaging|Ryosuke Tsumura; Toshifumi Tomioka; Yoshihiko Koseki; Kiyoshi Yoshinaka|2024|IEEE Robotics and Automation Letters|10.1109/lra.2024.3349959|安全/力控|F/I|安全接触力生成|最近邻|待人工|是（R09）|
|35|Smooth Path Planning and Dynamic Contact Force Regulation for Robotic Ultrasound Scanning|Le Zhang; Dapeng Yang; Baoshan Niu; Haonan Yang; Qi Huang; Li Jiang; Hong Liu|2025|IEEE Robotics and Automation Letters|10.1109/lra.2025.3604746|规划/动态力控|F/I/K|平滑路径+动态力调节|最近邻|待人工|否|
|36|EMG-based Hybrid Impedance-Force Control for Human-Robot Collaboration on Ultrasound Imaging|Teng Li; Hongjun Xing; Hamid D. Taghirad; Mahdi Tavakoli|2022|2022 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS)|10.1109/iros47612.2022.9981615|人机协同|F/E|EMG 调制混合阻抗-力控|路线代表|待人工|否|
|37|A Novel Ultrasound Robot With Force/Torque Measurement and Control for Safe and Efficient Scanning|Xianqiang Bao; Shuangyi Wang; Lingling Zheng; Richard James Housden; Joseph V. Hajnal; Kawal Rhode|2023|IEEE Transactions on Instrumentation and Measurement|10.1109/tim.2023.3239925|直接超声控制|F/K|力矩测量+恒力/深度接口|路线代表|待人工|否|
|38|Optimizing Variable Admittance Control for Remote Ultrasound Scanning Under Uncertain Environment|Songjie Xiao; Tengyue Wang; Yun Long; Liangjing Yang|2025|IEEE Access|10.1109/access.2025.3568150|在线变导纳|F/I/K|环境估计+变导纳|最近邻|待人工|否|
|39|A Unified Interaction Control Framework for Safe Robotic Ultrasound Scanning with Human-Intention-Aware Compliance|Xiangjie Yan; Shaqi Luo; Yongpeng Jiang; Mingrui Yu; Chen Chen; Senqiang Zhu; Gao Huang; Shiji Song; Xiang Li|2024|2024 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS)|10.1109/iros58592.2024.10801755|安全/人机交互|F/V/K|统一交互合规|路线代表|待人工|否|
|40|Passivity Preservation for Variable Impedance Control of Compliant Robots|Emmanouil Spyrakos-Papastavridis; Peter R. N. Childs; Jian S. Dai|2020|IEEE/ASME Transactions on Mechatronics|10.1109/tmech.2019.2961478|被动性/阻抗|K/F|被动性保持变阻抗|背景机制|待人工|否|
|41|Simultaneous Online Registration-Independent Stiffness Identification and Tip Localization of Surgical Instruments in Robot-Assisted Eye Surgery|Ali Ebrahimi; Shahriar Sefati; Peter Gehlbach; Russell H. Taylor; Iulian I. Iordachita|2023|IEEE Transactions on Robotics|10.1109/tro.2022.3201393|跨域刚度估计|F/V/K|在线刚度识别+状态估计|机制迁移背景|待人工|否|
|42|Diffusion-Based Impedance Learning for Contact-Rich Manipulation Tasks|Noah Geiger; Tamim Asfour; Neville Hogan; Johannes Lachner|2026|IEEE Transactions on Robotics|10.1109/tro.2026.3717255|学习阻抗|F/E|扩散+能量一致阻抗|路线代表|待人工|否|
|43|Quasi-Static Elastography-Driven Automated Robotic Ultrasound Screening and Localization|Hanying Liang; Shipeng Zhang; Guochen Ning; Jianwen Luo; Hongen Liao|2026|IEEE Transactions on Biomedical Engineering|10.1109/tbme.2025.3634552|弹性成像/任务|I/F/K|顺应力控+闭环定位|最近邻|待人工|否|
|44|Acoustic Shadowing Aware Robotic Ultrasound: Lighting up the Dark|Viviana Sutedjo; Maria Tirindelli; Christine Eilers; Walter Simson; Benjamin Busam; Nassir Navab|2022|IEEE Robotics and Automation Letters|10.1109/lra.2022.3141451|图像/路径优化|I/K|超声姿态优化|背景/任务收益|待人工|否|
|45|A Passive Variable Impedance Control Strategy with Viscoelastic Parameters Estimation of Soft Tissues for Safe Ultrasonography|Luca Beber; Edoardo Lamon; Davide Nardi; Daniele Fontanelli; Matteo Saveriano; Luigi Palopoli|2024|2024 IEEE International Conference on Robotics and Automation (ICRA)|10.1109/icra57147.2024.10610167|软组织/变阻抗|F/I/K|粘弹参数估计+QP变阻抗|最近邻/路线代表|待人工|是（Beber ICRA2024）|
|46|DAISS: Phase-Aware Imitation Learning for Dual-Arm Robotic Ultrasound-Guided Interventions|Feng Li; Pei Liu; Shiting Wang; Ning Wang; Zhongliang Jiang; Nassir Navab; Yuan Bi|2026|arXiv preprint|arXiv:2603.07663v1|学习运动/介入|I/V/E|阶段感知双臂动作|路线代表|arXiv已读|否|
|47|Toward Automating Transrectal Prostate Ultrasound Procedure: A Robotic Platform With Safe Probe Control Methods|Seungmin Lee; Joanne Yoon; June Young Seo; Sung-Joon Ye; Young-bong Bang|2026|IEEE Robotics and Automation Letters|10.1109/lra.2026.3746113|任务/安全|F/K|安全探头控制|发现层|发现层不读|否|

## D. 已读 arXiv 全文的机制要点与影响

1. **Contact-Aware Impedance Controller（edQRCdQdV8AJ）**：做到从 B-mode 图像检测声学接触、估计曲面探头有效接触点，并更新 contact frame 与 wrench application point；依赖图像能稳定反映接触、已知探头几何和琼脂仿真体条件。它把“接触状态/接触几何”补进控制链，但尚未在线优化阻抗参数或直接证明任务级收益。
2. **Dual-Rate Force-Image Control（1qrACtWl090J）**：把 1 kHz 力回路与 30 Hz、有延迟图像回路的耦合显式化，用局部刚度、力回路增益、旋转到力增益上界、预算和预测时域给出可执行姿态速率上限；依赖关节力矩力估计、接触刚度辨识和局部模型。主要推进是“模型/约束设计在线运动接口”，而非在线预测阻抗。
3. **ForceU-VLA（cbZN0njoA9sJ）**：用 FUSFM 联合建模超声与 6D 力/力矩，用 SAMM 按扫描阶段改变模态权重，再输出连续动作；依赖同步的视觉、超声、力和动作示范数据，以及训练分布对接触阶段的覆盖。它代表学习路线把实际力反馈并入决策，但阻抗参数和物理接触模型仍不是显式在线优化变量。
4. **DAISS（tzrKjhxQRsMJ）**：用 NDI leader-follower、双仿真体、多模态示范和阶段感知模仿学习处理探头—针具的非对称双臂运动；依赖专家示范、阶段划分、实时超声/外部视觉和匹配仿真体。它支持“阶段化运动策略改善任务执行”的机制证据，但不是在线刚度更新，也没有力反馈主导的阻抗适应。
5. **MGA（4PMWKo1YeV8J）**：在短时域运动—阻抗序列空间中将 RL 行为先验、模型滚动评估、约束流形修正和接触残差反馈合并；阻抗以 log-stiffness 决策并经矩阵指数保持正定。依赖可用的低层阻抗接口、动力学/接触模型、约束几何和有限滚动预算。它是当前“模型/优化路线”和“学习阻抗路线”的直接混合候选，但全文证据仍来自预印本中的接触任务，尚未证明超声、呼吸运动或本文任务收益。

**对三条路线和缺口的机制层影响：**

- 模型/优化路线不应只写成“根据力估计求参数”：已有证据显示可把局部刚度、回路带宽、时延和接触几何转成速率/阻抗可行域，并可用滚动模型和接触反馈重新评估候选。
- 学习阻抗路线不应只写成“端到端输出阻抗”：当前已读工作分别展示了力—超声融合、阶段调制、阶段化动作和扩散/策略先验；但直接、稳定、可解释的在线阻抗参数更新仍主要由模型/能量/约束层承担。
- 本文缺口应收窄为：在实际在线控制值约束下，把当前接触/刚度/运动信息绑定到可执行的阻抗或运动决策，并用扫描/图像/接触安全等任务收益验证；不能仅以“整段采集更优”或“离线策略精度”替代该因果绑定。
- 仍未完成的机制判断：37 篇人工全文的真实方法和结果、Field Map 的独立 coverage reviewer 结论，以及各路线在同一任务收益和计算预算下的可比性。

## E. 人工下载清单

Controller 的 `human_fulltext_request` 精确包含以下 37 篇；不要下载出版方之外的替代版本，也不要将 arXiv 已读的 5 篇重复提交。每个文件放在 `D:/pip_env/temp/scan_impedance_literature_20261031/aris_run/source-materials/`，建议命名为 `<paper_id>.pdf`。提交 manifest 时必须一次性、每个 `paper_id` 恰好一次，`source_path` 必须是 `source-materials/` 下的相对路径，文件必须存在且非空；`media_type` 可填 `application/pdf`，Controller 会自行计算 SHA-256。

|完整题名|作者|年份|venue|DOI|paper_id|目录/文件名|
|---|---|---:|---|---|---|---|
|A Survey of Autonomous Robotic Ultrasound Scanning Systems|Khushboo Munir; Abdullah F. Al-Battal; Ammar Alsheghri; Harald Becher; Michelle Noga; Kumaradevan Punithakumar|2025|IEEE Access|10.1109/access.2025.3574464|jkUy-uxymaMJ|source-materials/jkUy-uxymaMJ.pdf|
|Toward Clinical Applications of Intelligent Robotic Ultrasound Systems|Taiyu Han; Guochen Ning; Hanying Liang; Zihan Li; Zhongliang Jiang; Fang Chen; Yan Kang; Jianwen Luo; Hongen Liao|2026|IEEE Reviews in Biomedical Engineering|10.1109/rbme.2025.3610605|AdPA6qQA2QMJ|source-materials/AdPA6qQA2QMJ.pdf|
|A Review and Perspective of Techniques for Autonomous Robotic Ultrasound Acquisitions|Yanding Qin; Lele Dang; Fan Ren; Zhuomao Li; Lijun Duan; Hongpeng Wang; Jianda Han|2026|Sensors|10.3390/s26072081|m8UoVKbdM_MJ|source-materials/m8UoVKbdM_MJ.pdf|
|A review of robot-assisted ultrasound examination: Systems and technology|Haiyan Du; Xinran Zhang; Yongde Zhang; Fujun Zhang; Letao Lin; Tao Huang|2024|The International Journal of Medical Robotics and Computer Assisted Surgery|10.1002/rcs.2660|vKkpCit7Dc8J|source-materials/vKkpCit7Dc8J.pdf|
|Robot-Assisted Medical Imaging: A Review|Septimiu E. Salcudean; Hamid Moradi; David G. Black; Nassir Navab|2022|Proceedings of the IEEE|10.1109/jproc.2022.3162840|P4vOsQKTB1oJ|source-materials/P4vOsQKTB1oJ.pdf|
|A shared-control ultrasound-guided robotic system for liver puncture under respiratory motion|Wenlong Zhang; Weihao Lu; Qijie Lu; Wentuo Yang; Le Xie; Yuchen Yang|2026|Scientific Reports|10.1038/s41598-026-56017-3|SS4c7arItfoJ|source-materials/SS4c7arItfoJ.pdf|
|Task Space Compliant Control and Six-Dimensional Force Regulation Toward Automated Robotic Ultrasound Imaging|Junchen Wang; Chunheng Lu; Yifei Lv; Siqin Yang; Mingbo Zhang; Yu Shen|2024|IEEE Transactions on Automation Science and Engineering|10.1109/tase.2023.3282974|f76lmG78c6cJ|source-materials/f76lmG78c6cJ.pdf|
|A Lightweight Robotic Ultrasound System With Hybrid Force Control and Modular Force-Range Switching|Zehao Wu; Weijian Zhang; Xianli Wang; I-Ming Chen; Qingsong Xu|2026|IEEE Transactions on Robotics|10.1109/tro.2026.3734315|wftPH5bwT0gJ|source-materials/wftPH5bwT0gJ.pdf|
|Design of a Cost-Effective Ultrasound Force Sensor and Force Control System for Robotic Extra-Body Ultrasound Imaging|Yixuan Zheng; Hongyuan Ning; Eason Rangarajan; Aban Merali; Adam Geale; Lukas Lindenroth; Zhouyang Xu; Weizhao Wang; Philipp Kruse; Steven Morris; Liang Ye; Xinyi Fu; Kawal Rhode; Richard James Housden|2025|Sensors|10.3390/s25020468|-fQhXbjRKkwJ|source-materials/-fQhXbjRKkwJ.pdf|
|Autonomous Robotic Ultrasound Scanning for Liver Interventions: an Integrated Redundancy Resolution and Motion Control Framework|Peng Chen; Hiroshi Yokoi; Bo Zhang|2026|Intelligent Systems & Robotic Mechanics|10.62836/isrm.v3i1.0001|OWn7BXc3ULAJ|source-materials/OWn7BXc3ULAJ.pdf|
|Learning Variable Impedance Control for Contact Sensitive Tasks|Miroslav Bogdanovic; Majid Khadiv; Ludovic Righetti|2020|IEEE Robotics and Automation Letters|10.1109/lra.2020.3011379|dqYs440MKMMJ|source-materials/dqYs440MKMMJ.pdf|
|On the Explicit Robust Force Control via Disturbance Observer|Emre Sariyildiz; Kouhei Ohnishi|2015|IEEE Transactions on Industrial Electronics|10.1109/tie.2014.2361611|pNSC0fyjM9wJ|source-materials/pNSC0fyjM9wJ.pdf|
|Minimally Model-Based Trajectory Tracking and Variable Impedance Control of Flexible-Joint Robots|Emmanouil Spyrakos-Papastavridis; Jian S. Dai|2021|IEEE Transactions on Industrial Electronics|10.1109/tie.2020.2994886|W_DvVyMz0QAJ|source-materials/W_DvVyMz0QAJ.pdf|
|Tactile-Guided Robotic Ultrasound: Mapping Preplanned Scan Paths for Intercostal Imaging|Yifan Zhang; Dianye Huang; Nassir Navab; Zhongliang Jiang|2025|2025 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS)|10.1109/iros60139.2025.11246299|txpI_b7Ta1IJ|source-materials/txpI_b7Ta1IJ.pdf|
|AuSoScan: Automatic Scoliosis Assessment by Ultrasound Scanning With Soft Contact Control|Anqing Duan; Wanli Liuchen; Peng Zhou; Dezhen Song; Chenguang Yang; Yongping Zheng; David Navarro-Alarcon|2025|IEEE/ASME Transactions on Mechatronics|10.1109/tmech.2025.3583041|QPuQ4sfxSTIJ|source-materials/QPuQ4sfxSTIJ.pdf|
|Ultrasound-Guided Assistive Robots for Scoliosis Assessment With Optimization-Based Control and Variable Impedance|Anqing Duan; Maria Victorova; Jingyuan Zhao; Yuxiang Sun; Yongping Zheng; David Navarro-Alarcon|2022|IEEE Robotics and Automation Letters|10.1109/lra.2022.3186504|c2nzgpQiQ0IJ|source-materials/c2nzgpQiQ0IJ.pdf|
|Inverse-Reinforcement-Learning-Based Robotic Ultrasound Active Compliance Control in Uncertain Environments|Guochen Ning; Hanying Liang; Xinran Zhang; Hongen Liao|2024|IEEE Transactions on Industrial Electronics|10.1109/tie.2023.3250767|YriQFb_sHqEJ|source-materials/YriQFb_sHqEJ.pdf|
|Compliance Control of a Robotic Breast Ultrasound System Based on Variable Admittance Control|Yongde Zhang; Dong Guo; Xuequan Huang; Chuang He|2026|The International Journal of Medical Robotics and Computer Assisted Surgery|10.1002/rcs.70194|jo-MpUbeRIYJ|source-materials/jo-MpUbeRIYJ.pdf|
|Explicit force-based control strategies for compliant robot-assisted ultrasound imaging|Adrian Piedra; R Brooke Jeffrey; Oussama Khatib|2026|The International Journal of Robotics Research|10.1177/02783649261441635|BIP_InW4FSIJ|source-materials/BIP_InW4FSIJ.pdf|
|Dual-Agent-Based Robotic Ultrasound Path Planning and Interaction Control in External-Vision-Independent Environments|Xinye Wang; Zhiyuan He; Peng Chen; Zhe Wang; Tao Sun|2025|IEEE Transactions on Automation Science and Engineering|10.1109/tase.2025.3614996|D8AVZtwhD6gJ|source-materials/D8AVZtwhD6gJ.pdf|
|Breaking the Latency Barrier: Synergistic Perception and Control for High-Frequency 3D Ultrasound Servoing|Yizhao Qian; Yujie Zhu; Jiayuan Luo; Li Liu; Yixuan Yuan; Hongen Liao; Guochen Ning|2026|2026 IEEE International Conference on Robotics and Automation (ICRA)|10.1109/icra57385.2026.11697234|P0yee_3nRaAJ|source-materials/P0yee_3nRaAJ.pdf|
|Impedance Control on Arbitrary Surfaces for Ultrasound Scanning Using Discrete Differential Geometry|Michael Dyck; Arne Sachtler; Julian Klodmann; Alin Albu-Schaffer|2022|IEEE Robotics and Automation Letters|10.1109/lra.2022.3184800|i4wx6VF4kN0J|source-materials/i4wx6VF4kN0J.pdf|
|Optimization-Based Variable Impedance Control of Robotic Manipulator for Medical Contact Tasks|Junling Fu; Ilaria Burzo; Elisa Iovene; Jianzhuang Zhao; Giancarlo Ferrigno; Elena De Momi|2024|IEEE Transactions on Instrumentation and Measurement|10.1109/tim.2024.3372209|GBHeodrYfZQJ|source-materials/GBHeodrYfZQJ.pdf|
|Precise Repositioning of Robotic Ultrasound: Improving Registration-Based Motion Compensation Using Ultrasound Confidence Optimization|Zhongliang Jiang; Nehil Danis; Yuan Bi; Mingchuan Zhou; Markus Kroenke; Thomas Wendler; Nassir Navab|2022|IEEE Transactions on Instrumentation and Measurement|10.1109/tim.2022.3200360|vHYta8tPiJ8J|source-materials/vHYta8tPiJ8J.pdf|
|Dynamic Arm Impedance Regulation: Experiences From Medical Ultrasound Scanning Experts|Baoshan Niu; Dapeng Yang; Tao Geng; Le Zhang; Boyu Ma; Qi Huang; Hong Liu|2024|IEEE Transactions on Instrumentation and Measurement|10.1109/tim.2024.3428622|4SmWT3oeIIgJ|source-materials/4SmWT3oeIIgJ.pdf|
|Safe Contact Force Generation for Robotic Thyroid Ultrasound Imaging|Ryosuke Tsumura; Toshifumi Tomioka; Yoshihiko Koseki; Kiyoshi Yoshinaka|2024|IEEE Robotics and Automation Letters|10.1109/lra.2024.3349959|Ct8Urnw7rssJ|source-materials/Ct8Urnw7rssJ.pdf|
|Smooth Path Planning and Dynamic Contact Force Regulation for Robotic Ultrasound Scanning|Le Zhang; Dapeng Yang; Baoshan Niu; Haonan Yang; Qi Huang; Li Jiang; Hong Liu|2025|IEEE Robotics and Automation Letters|10.1109/lra.2025.3604746|CkE-5mElf_IJ|source-materials/CkE-5mElf_IJ.pdf|
|EMG-based Hybrid Impedance-Force Control for Human-Robot Collaboration on Ultrasound Imaging|Teng Li; Hongjun Xing; Hamid D. Taghirad; Mahdi Tavakoli|2022|2022 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS)|10.1109/iros47612.2022.9981615|oLEp1P_qPKQJ|source-materials/oLEp1P_qPKQJ.pdf|
|A Novel Ultrasound Robot With Force/Torque Measurement and Control for Safe and Efficient Scanning|Xianqiang Bao; Shuangyi Wang; Lingling Zheng; Richard James Housden; Joseph V. Hajnal; Kawal Rhode|2023|IEEE Transactions on Instrumentation and Measurement|10.1109/tim.2023.3239925|bGqa7m5UjxwJ|source-materials/bGqa7m5UjxwJ.pdf|
|Optimizing Variable Admittance Control for Remote Ultrasound Scanning Under Uncertain Environment|Songjie Xiao; Tengyue Wang; Yun Long; Liangjing Yang|2025|IEEE Access|10.1109/access.2025.3568150|memfv6bMNssJ|source-materials/memfv6bMNssJ.pdf|
|A Unified Interaction Control Framework for Safe Robotic Ultrasound Scanning with Human-Intention-Aware Compliance|Xiangjie Yan; Shaqi Luo; Yongpeng Jiang; Mingrui Yu; Chen Chen; Senqiang Zhu; Gao Huang; Shiji Song; Xiang Li|2024|2024 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS)|10.1109/iros58592.2024.10801755|AeZVQ15NBVsJ|source-materials/AeZVQ15NBVsJ.pdf|
|Passivity Preservation for Variable Impedance Control of Compliant Robots|Emmanouil Spyrakos-Papastavridis; Peter R. N. Childs; Jian S. Dai|2020|IEEE/ASME Transactions on Mechatronics|10.1109/tmech.2019.2961478|tZYsUFlRIUIJ|source-materials/tZYsUFlRIUIJ.pdf|
|Simultaneous Online Registration-Independent Stiffness Identification and Tip Localization of Surgical Instruments in Robot-Assisted Eye Surgery|Ali Ebrahimi; Shahriar Sefati; Peter Gehlbach; Russell H. Taylor; Iulian I. Iordachita|2023|IEEE Transactions on Robotics|10.1109/tro.2022.3201393|HokaNUBl3lAJ|source-materials/HokaNUBl3lAJ.pdf|
|Diffusion-Based Impedance Learning for Contact-Rich Manipulation Tasks|Noah Geiger; Tamim Asfour; Neville Hogan; Johannes Lachner|2026|IEEE Transactions on Robotics|10.1109/tro.2026.3717255|kexwZicrUfAJ|source-materials/kexwZicrUfAJ.pdf|
|Quasi-Static Elastography-Driven Automated Robotic Ultrasound Screening and Localization|Hanying Liang; Shipeng Zhang; Guochen Ning; Jianwen Luo; Hongen Liao|2026|IEEE Transactions on Biomedical Engineering|10.1109/tbme.2025.3634552|Vzv82t4dJgsJ|source-materials/Vzv82t4dJgsJ.pdf|
|Acoustic Shadowing Aware Robotic Ultrasound: Lighting up the Dark|Viviana Sutedjo; Maria Tirindelli; Christine Eilers; Walter Simson; Benjamin Busam; Nassir Navab|2022|IEEE Robotics and Automation Letters|10.1109/lra.2022.3141451|11nhTera-oQJ|source-materials/11nhTera-oQJ.pdf|
|A Passive Variable Impedance Control Strategy with Viscoelastic Parameters Estimation of Soft Tissues for Safe Ultrasonography|Luca Beber; Edoardo Lamon; Davide Nardi; Daniele Fontanelli; Matteo Saveriano; Luigi Palopoli|2024|2024 IEEE International Conference on Robotics and Automation (ICRA)|10.1109/icra57147.2024.10610167|ODwWwSVciuYJ|source-materials/ODwWwSVciuYJ.pdf|

提交 manifest（manifest JSON 必须列出以上 37 个 `paper_id` 各一次）：

```powershell
python -X utf8 -m arisctl --root "D:/pip_env/temp/scan_impedance_literature_20261031/aris_run" submit-human-fulltext-batch recent-related-work-20261031 "D:/pip_env/temp/scan_impedance_literature_20261031/aris_run/human_fulltext_manifest.json"
```

manifest 每项最小形式为：

```json
{
  "papers": [
    {"paper_id": "jkUy-uxymaMJ", "source_path": "source-materials/jkUy-uxymaMJ.pdf", "media_type": "application/pdf"}
  ]
}
```

不要在人工 manifest 中加入 5 篇已读 arXiv 论文；提交后 Controller 才会回到 `PAPER_READING`，不会自动制造 Evidence。
