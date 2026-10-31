# 本项目研究现状的近邻与同期工作补查

日期：2026-10-31。effort：balanced。范围：2015–2026；按用户要求额外核读 F10（2014）。输出为中文临时调研备忘录；**没有修改项目文件、既有文献索引或正式 ARIS 状态，也未提交/推送。** 本轮用户“不要修改项目文件”的要求优先，因此没有运行会更新项目工作树的 pull 或 research_sync。

采用用户指定的 ARIS `research-lit`、`novelty-check`，以及 `openalex`、`semantic-scholar` 子技能。属于 **ad-hoc combined**（问题与方法分别做初步风险筛查），不是 Controller 管理的正式 landscape/novelty Gate。用户明确允许 OpenAlex/S2，按此检索；未把普通网页搜索冒称 Google Scholar。网页发现与摘要比较均为暂定判断，不能直接晋升正式 ARIS Evidence。

已读项目 `research/paper_draft_zh.md` 第一、二节、`research/field_state.md`、`literature/literature_index.md`，并查看交接与文献阶段要求。N01–N04、M01–M03、C10、S01/S03/S04、G13 作为既有核读锚点；不把重新检出的锚点计成新增发现。项目中 `wijayarathne.pdf` 实际是 **T-RO2023**，不能因作者相同误记为 N04 ICRA2020。

**先说结论。** 文献已覆盖“在线刚度更新＋未知生理运动力补偿”“无视觉的力反馈预测控制”和“IMM 生理运动预测”这三层能力。当前草稿不能把它们的宽泛组合写成尚未解决。较有希望保留、但仍需最接近原文验证的增量是：**在共用探头力觉与位置接口下，限定漂移子空间和有效激励后识别局部刚度；区分运动量测与反馈所用刚度；处理量测误差、伺服时延和增益更新对闭环的耦合。** 本轮没有找到足够全文证据证明这一完整具体机制已被等价覆盖，也没有证据认证其新颖。

## A. 候选清单

类别：①在线刚度＋环境运动/表面位置估计；②力觉生理运动估计/预测补偿；③2023–2026超声运动下力控；④AOB/Kalman/IMM/多模型；⑤变刚度增益调度/归一化与稳定性。多类可并列。“威胁新颖性”指需要正面对照的候选风险，**不代表已经由全文确认覆盖本文**。人机协作的意图/阻抗工作只作外域背景，不放大为超声直接先例。

证据级别按实际阅读区分：“摘要”是数据库或原始出版/机构页面上的实际摘要；已有本地全文、此前已核读与本轮读全文分别说明。只有 F10 本轮完成整篇核读。未读全文者，摘要未提及某能力不等于该能力不存在。表内 OpenAlex 链接是身份来源；机制解释另引原始摘要或本地原文。

| 编号 | 题目 | 作者 | 年份 | venue | DOI | 类别 | 与本文关系 | 是否已在本地（本项目） | 证据级别 | 来源 |
|---|---|---|---|---|---|---|---|---|---|---|
| L01 | Viscoelastic model based force control for soft tissue interaction and its application in physiological motion compensation | Pedro Moreira、Nabil Zemiti、Chao Liu、Philippe Poignet | 2014（指定例外） | Computer Methods and Programs in Biomedicine | [10.1016/j.cmpb.2014.01.017](https://doi.org/10.1016/j.cmpb.2014.01.017) | 1/2/4/5 | 威胁新颖性：已有在线刚度＋未知运动补偿 | 是，F10 | 已读全文，本轮PDF核读 | [OpenAlex](https://openalex.org/W2019057824) |
| L02 | Robot Force Control on a Beating Heart | Rui Cortesão、Michel Dominici | 2017 | IEEE/ASME Transactions on Mechatronics | [10.1109/TMECH.2017.2696259](https://doi.org/10.1109/TMECH.2017.2696259) | 2/4 | 威胁新颖性：双AOB力控及自主运动补偿 | 否 | 摘要；全文待核 | [OpenAlex](https://openalex.org/W2609790749) |
| L03 | Correction to “Robot Force Control on a Beating Heart” | Rui Cortesão、Michel Dominici | 2017 | IEEE/ASME Transactions on Mechatronics | [10.1109/TMECH.2017.2765798](https://doi.org/10.1109/TMECH.2017.2765798) | 4 | 背景：L02勘误，必须随原文核查 | 否 | 摘要/勘误元数据 | [OpenAlex](https://openalex.org/W2772547471) |
| L04 | Cascade force control for autonomous beating heart motion compensation | Michel Dominici、Rui Cortesão | 2015 | Control Engineering Practice | [10.1016/j.conengprac.2014.12.012](https://doi.org/10.1016/j.conengprac.2014.12.012) | 2/4 | 威胁新颖性：力反馈、无视觉、MPC–AOB | 否 | 出版方实际摘要＋highlights | [OpenAlex](https://openalex.org/W2049050176) |
| L05 | Towards physiological motion compensation for flexible needle interventions | Pedro Moreira、Momen Abayazid、Sarthak Misra | 2015 | IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS) | [10.1109/IROS.2015.7353468](https://doi.org/10.1109/IROS.2015.7353468) | 2/4 | 威胁新颖性：从接触力估生理运动 | 否 | 摘要；KF细节未核 | [OpenAlex](https://openalex.org/W2199871332) |
| L06 | Robot-Assisted Ultrasound Probe Force Control Under Respiration-Induced Motion | Kai-Shiang Yuan、Yu‐Hsiu Lee | 2024 | IEEE International Conference on Advanced Intelligent Mechatronics (AIM) | [10.1109/AIM55361.2024.10637095](https://doi.org/10.1109/AIM55361.2024.10637095) | 2/3 | 威胁新颖性：超声呼吸预测补偿直接近邻 | 否 | 摘要；信号和接口待核 | [OpenAlex](https://openalex.org/W4401750884) |
| L07 | Beating Heart Motion Accurate Prediction Method Based on Interactive Multiple Model: An Information Fusion Approach | Fan Liang、Weihong Xie、Yang Yu | 2017 | BioMed Research International | [10.1155/2017/1279486](https://doi.org/10.1155/2017/1279486) | 4 | 威胁新颖性：IMM生理运动预测已有先例 | 否 | 摘要；并非全文 | [OpenAlex](https://openalex.org/W2764129983) |
| L08 | Impedance estimation for robot contact with uncalibrated environments | Wenrui Wang、Qinwen Li、Chenghua Lu 等 | 2021 | Mechanical Systems and Signal Processing | [10.1016/j.ymssp.2021.107819](https://doi.org/10.1016/j.ymssp.2021.107819) | 1 | 威胁新颖性：联合阻抗与环境位置估计 | 否 | 出版方摘要；持续运动能力待核 | [OpenAlex](https://openalex.org/W3145470387) |
| L09 | An interaction controller formulation to systematically avoid force overshoots through impedance shaping method with compliant robot base | Loris Roveda、Nicola Pedrocchi、Federico Vicentini 等 | 2016 | Mechatronics | [10.1016/j.mechatronics.2016.08.001](https://doi.org/10.1016/j.mechatronics.2016.08.001) | 1/4/5 | 近邻：EKF刚度＋KF底座运动，改变的是机器人侧 | 否 | 出版方摘要 | [OpenAlex](https://openalex.org/W2512251902) |
| L10 | Hunt–Crossley model based force control for minimally invasive robotic surgery | A. Pappalardo、A. Albakri、C. Liu 等 | 2016 | Biomedical Signal Processing and Control | [10.1016/j.bspc.2016.05.003](https://doi.org/10.1016/j.bspc.2016.05.003) | 4/5 | 近邻：HC黏弹模型＋Kalman AOB力控 | 否 | 出版方摘要 | [OpenAlex](https://openalex.org/W2399731623) |
| L11 | Flexible Needle Steering in Moving Biological Tissue With Motion Compensation Using Ultrasound and Force Feedback | Jason Chevrie、Navid Shahriari、Marie Babel、Alexandre Krupa、Sarthak Misra | 2018 | IEEE Robotics and Automation Letters | [10.1109/LRA.2018.2809484](https://doi.org/10.1109/LRA.2018.2809484) | 2/4 | 近邻：力反馈运动补偿与US/EM针跟踪融合 | 否 | 摘要 | [OpenAlex](https://openalex.org/W2788264179) |
| L12 | Bayesian Estimation of Human Impedance and Motion Intention for Human–Robot Collaboration | Xinbo Yu、Wei He、Yanan Li 等 | 2021刊期（OA记2019） | IEEE Transactions on Cybernetics | [10.1109/TCYB.2019.2940276](https://doi.org/10.1109/TCYB.2019.2940276) | 1/4 | 背景：人臂阻抗/意图联合Bayesian估计，非扫描直接先例 | 否 | 摘要 | [OpenAlex](https://openalex.org/W2980410359) |
| L13 | Simultaneous Estimation of Human Motion Intention and Time-Varying Arm Stiffness for Enhanced Human–Robot Interaction | Huayang Wu、Chengzhi Zhu、Long Cheng 等 | 2025刊期/2024在线 | IEEE Transactions on Cognitive and Developmental Systems | [10.1109/TCDS.2024.3480854](https://doi.org/10.1109/TCDS.2024.3480854) | 1/4 | 背景：握力刚度先验＋LQE联合估意图 | 否 | 摘要 | [OpenAlex](https://openalex.org/W4403446657) |
| L14 | Respiratory Motion-Robust Robotic Ultrasound Acquisitions via Vision-Haptic Fusion Control and 3-D Compensation | Guanglin Cao、Danqian Cao、Hongbin Liu | 2025 | IEEE Transactions on Instrumentation and Measurement | [10.1109/TIM.2025.3612625](https://doi.org/10.1109/TIM.2025.3612625) | 3 | 近邻：已有独立视觉预测＋导纳参考补偿 | 是，G13 | 本轮摘要；项目已有原文核读 | [OpenAlex](https://openalex.org/W4414603892) |
| L15 | Robotic Respiratory Liver Motion Compensation for Stable Ultrasound Imaging of Small Lesions | Long Lei、Yingying Han、Jinhang Li 等 | 2025 | IEEE Sensors Journal | [10.1109/JSEN.2025.3598958](https://doi.org/10.1109/JSEN.2025.3598958) | 3 | 近邻：US图像ReME-Net估呼吸运动＋变导纳 | 否 | 摘要 | [OpenAlex](https://openalex.org/W4413371345) |
| L16 | Safe Variable-Force Control for Robotic Ultrasound Scanning via ℓ 2 -Regularized Fuzzy Broad Offline Reinforcement Learning | Liang Yao、Yun Zhang、Yeying Qin 等 | 2026 | IEEE Transactions on Medical Robotics and Bionics | [10.1109/TMRB.2026.3699233](https://doi.org/10.1109/TMRB.2026.3699233) | 3 | 近邻：在体呼吸条件下学习阻抗适应率 | 否 | 摘要；尚不确认量测与接口 | [OpenAlex](https://openalex.org/W7163132243) |
| L17 | Motion frequency exploration based force separator for surgical robots interacting with a beating heart | Yanran Wei、Wenshuo Li、Jiayin Wang 等 | 2024 | Neurocomputing | [10.1016/j.neucom.2024.128079](https://doi.org/10.1016/j.neucom.2024.128079) | 2/4 | 近邻：分离心跳运动耦合力与机器人状态耦合力 | 否 | 出版方摘要；不等于估局部组织刚度 | [OpenAlex](https://openalex.org/W4399771792) |
| L18 | Switched-Impedance Control of Surgical Robots in Teleoperated Beating-Heart Surgery | Lingbo Cheng、Mahdi Tavakoli | 2018 | Journal of Medical Robotics Research | [10.1142/S2424905X18410039](https://doi.org/10.1142/S2424905X18410039) | 2/4/5 | 近邻：接触/非接触切换阻抗与心跳补偿 | 否 | 摘要 | [OpenAlex](https://openalex.org/W2790279171) |
| L19 | Stability Considerations for Variable Impedance Control | Klas Kronander、Aude Billard | 2016 | IEEE Transactions on Robotics | [10.1109/TRO.2016.2593492](https://doi.org/10.1109/TRO.2016.2593492) | 5 | 背景：时变K/D充分稳定条件 | 是，H04 | 本轮摘要＋§II–III原文复核 | [OpenAlex](https://openalex.org/W2517641405) |
| L20 | Passivity Preservation for Variable Impedance Control of Compliant Robots | Emmanouil Spyrakos-Papastavridis、Peter Childs、Jian S. Dai | 2020刊期/2019在线 | IEEE/ASME Transactions on Mechatronics | [10.1109/TMECH.2019.2961478](https://doi.org/10.1109/TMECH.2019.2961478) | 5 | 背景：无源保持；须连读L33–L34争论 | 否 | 摘要 | [OpenAlex](https://openalex.org/W2997161337) |
| L21 | Implicit force control for an industrial robot based on stiffness estimation and compensation during motion | Roberto Rossi、Luca Fossali、Alberto Novazzi 等 | 2016 | Proceedings - IEEE International Conference on Robotics and Automation/Proceedings | [10.1109/ICRA.2016.7487243](https://doi.org/10.1109/ICRA.2016.7487243) | 1/5 | 威胁新颖性：位置接口、运动中刚度辨识与增益缩放 | 是，N01，manual_inbox | 既有全文核读；本轮局部原文复核 | [OpenAlex](https://openalex.org/W2416622902) |
| L22 | Identification of Compliant Contact Parameters and Admittance Force Modulation on a Non-stationary Compliant Surface | Lasitha Wijayarathne、Frank L. Hammond | 2020 | Proceedings - IEEE International Conference on Robotics and Automation/Proceedings | [10.1109/ICRA40945.2020.9196897](https://doi.org/10.1109/ICRA40945.2020.9196897) | 1/5 | 威胁新颖性：机械探测＋非静止环境位置型自适应力控 | 既有N04核读记录；本轮项目PDF目录未定位 | 本轮摘要；此前arXiv核读记录 | [OpenAlex](https://openalex.org/W3019009855) |
| L23 | Computed-Torque Control for Robotic-Assisted Tele-Echography Based on Perceived Stiffness Estimation | Luís Santos、Rui Cortesão | 2018 | IEEE Transactions on Automation Science and Engineering | [10.1109/TASE.2018.2790900](https://doi.org/10.1109/TASE.2018.2790900) | 1/4/5 | 威胁新颖性：感知刚度更新力控增益/观察器 | 是，manual_inbox | 本轮摘要；项目已有原文核读 | [OpenAlex](https://openalex.org/W2793044430) |
| L24 | Explicit force-based control strategies for compliant robot-assisted ultrasound imaging | Adrian Piedra、R. Brooke Jeffrey、Oussama Khatib | 2026 | The International Journal of Robotics Research | [10.1177/02783649261441635](https://doi.org/10.1177/02783649261441635) | 3 | 近邻：力矩接口、移动模体显式力/力矩控制 | 是，C10 | 本轮摘要；项目已有原文核读 | [OpenAlex](https://openalex.org/W7160298333) |
| L25 | Observer-based sliding mode admittance control for ultrasound robot force-tracking in complex interaction environments | Sen Li 等（完整作者以原文为准） | 2026 | Control Engineering Practice | [10.1016/j.conengprac.2026.106764](https://doi.org/10.1016/j.conengprac.2026.106764) | 3/4 | 近邻：DOSMAC集总扰动估计与滑模导纳 | 是，S01 | 本轮摘要；项目已有原文核读 | [OpenAlex](https://openalex.org/W7123713864) |
| L26 | Autonomous robotic ultrasound scanning on complex surfaces: Automatic trajectory generation and disturbance observer-based model predictive admittance control | Jia'en Shen、Yongfei Feng、Qiqi Pan 等 | 2026在线/2027卷期 | Control Engineering Practice | [10.1016/j.conengprac.2026.107240](https://doi.org/10.1016/j.conengprac.2026.107240) | 3/4 | 近邻：DO-MPAC，ESO＋预测虚拟动作 | 是，S03 | 本轮首页摘要；项目已有原文核读 | [OpenAlex](https://openalex.org/W7208705925) |
| L27 | Trajectory-aware feedforward admittance control for contact force tracking with time-varying reference trajectory in a teleoperated robotic ultrasound system | Sen Li、Yufan Pi、Weihua Li 等 | 2026刊期/2025在线 | Control Engineering Practice | [10.1016/j.conengprac.2025.106674](https://doi.org/10.1016/j.conengprac.2025.106674) | 3/5 | 近邻：已知参考速度/加速度前馈与时延补偿 | 是，S04 | 本轮摘要；项目已有原文核读 | [OpenAlex](https://openalex.org/W4417070314) |
| L28 | A Passive Variable Impedance Control Strategy with Viscoelastic Parameters Estimation of Soft Tissues for Safe Ultrasonography | Luca Beber、Edoardo Lamon、Davide Nardi 等 | 2024 | Proceedings - IEEE International Conference on Robotics and Automation/Proceedings | [10.1109/ICRA57147.2024.10610167](https://doi.org/10.1109/ICRA57147.2024.10610167) | 1/3/5 | 近邻：在线黏弹参数＋QP变阻抗＋能量罐 | 是，beber.pdf | 本轮摘要；项目已有核读 | [OpenAlex](https://openalex.org/W4401415597) |
| L29 | Contact-Aware Impedance Controller for Robot-Assisted Ultrasound Imaging | Md Miraj Arefin、M Efe Tiryaki | 2026-09-28 | arXiv (Cornell University) | [10.48550/arXiv.2609.34400](https://doi.org/10.48550/arXiv.2609.34400) | 3 | 近邻：US接触定位/控制坐标更新，预印本 | 否 | 原站实际摘要；helper verify_pending | [OpenAlex](https://openalex.org/W7214772732) |
| L30 | Dual-Rate Force-Image Control with Model-Based Orientation Limits for Robotic Ultrasound | Tyler Foster、Qiang Zhang、A B M Tahidul Haque 等 | 2026-10-05 | arXiv (Cornell University) | [10.48550/arXiv.2610.05839](https://doi.org/10.48550/arXiv.2610.05839) | 3/5 | 近邻：双速率力/图像与模型限速，预印本 | 否 | 原站实际摘要；helper verify_pending | [OpenAlex](https://openalex.org/W7220742231) |
| L31 | Probabilistic Approach to Online Stiffness Estimation for Robotic Tasks | Toshiaki Tsuji、Tsukasa Kusakabe | 2022 | IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS) | [10.1109/IROS47612.2022.9981146](https://doi.org/10.1109/IROS47612.2022.9981146) | 1/4 | 近邻：概率可靠度门控刚度辨识 | 否 | 摘要 | [OpenAlex](https://openalex.org/W4312443790) |
| L32 | Stabilizing Switching Force Control for the Hunt-Crossley Model | Yun-jae Yang、Jeong-Un Nam、Tae-Woong Yoon | 2024 | International Journal of Control Automation and Systems | [10.1007/s12555-023-0426-3](https://doi.org/10.1007/s12555-023-0426-3) | 5 | 背景：HC接触/非接触切换稳定性 | 否 | 出版方/机构库摘要 | [OpenAlex](https://openalex.org/W4396779937) |
| L33 | Counterexample to Passivity Preservation for Variable Impedance Control of Compliant Robots | Gustavo A. Medrano‐Cerda | 2023刊期（OA记2020） | IEEE/ASME Transactions on Mechatronics | [10.1109/TMECH.2020.3019651](https://doi.org/10.1109/TMECH.2020.3019651) | 5 | 背景：针对L20理论的反例，争论未由本轮裁决 | 否 | 摘要 | [OpenAlex](https://openalex.org/W3080632296) |
| L34 | Confutation of the “Counterexample to Passivity Preservation for Variable Impedance Control of Compliant Robots” | Emmanouil Spyrakos-Papastavridis、Peter Childs、Jian S. Dai | 2024 | IEEE/ASME Transactions on Mechatronics | [10.1109/TMECH.2024.3382261](https://doi.org/10.1109/TMECH.2024.3382261) | 5 | 背景：对L33反例的作者反驳 | 否 | 摘要 | [OpenAlex](https://openalex.org/W4394938941) |
| L35 | Towards Robotised Palpation for Cancer Detection through Online Tissue Viscoelastic Characterisation with a Collaborative Robotic Arm | Luca Beber、Edoardo Lamon、Giacomo Moretti 等 | 2024 | IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS) | [10.1109/IROS58592.2024.10802582](https://doi.org/10.1109/IROS58592.2024.10802582) | 1/4 | 威胁新颖性：EKF联合压入量和黏弹性质 | 否；本轮arXiv下载未完成 | 摘要；全文待核 | [OpenAlex](https://openalex.org/W4405787522) |


出版信息注意：L12、L13、L20、L26、L27、L33 的在线年与卷期年不同，不能用 DOI 中的年份机械决定发表年份。L26 本地PDF首页为 CEP178(2027)，DOI及在线记录为2026，因此作为2026已可见同期工作保留。L29/L30仅为预印本；L29自述“submitted to ICRA2027”，并非已被ICRA2027接收。二者 arXiv 标识与题名已在原站核对，ARIS DOI helper 对其 `verify_pending`，不伪装成 Crossref 已核验的正式期刊论文。

### 五类覆盖的实质判断

| 检索重点 | 这次补到的主要能力 | 仍不能据此断言 |
|---|---|---|
| ① 刚度与运动/位置 | F10在线β与未知运动补偿；L08联合阻抗/接触位置；L35联合压入量/黏弹性质；L09估刚度与机器人底座运动 | L08/L35是否在未知持续运动中实现可靠分离，需全文；“估压入量”不能自动写成“估自由表面运动” |
| ② 以力估生理运动 | CEP2015无视觉MPC–AOB；TMECH2017双AOB；IROS2015接触力估运动；AIM2024呼吸预测补偿 | 不能再说运动预测只能依赖视觉/声纳/独立位移传感器；是否与同一扫描触点共用测力须逐篇查 |
| ③ 2023–2026超声 | AIM2024、Lei2025图像运动＋变导纳、FLAIC2026在体呼吸力控；保留G13/C10/S01/S03/S04与2026秋季预印本 | 图像配准/重建改善不等于力控改善；外部无跟踪设备不等于只用力与编码器 |
| ④ AOB/KF/IMM | Moreira/Cortesão/Dominici已有力觉补偿；Liang2017已用IMM区分常态和异常心动；Tsuji2022可靠度筛选刚度数据 | 不能把Kalman/IMM本身、两个过程噪声或“可信才更新”作为首次机制 |
| ⑤ 时变/切换稳定性 | Rossi运动中刚度误差鲁棒性；H04的时变K/D条件；Beber能量罐；HC接触/非接触稳定切换；PPC争论链 | 冻结刚度稳定不推出时变增益、延迟预测和恢复切换的完整稳定；不能单凭争论任一方摘要裁决 |

L17 的“分离”对象是运动耦合力与机器人动态状态耦合力，不是本文的刚度与组织运动。L09 分离的运动来自机器人柔顺底座，不是患者侧；这两篇是机制核查线索，不能用名称相似宣布等价。L11 同时使用US目标与EM针追踪，不能描述为整个系统无额外传感器。L15 的“without tracking device”仍使用US图像估运动。[L09出版方摘要](https://www.sciencedirect.com/science/article/pii/S0957415816300769)、[L17出版方摘要](https://www.sciencedirect.com/science/article/pii/S0925231224008506)、[L11记录](https://openalex.org/W2788264179)、[L15记录](https://openalex.org/W4413371345)。

## B. 优先补引/全文核查的8篇及具体区别

以下是**相对当前正文的补引优先级**，不表示应降低已读 Rossi/N01、Wijayarathne/N04、TASE2018 的比较优先级。既有三篇仍是位置接口、主动辨识与刚度缩放的直接基准；它们不需要作为“新文献”重复下载。

### 1. F10：Moreira等，CMPB2014（L01）——必须进入现状主链

原文已经把在线β-KF更新、AOB状态/误差补偿、增益更新和未知运动下的力控连接起来。其静态与运动实验并非只用固定离线刚度。它直接压缩“已有方法只取刚度/运动之一”的宽泛缺口表述。

与本文的可核区别：F10以黏弹状态反馈和active state处理合并误差，未给针对运动子空间的不变刚度估计或显式等效运动量测；采用动力学逆/转矩执行。本文用局部漂移投影、双用途刚度和位置指令预测前馈。**区别是否带来任务收益仍待验证，接口不同本身不是算法创新。** 详细原文证据见C。无需下载。

### 2. Cortesão与Dominici，TMECH2017＋勘误（L02/L03）——最重要的AOB近邻

摘要明确提出双AOB：一支调期望交互力，一支自主补偿心脏运动；离体心脏三维实验给出原始数据RMSE低于0.56 N。足以否定“力控与运动估计/补偿在既有文献中始终分开”的笼统判断。[OpenAlex实际摘要](https://openalex.org/W2609790749)、[S2身份交叉核验](https://www.semanticscholar.org/paper/56ae5ddb9e70c5ecaf6d5259cfcf54858bb5ce92)。

尚待核清：第二AOB实际估计的是环境位置/速度还是等效扰动；参数如何适应；执行接口；勘误改了什么。**未读全文，不能给它贴“固定刚度”“无延迟补偿”的标签。** 本文需要比较的是短窗投影和刚度标定量测相对于双AOB的残余改进，而非“多了运动观察器”。需用户下载原文及勘误。

### 3. Dominici与Cortesão，CEP2015（L04）——力反馈预测参考已有实现

出版方明确：无需视觉或心脏/呼吸运动先验；内层Kalman AOB作模型参考自适应，外层MPC生成心动补偿的控制参考。这里“预测＋补偿参考”已通过接触力反馈实现，不能继续把运动预测前馈总体归为依赖独立运动传感器的一类。[原始摘要与highlights](https://www.sciencedirect.com/science/article/abs/pii/S0967066114002858)。

与本文的潜在区别是：其级联参考/扰动架构与本文的显式刚度标定等效运动量测、局部漂移不变性和保守/准确刚度分工是否不同。没有全文，尚不能认定它未处理刚度变化，也不能把MPC的预测对象自动等同于患者物理位移。需用户下载。

### 4. Yuan与Lee，AIM2024（L06）——同任务同期预测补偿

六自由度臂、六轴F/T与超声探头；以导纳模式定位，以利用呼吸准周期性的自适应模式预测并消除扰动，在海绵腹部模体上优于基本导纳。任务、接触传感和预测补偿三方面都很近。[作者机构库实际摘要](https://scholars.lib.ntu.edu.tw/entities/publication/8f89fceb-12b6-4cd1-8fa4-920732e7a913)。

本文针对呼吸与突发推拉、短窗刚度辨识及跨刚度响应；但**仅凭摘要不能确认对方只用力/编码器、是否估刚度、采用何种接口或是否只适合周期运动**。摘要支持“利用准周期性”，不支持“非周期时必然失效”。需用户下载；这是优先于新增普通静态力控论文的全文缺口。

### 5. Moreira、Abayazid与Misra，IROS2015（L05）——力觉可以成为生理运动观测

摘要明确用保持接触的力传感器估计组织运动，并用于柔性针运动补偿；肺、肾、心脏三类扰动的15次试验，平均靶向误差1.05±0.41 mm。它是“从力获得生理运动信息”必须引用的先例。[OpenAlex实际摘要](https://openalex.org/W2199871332)、[会议官方摘要，PDF内检索标题](https://tams.informatik.uni-hamburg.de/publications/2015/IROS_Digest_WWW.pdf)。

已知与本文的区别：该力触点用于针介入的运动观测，与“同一扫描探头边调力边辨识/观测运动”有不同信息结构。KF如何用刚度换算、是否需已知接触位置/预标定、针/靶追踪传感如何参与，需正文确认。**本轮只读摘要，不把检索片段中的KF或标定细节写成已核机制。** 需用户下载。

### 6. Liang、Xie与Yu，BioMed Research International2017（L07）——IMM必须承认既有应用

实际摘要已经说明：正常心动与心律异常分别对应不同动力学模式，IMM作自适应模式融合，并按信号质量修改转换概率，使用四个数据集评价。故“首次把IMM引入生理运动预测”“两模式应对平静/异常变化”都不是可保留的一般创新。[实际摘要](https://openalex.org/W2764129983)。

其摘要并未建立力觉观测、在线刚度辨识或接触闭环；本文的可讨论差异是**刚度依赖的等效运动观测及其控制误差结构**，而不是IMM命名或模块出现。具体输入、预测时域、模式方程仍需全文，不用其摘要认证本文具体IMM实现不同。需用户下载，虽然它有公开全文，也按本任务集中交接。

### 7. Beber等，IROS2024（L35）——不能漏掉联合压入量/参数的EKF

这篇不同于本地的 **Beber ICRA2024 能量罐变阻抗论文**。摘要报告：协作臂触诊中，以降维Hunt–Crossley接触模型构造EKF，在线联合估计压入量与组织黏弹性质；可从机器人控制器数据工作而不用独立F/T，实验包含硬夹杂硅胶。[OpenAlex实际摘要](https://openalex.org/W4405787522)。

这对“位移和刚度只能各管一半”的叙述是直接风险。压入量是相对形变，不等于自由表面位置；尚未核清它是否允许组织底面未知运动、如何可观、是否有主动激励，也不能把“无独立力传感器”理解为只用编码器。本文需要展示对运动子空间消除后的刚度识别条件，而不能只说“同时估刚度和位移”。arXiv2404.09542下载尝试未完成，本轮仍为摘要级，需用户下载。

### 8. Yao等，TMRB2026（L16）——最新在体呼吸条件的强任务近邻

FLAIC以有限示范训练离线RL，在线调节阻抗适应率；摘要报告仿真、模体及呼吸下人体上肢静脉/腹部扫描，改善过冲与跟踪，无需逐患者重新调参。它不只是静态表面演示。[实际摘要](https://openalex.org/W7163132243)。

本文是显式可解释的辨识—量测—预测路线；对方从示范学适应率。这是目前摘要支持的机制差别。全文未核：状态输入是否含额外观测、如何处理呼吸、是否显式估刚度/运动、位置还是力矩接口。不能因对方用了学习就断言它需要更多在线传感，也不能用本项目仿真宣称优于其在体能力。需用户下载。

### 紧随其后的核查：联合表面位置、稳定性与近期预印本

**L08（Wang MSSP2021）**用L-BFGS联合估环境阻抗和未标定接触位置，直接约束“未知表面位置＋刚度联合估计”的广泛首创表述；正文需核查表面是否持续运动、时间窗假设及唯一性条件。[出版方摘要](https://www.sciencedirect.com/science/article/abs/pii/S0888327021002144)。**L31（Tsuji IROS2022）**已有概率噪声预测、只用可靠数据更新刚度的思想；本项目的“可信窗口才更新”应与“对指定漂移代数不变”区分。[实际摘要](https://openalex.org/W4312443790)。

**L19（H04）**原文§III、Theorem1针对恒定目标惯量H、对称正定且连续可微K(t)、D(t)，给出存在α>0使 αH−D(t)及 K̇(t)+αḊ(t)−2αK(t)负半定的条件；严格负定进一步给渐近性（刊p.1300）。因此本项目同时缩放导纳M/D/K、含位置伺服滞后和接触恢复时，不能原样套用。**L21（Rossi）**已经给运动中估计误差鲁棒性，必须保留为直接理论比较；投影代替带通的差异要落到不变条件与残差界，而不是“我们考虑了估计误差”。

**L32**的摘要将接触内PID与非接触位置/速度控制分别构造Lyapunov函数，再验证切换所需下降关系，是接触恢复完整分析的参考；不等于本文混合系统自动稳定。[作者机构摘要](https://pure.korea.ac.kr/en/publications/stabilizing-switching-force-control-for-the-hunt-crossley-model/)。**L20/L33/L34**必须作为“原方法—反例—作者反驳”一组对待；本轮没有核读证明，既不把PPC当成无争议保证，也不把反例标题当作已推翻其全部结果。[原论文机构记录](https://kclpure.kcl.ac.uk/portal/en/publications/passivity-preservation-for-variable-impedance-control-of-complian/)、[反例身份/摘要](https://openalex.org/W3080632296)、[作者反驳机构记录](https://kclpure.kcl.ac.uk/portal/en/publications/confutation-of-the-counterexample-to-passivity-preservation-for-v/)。

2026秋季的 **L29** 是US图像驱动的有效接触点/坐标系更新，**L30** 是局部刚度进入旋转速率界和双速率力/图像控制；都不是已证实的“力/编码器联合刚度—组织运动分离”。L30还在摘要中报告：所测时延下重建收益不显著，固定限速与模型上限的力误差无可分辨差异。这提醒本项目预测与更新收益须由对照证明，不能仅因增加预测模块就推定增益。[L29原站摘要](https://arxiv.org/abs/2609.34400)、[L30原站摘要](https://arxiv.org/abs/2610.05839)。

### 需用户下载清单（集中一次，不含已在本地的F10等）

第一批为上述8项中的7项缺原文＋1篇配套勘误；F10无需下载。

| 优先 | 编号 | 需下载原文 | DOI | 原文要解决的具体疑问 |
|---|---|---|---|---|
| 1 | L02 | Robot Force Control on a Beating Heart | 10.1109/TMECH.2017.2696259 | 双AOB观测对象、刚度更新与接口 |
| 1配套 | L03 | Correction to “Robot Force Control on a Beating Heart” | 10.1109/TMECH.2017.2765798 | 纠正内容是否影响机制/稳定性 |
| 2 | L04 | Cascade force control for autonomous beating heart motion compensation | 10.1016/j.conengprac.2014.12.012 | 由力到运动参考的实际状态和参数依赖 |
| 3 | L06 | Robot-Assisted Ultrasound Probe Force Control Under Respiration-Induced Motion | 10.1109/AIM55361.2024.10637095 | 输入、准周期预测机制、刚度/接口 |
| 4 | L35 | Towards Robotised Palpation for Cancer Detection through Online Tissue Viscoelastic Characterisation with a Collaborative Robotic Arm | 10.1109/IROS58592.2024.10802582；arXiv2404.09542 | 相对压入量和材料参数如何联合可观，运动边界 |
| 5 | L05 | Towards physiological motion compensation for flexible needle interventions | 10.1109/IROS.2015.7353468 | 力运动换算、标定与额外触点 |
| 6 | L07 | Beating Heart Motion Accurate Prediction Method Based on Interactive Multiple Model: An Information Fusion Approach | 10.1155/2017/1279486 | 观测来源、模式/噪声结构、先验 |
| 7 | L16 | Safe Variable-Force Control for Robotic Ultrasound Scanning via ℓ2-Regularized Fuzzy Broad Offline Reinforcement Learning | 10.1109/TMRB.2026.3699233 | 在体呼吸补偿条件、状态输入和控制接口 |

第二批为按具体论证需要再下载，不要求为了本轮报告全部凑齐：L08（10.1016/j.ymssp.2021.107819）、L09（10.1016/j.mechatronics.2016.08.001）、L31（10.1109/IROS47612.2022.9981146）；理论若采用PPC，连读L20/L33/L34三篇（10.1109/TMECH.2019.2961478；10.1109/TMECH.2020.3019651；10.1109/TMECH.2024.3382261）；恢复切换参考L32（10.1007/s12555-023-0426-3）。图像运动任务比较可补L15（10.1109/JSEN.2025.3598958）。L29/L30是非优先预印本，可从原站下载；不应优先于首批机制近邻。

## C. F10 原文核读结论

原文：Pedro Moreira、Nabil Zemiti、Chao Liu、Philippe Poignet，*Viscoelastic model based force control for soft tissue interaction and its application in physiological motion compensation*，CMPB116(2):52–67，2014，DOI10.1016/j.cmpb.2014.01.017。本地路径：[F10.pdf](E:/github_project/scan-impedance-update/literature/papers/field_supplements/F10.pdf)。以下PDF页码自1计，刊页52–67。独立paper_reader读原文并渲染核对关键公式，主代理与本地文本交叉核查。

### C1. 它已经做到什么

1. **黏弹建模。** 比较多个组织模型，采用Kelvin–Boltzmann模型，以压入松弛试验的机器人运动和测力离线最小二乘辨识参数。组织模型中的x明确是变形量（deformation）。见§3.1–3.3，PDF4–5/刊55–56；Table1。
2. **状态反馈＋Kalman AOB。** 将机器人动力学补偿后与接触模型组成状态空间，用状态反馈实现期望力响应；AOB附加随机游走active state，以补偿状态估计误差以及实验中体现的运动/模型扰动。见§4.1–4.2，PDF7–10/刊58–61，式(12)–(22)。
3. **在线刚度不是遗漏项。** §5.1明确β由Kalman滤波在线估计，α、γ保持离线辨识值。在线β更新AOB矩阵和反馈增益。作者解释，平稳交互可能不足以激励全部黏弹参数，且β失配更影响稳定性，所以只在线更新β。见PDF12/刊63。**不能把F10归为固定刚度的普通扰动补偿。**
4. **未知运动下调力。** D2M2执行接触力控；Viper S650搭载PVC仿体产生呼吸/心跳扰动，D2M2无其运动先验。呼吸为0.2Hz、幅值5mm正弦；心跳为已记录在体运动的回放。见§5.3，PDF12–14/刊63–65，Fig.12–16。这里“运动未知”是控制器不接收其先验，不代表实际人体临床验证。

### C2. 刚度、运动从哪里来，依赖什么

| 信息/接口 | 原文实际给出的来源 | 应保留的边界 |
|---|---|---|
| 接触力 | 力测量，作为状态空间输出y；AOB据此校正 | 力噪声限制AOB调得多快 |
| 组织模型参数 | α、γ离线压入试验；β在线KF | 未完整给出β-KF的状态/量测方程与运动下可辨识性 |
| 环境运动 | 第二机器人产生；控制器借测力与active state抑制，不接收运动先验 | active state不是已证实的物理表面位置/速度；不能与显式运动量测画等号 |
| 控制执行 | 动力学逆与关节广义力/转矩；式(4)τc=JᵀFc、式(6)笛卡尔动力学补偿 | 不是只给商业位置控制器发送位置增量 |
| “1kHz位置控制” | §5.3说的是产生扰动的Viper | 不能移植为D2M2的力控制接口结论 |
| 接触模型中的x | §3定义为组织变形/相对压入；§4机器人模型又用单轴机器人坐标x | 移动环境下如何取得相对压入、提供给β-KF，原文交代不足；不要替作者假设外部位移传感器，也不要假设已用编码器完美重建 |

**AOB估计的对象。** 原文用p_k=L e_k表示状态估计误差经反馈增益后的作用，p按随机游走演化，并由Kalman增益校正（式(17)–(22)，PDF9/刊60）。它处理的是反馈中的误差/扰动作用，不是“k、表面位置、速度”联合物理状态的公开推导。β-KF与AOB并存，也不自动等于刚度和运动被无混淆分离。

**原文记号冲突。** Table1、Fig.5/式(11)对应 f+γ ḟ=βx+αẋ；§5.1却印为 f=−α ḟ+γẋ+β_kx。独立reader直接渲染核对，确认是原文中的α/γ位置不一致，而非仅文本抽取错误。可确定β在线，另外两参数离线；复现前应解决这一符号冲突。本报告不据此否定其实验，也不把它作为本文创新依据。

### C3. 量化证据到哪里

| 实验 | 原文结果 | 能支持/不能支持 |
|---|---|---|
| 静态离体组织，§5.2，Fig.10–11 | 参数α=34.3、β=200、γ=0.024；期望极点−10、设计上升时间213ms；作者报告无超调/静差，给出β估计和协方差曲线 | 支持该实例的力跟踪及在线β变化；没有独立刚度真值/误差指标 |
| 运动PVC仿体，§5.3–5.4，Fig.13–16 | 列参数α=5、β=290、γ=0.02；设计上升时间80ms；目标力2N；作者报告呼吸补偿率87%、心跳79% | 支持未知回放运动下抑制力波动；不是物理运动恢复精度，不是跨数十倍刚度扫描验证 |
| 呼吸/心跳误差描述 | 呼吸“无控峰值误差约0.9N、控制后峰峰约0.1N”；心跳约1.4N降至0.3N | 原文的峰值/峰峰口径并不完全统一，不用近似值重算更精确补偿率 |

运动实验没有对应β估计曲线、真值或“运动泄漏到刚度估计”的对照。Fig.11属于静态实验；不能拿它证明运动中刚度准确。另一方面，闭环能抑制运动是实际已有能力，不能因缺显式分离而把它描述成不能处理运动。

### C4. 与本文的具体比较

| 比较维度 | F10 | 本文当前方法及需证明的增量 |
|---|---|---|
| 刚度认识 | 在线β-KF，其他黏弹参数离线；未给运动分离条件 | 投影消去指定低阶漂移子空间；需给非零剩余激励及模型外残差影响 |
| 运动处理 | active state合并误差/扰动，无显式表面运动状态 | 以df/k−dy形成等效运动增量，包含刚度变化/黏弹项；不得称物理运动恢复 |
| 量测与反馈用刚度 | β更新观察器/反馈，未见本文这种双用途刚度划分 | 量测取未加裕量k̂、反馈取保守k̄；需量化二者误差、相关噪声与闭环后果 |
| 预测 | AOB随机游走；本篇未展示IMM/按位置伺服时延的显式运动预测 | IMM跟踪与时延前馈；IMM/预测本身不新，重点是刚度依赖量测及延迟耦合 |
| 接口 | 模型补偿/转矩控制 | 商业位置接口；需包含内环滞后/延迟验证，接口限制提供任务价值而非自动创新 |
| 稳定性 | 特定模型及失配的稳定/裕度分析 | 目前冻结刚度分析，变增益/恢复切换仍未证明；不能借F10或H04替代 |

**F10一句话定位：** 它是“在线更新组织刚度并用力反馈补偿未知生理运动”的直接先例；本文要进一步回答的是有限窗口、共用探头测力和位置接口下的辨识条件与估计—预测—反馈耦合，而不是重新提出刚度适应或未知运动补偿。

## D. 初步查新判断、检索边界与复核记录

### D1. 问题与方法分开判断

**问题新颖性。** 运动与变刚度导致接触力控制困难，是既有问题；“在线刚度＋生理运动补偿”已有F10，“无视觉力预测补偿”已有CEP2015。当前问题值得研究，但宽泛问题框架的新颖性低。残余目标应保留完整任务约束，同时说明这些约束引入的真实辨识/闭环困难，不能仅缩到一个接口就宣称新问题。

**方法新颖性。** 投影剔除指定漂移、刚度标定的等效运动观测与两种刚度分工有待对近邻全文核查；IMM、KF、主动激励、按柔度缩放和预测前馈都是可继承工具。df/k−dy在冻结局部弹性模型下首先是接触关系的代数结果，科学增量应落在其量测误差/控制耦合及可验证能力。成熟模块的组合不自动失败，也不自动新颖。**当前结论：UNCERTAIN／PROCEED WITH CAUTION；不是NOVEL认证。**

建议将草稿的宽泛缺口改为：

> 已有研究能够通过接触力反馈补偿未知生理运动，并在线更新接触参数。本文关注共用探头力觉和位置接口下，有限窗口局部刚度辨识与等效运动补偿的条件及误差耦合；重点检验漂移投影、量测与反馈刚度分工能否在跨刚度和时延条件下带来额外闭环收益。

这段是本轮比较推断，不是某篇文献的原话。尚需首批全文才能进一步界定相对已有方法的增量。

### D2. 实际检索与未完成范围

- 24条主/细化OpenAlex查询，保存360条返回记录、去重321条身份；包括五类任务、Cortesão/Moreira名字扩展、联合位置估计、HC模型、时变/切换稳定性。**这些是检索返回数，不是321篇已核读或全部摘要纳入判定。** 本报告筛选35条候选，保留原始查询JSON以便追查。
- 另外按精确题名/DOI补查Moreira2015、Beber2024、切换稳定性、PPC反例/反驳和2026秋季预印本。OpenAlex对反例标题检索漏回原条目，Crossref精确题名返回DOI后才由OpenAlex按DOI找到，说明不能把一次零/低相关结果当作没有该文献。
- S2初始主题搜索返回HTTP429；后续按TMECH2017的DOI重试成功。S2该条摘要由出版方屏蔽，身份可核，机制来自OpenAlex实际摘要。未宣称完成S2系统覆盖。
- ARIS `verify_papers.py` 经Crossref核实26条第一组、8条补充组的正式出版身份（含1篇未列入A表的MedIA综述）；两条秋季arXiv由原站和OpenAlex核实，但helper返回`verify_pending`，保留这一限制。
- 本地存在性按manifest、manual_inbox及PDF前两页题名/DOI逐项核对。N04在现状中有已读arXiv版本记录，但项目当前PDF目录未定位，明确标记，未误用同作者T-RO2023 PDF充数。未扩大搜索到用户整台磁盘。
- 非arXiv新候选没有自动抓取全文，集中列入用户下载清单；Beber IROS2024的arXiv直接下载长时间未完成，已中止，仍按摘要级处理。未把网页片段当作核读全文。
- T-RO/IJRR/TMECH/RA-L/ICRA/IROS/TBME/MedIA/Mechatronics/CEP/TASE为优先，不排除AIM、MSSP、TMRB等可能改变最近邻判断的来源。尚未做每个期刊逐卷逐期查全、完整前向引文遍历或系统综述式全部摘要筛选，**不宣布领域饱和**。最关键缺口已具体化为首批7篇原文和1篇勘误。

### D3. 独立复核

按novelty-check运行了一个无历史上下文的新子代理，分别审查P1/P2与M1–M4。技能中历史指定的`gpt-5.6-sol`在本宿主不可用，实际采用继承模型、xhigh，结果记录为**provisional same-family**；没有冒称跨模型认证。独立意见为：宽泛问题框架重合风险高；投影不变性/双刚度设计增量不确定；IMM与归一化的一般机制已有明确先例；应继续保留Rossi与Wijayarathne作为首要基准。

该复核完成后又补检出L35及PPC争论链，因此这几项的具体判断属于主代理新增、未冒称经过该次独立复核。review prompt/response及元数据保存在本scratchpad的`.aris/traces/novelty-check/2026-10-31_run01/`，项目内无写入。
