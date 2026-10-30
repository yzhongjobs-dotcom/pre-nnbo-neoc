# Active Field Map

coverage_status: SUFFICIENT
coverage_record:
research_effort_budget: {"queries": 32, "accepted_evidence_cards": 27}
stopping_reason: 已核读当前三篇人工补充的决定性直接近邻并修正路线边界；定向精确DOI复查无新增候选。结合此前机制骨干、近邻及引用链轮次，申请独立判断本范围的工作饱和，不声称无遗漏。

## Field Core Purposes

```json
[
  "在运动且接触性质变化的组织扫描中维持有用的接触力、运动和柔顺性；理解环境信息如何进入参考、增益与执行。调研本身不选择新方法。"
]
```

## Typical Tasks And Scenarios

```json
[
  "恒力超声扫描为任务核心；未知环境制造/交互中的直接控制机制和触诊辨识为近邻；纯重建/规划及一般优化留作背景。"
]
```

## Core Bottlenecks

```json
[
  {
    "id": "B1",
    "text": "机器人/环境不确定性和观测条件影响力调节及模型适应。环境辨识与机器人阻抗调节必须分开。"
  },
  {
    "id": "B2",
    "text": "运动精度、柔顺性和接触响应之间的任务权衡。"
  },
  {
    "id": "B3",
    "text": "参数时变、交互模型及执行接口决定稳定性保证的适用条件。"
  },
  {
    "id": "B4",
    "text": "学习样本、计算/训练负担、表达限制和泛化之间的权衡；不据综述认定本项目的实际计算瓶颈。"
  }
]
```

## Method Families

```json
[
  {
    "id": "M1",
    "mechanism": "误差自适应、结构参数前馈与鲁棒反馈直接调节力/运动。 DOSMAC将刚度、位置误差和输入增益变化合并为扰动后补偿，原实物共同变化来自固定曲面扫描。 R11跨周期适应前馈/K/D；DO-MPAC以ESO补偿把力误差近似转成双积分器，再解析选择虚拟动作，m/b固定、k=0。 显式操作空间力/力矩调节通过力矩执行绕过位置导纳环节，硬件条件需分开。",
    "source_ids": [
      "local-E32",
      "local-G22",
      "jh4Hbf-KuS8J",
      "local-S01",
      "mMyBOC13ISAJ",
      "local-S03",
      "BIP_InW4FSIJ"
    ]
  },
  {
    "id": "M2",
    "mechanism": "在线环境/预测参数辨识进入增益规则、参考QP或预测运动决策；区分估计变量、增益变量、动作和目标。 MHE-VIC联合估计环境阻抗和机器人扰动，按Kr=估计Ke及阻尼规则更新，法向表面预先测量并视为静态。 E41直接把增益作为优化变量，OS-QP与MS非线性求解的计算和约束条件不同。 E45明确离线OPFB-ADP训练后执行历史反馈策略；环境力学未知常量、位置子系统变化，不能称在线材料参数更新。 C11为相邻的预标定组织模型/几何驱动目标力分支，点间由图像继续更新目标，不归作在线增益选择。",
    "source_ids": [
      "GkLovpurDJwJ",
      "local-H01",
      "local-G02",
      "Ae8L1SwFxN4J",
      "local-S02",
      "local-E41",
      "local-S03",
      "tS1rA372YsUJ",
      "CkE-5mElf_IJ"
    ]
  },
  {
    "id": "M3",
    "mechanism": "示教或策略学习形成接触行为；超声图像收益可以驱动目标力增量，不自动意味着阻抗增益学习。 E27学习增益与参考增量，G23联合输出位置与增益；ACT-US则学习名义轨迹并用图像评分选择目标力，呼吸补偿由导纳执行。 跨周期适应、离线值迭代与预训练策略是不同学习方式。",
    "source_ids": [
      "local-G22",
      "rYX23DX-cGMJ",
      "local-E27",
      "local-G23",
      "zws_WiLclY0J",
      "mMyBOC13ISAJ",
      "tS1rA372YsUJ"
    ]
  },
  {
    "id": "M4",
    "mechanism": "运动预测补偿参考，或主动/被动机械柔顺结构隔离扰动。 ACT-US摘要视觉前馈与方法段落力信号呼吸预测的连接尚不充分，不据此补造接口。 TAFAC补偿已知给定参考的速度/加速度，不等于估计未知组织运动。",
    "source_ids": [
      "local-G13",
      "wftPH5bwT0gJ",
      "zws_WiLclY0J",
      "local-S04"
    ]
  },
  {
    "id": "M5",
    "mechanism": "几何/材料联合估计及稳定执行为相关控制基础，二者不等同为同一控制器。",
    "source_ids": [
      "local-G21",
      "531Glpcm7LgJ",
      "QC0NU0Ejk2UJ",
      "Ns0tgl-gOSMJ"
    ]
  }
]
```

## Family Development Traces

```json
[
  {
    "transition_id": "T-REFERENCE-CONSISTENCY",
    "family": "M2",
    "previous_problem_or_bottleneck": "Lin2021作者指出，未知环境中的给定运动/力参考不一定同时可行，且切换时存在运动失配。这里记录作者的设计动机，不把其文献概括当作全部先前方法的独立复核结论。",
    "progress_and_conditions": "其环境估计驱动的参考QP在稳定单接触、解耦KV与有界机械臂误差条件下重规划参考；Panda转移和斜面实验支持该实现的接触/跟踪能力。",
    "residual_or_new_bottleneck": "RLSE激励要求和首次接触位置处理仍有限；目标矩阵按预设规则变化，尚不是任务收益驱动的在线最优增益选择。",
    "research_question_shift": "作者将设计重心从直接给定目标误差模型推进到未知环境约束下可执行参考与内层模型的协调。不是声称整个领域按同一时间线演进。",
    "subsequent_direction": "作者在结论提出用强化学习自动选择目标模型，属于未来工作陈述；后续是否实现须由前向引用链核查。",
    "transition_problem_status": "partially_addressed",
    "evidence_ids": [
      "GkLovpurDJwJ"
    ]
  },
  {
    "transition_id": "T-LUMPED-COMPENSATION",
    "family": "M1",
    "previous_problem_or_bottleneck": "未知刚度与接触位置变化时，单纯给定导纳及误差积累存在力偏差。",
    "progress_and_conditions": "DOSMAC把环境变化写入总扰动，用ESO与滑模补偿，在单解耦低频模型及软硬/固定曲面扫描中改善原实验RMSE。",
    "residual_or_new_bottleneck": "总扰动处理支持调力但不提供独立材料参数用于增益选择；实物独立呼吸及真实位置误差作用尚未由所核实验建立。",
    "research_question_shift": "从逐项恢复物理量扩展到补偿其对控制的合并影响；这里只描述该论文形成的替代路线。",
    "subsequent_direction": "精确材料估计与总扰动调力各有用途，当前证据不要求领域统一转向其中之一。",
    "transition_problem_status": "partially_addressed",
    "evidence_ids": [
      "local-S01"
    ]
  },
  {
    "transition_id": "T-FEEDFORWARD-IMPEDANCE",
    "family": "M3",
    "previous_problem_or_bottleneck": "重复交互中只学习前馈无法稳定位置/速度发散作用。",
    "progress_and_conditions": "R11原文Fig5比较支持前馈与K/D共适应的分工，周期参考、机器人模型及有界集合收敛条件明确；后来的E27和G23分别使用预训练/历史修正和状态策略输出增益与参考。",
    "residual_or_new_bottleneck": "周期适应、训练分布与动作接口各有限制；学出增益不等于识别材料，也不自动证明单次动态组织任务收益。",
    "research_question_shift": "从只抵消重复负载转向同时塑造反馈响应，再形成不同训练与联合动作表示；这是证据支持的机制关系，不断言作者之间未核准的直接发现顺序。",
    "subsequent_direction": "以在线模型/历史或策略联合选择参考与增益，并单独报告训练成本和任务收益。",
    "transition_problem_status": "partially_addressed",
    "evidence_ids": [
      "mMyBOC13ISAJ",
      "local-E27",
      "local-G23"
    ]
  },
  {
    "transition_id": "T-LUMPED-PREDICTIVE",
    "family": "M2",
    "previous_problem_or_bottleneck": "DO-MPAC原文直接引用DOSMAC：总扰动补偿后仍需协调瞬态速度、超调与抖振。",
    "progress_and_conditions": "DO-MPAC在相同法向固定导纳m/b,k=0基础上，用ESO及反馈线性化近似得到双积分器，离线解析预测反馈减少原测试误差。原文对照支持本实例推进，不声称所有SMC都会失效。",
    "residual_or_new_bottleneck": "噪声/ESO带宽、突变滞后和稀有超时仍存在；省去显式MPC约束换取解析速度，独立患者运动有待实物检验。",
    "research_question_shift": "以补偿后的力误差动态选择预测虚拟动作，改变瞬态决策；这条替代路线没有要求逐项恢复材料/运动。",
    "subsequent_direction": "该文提出在线约束、多轴与动态点云更新；是作者后续方向而非已验证能力。",
    "transition_problem_status": "partially_addressed",
    "evidence_ids": [
      "local-S01",
      "local-S03"
    ]
  }
]
```

## Problem Method Matrix

```json
[
  {
    "problem": "B1",
    "method": "M1",
    "mechanism": "误差自适应、结构参数前馈与鲁棒反馈直接调节力/运动。 DOSMAC将刚度、位置误差和输入增益变化合并为扰动后补偿，原实物共同变化来自固定曲面扫描。 R11跨周期适应前馈/K/D；DO-MPAC以ESO补偿把力误差近似转成双积分器，再解析选择虚拟动作，m/b固定、k=0。",
    "source_ids": [
      "local-E32",
      "local-G22",
      "jh4Hbf-KuS8J",
      "local-S01",
      "mMyBOC13ISAJ",
      "local-S03"
    ]
  },
  {
    "problem": "B1",
    "method": "M2",
    "mechanism": "在线环境/预测参数辨识进入增益规则、参考QP或预测运动决策；区分估计变量、增益变量、动作和目标。 MHE-VIC联合估计环境阻抗和机器人扰动，按Kr=估计Ke及阻尼规则更新，法向表面预先测量并视为静态。 E41直接把增益作为优化变量，OS-QP与MS非线性求解的计算和约束条件不同。 E45明确离线OPFB-ADP训练后执行历史反馈策略；环境力学未知常量、位置子系统变化，不能称在线材料参数更新。",
    "source_ids": [
      "GkLovpurDJwJ",
      "local-H01",
      "local-G02",
      "Ae8L1SwFxN4J",
      "local-S02",
      "local-E41",
      "local-S03",
      "tS1rA372YsUJ"
    ]
  },
  {
    "problem": "B2",
    "method": "M3",
    "mechanism": "示教或策略学习形成接触行为；超声图像收益可以驱动目标力增量，不自动意味着阻抗增益学习。 E27学习增益与参考增量，G23联合输出位置与增益；ACT-US则学习名义轨迹并用图像评分选择目标力，呼吸补偿由导纳执行。 跨周期适应、离线值迭代与预训练策略是不同学习方式。",
    "source_ids": [
      "local-G22",
      "rYX23DX-cGMJ",
      "local-E27",
      "local-G23",
      "zws_WiLclY0J",
      "mMyBOC13ISAJ",
      "tS1rA372YsUJ"
    ]
  },
  {
    "problem": "B1",
    "method": "M4",
    "mechanism": "运动预测补偿参考，或主动/被动机械柔顺结构隔离扰动。 ACT-US摘要视觉前馈与方法段落力信号呼吸预测的连接尚不充分，不据此补造接口。",
    "source_ids": [
      "local-G13",
      "wftPH5bwT0gJ",
      "zws_WiLclY0J"
    ]
  },
  {
    "problem": "B3",
    "method": "M5",
    "mechanism": "几何/材料联合估计及稳定执行为相关控制基础，二者不等同为同一控制器。",
    "source_ids": [
      "local-G21",
      "531Glpcm7LgJ",
      "QC0NU0Ejk2UJ",
      "Ns0tgl-gOSMJ"
    ]
  }
]
```

## Assumption Effectiveness Failure Matrix

```json
[
  {
    "family": "M1",
    "assumptions": [
      "有界结构参数和不确定性；内外环接口",
      "R11周期历史与机器人模型条件；DO-MPAC可微有界扰动变化率、名义Schur稳定。"
    ],
    "effective": [
      "ARIC气球变力接触实验",
      "DO-MPAC移动表面/时变刚度仿真与静态异质曲面五次重复实物。"
    ],
    "failure": [
      "模型K无特定物理意义，不等于真实刚度识别",
      "ESO高带宽放大噪声、突变滞后；硬件没有独立呼吸平台，执行饱和不等于约束MPC。"
    ],
    "source_ids": [
      "local-E32",
      "local-G22",
      "jh4Hbf-KuS8J",
      "local-S01",
      "mMyBOC13ISAJ",
      "local-S03"
    ]
  },
  {
    "family": "M2",
    "assumptions": [
      "可用接触模型、激励和估计/控制动态条件",
      "E45环境参数未知常量、线性位置子系统、端点参考、PE采集、完美内环；离线训练后执行。"
    ],
    "effective": [
      "原文参考/增益和预测控制实例"
    ],
    "failure": [
      "Roveda带宽分离、Lin阻尼激励、Wang力方向平衡速度近似各自限定适用性",
      "DO-MPAC通过离线解析增益节省在线计算，但没有在线增益优化或显式MPC约束。"
    ],
    "source_ids": [
      "GkLovpurDJwJ",
      "local-H01",
      "local-G02",
      "Ae8L1SwFxN4J",
      "local-S02",
      "local-E41",
      "local-S03",
      "tS1rA372YsUJ"
    ]
  },
  {
    "family": "M3",
    "assumptions": [
      "观测、训练分布、奖励代理和约束定义",
      "R11周期历史；E27随机预训练及前试次回归；G23状态策略；ACT-US多模态示教和力区间弱标签。"
    ],
    "effective": [
      "TASE脊柱模体目标力策略",
      "原文支持增益/参考学习，或图像目标力与呼吸补偿，两者不同。"
    ],
    "failure": [
      "期望累计约束不等于逐时刻硬保证",
      "ACT-US表格RMS、标准差与结论方差措辞不一致；不照搬综合百分比。"
    ],
    "source_ids": [
      "local-G22",
      "rYX23DX-cGMJ",
      "local-E27",
      "local-G23",
      "zws_WiLclY0J",
      "mMyBOC13ISAJ",
      "tS1rA372YsUJ"
    ]
  },
  {
    "family": "M4",
    "assumptions": [
      "视觉周期运动模型，或额外机械恒力机构"
    ],
    "effective": [
      "Cao呼吸平台和Wu软组织扫描"
    ],
    "failure": [
      "重拟合瞬态、延迟与残余机械振荡"
    ],
    "source_ids": [
      "local-G13",
      "wftPH5bwT0gJ",
      "zws_WiLclY0J"
    ]
  },
  {
    "family": "M5",
    "assumptions": [
      "先验几何/准静态模型，或能量/功率模型条件"
    ],
    "effective": [
      "CARE、Chalasani估计和装配/切割执行实例"
    ],
    "failure": [
      "估计不直接验证动态扫描；切割/装配不是超声控制的直接实证"
    ],
    "source_ids": [
      "local-G21",
      "531Glpcm7LgJ",
      "QC0NU0Ejk2UJ",
      "Ns0tgl-gOSMJ"
    ]
  },
  {
    "family": "M1",
    "assumptions": "C10：力矩执行、非奇异构型、动力学补偿及测量准确；柔顺接触允许较大增益。",
    "effective_conditions": "C10移动模体与异质曲面扫描支持力/力矩调节；1.6Hz目标力的EFC/调参IAAC滞后5.76/121度。",
    "failure_conditions": "不等同于位置接口可直接迁移；未知更硬接触性能、在线阻抗选择和临床诊断效用均未建立。",
    "source_ids": [
      "BIP_InW4FSIJ"
    ]
  },
  {
    "family": "M2",
    "assumptions": "C11：预标定力—应变、点云几何、相邻点归一化变形连续与图像置信度。",
    "effective_conditions": "三种乳房仿体、125Hz机器人/1kHz力采样、1mm/s扫描中调整目标力并作固定映射力反馈。",
    "failure_conditions": "独立组织运动和在线阻抗更新未验证；式11–12方向疑点不改写为推定代码错误。",
    "source_ids": [
      "CkE-5mElf_IJ"
    ]
  },
  {
    "family": "M4",
    "assumptions": "TAFAC：已知xin速度/加速度、线性KV与小控制周期，m/b固定。",
    "effective_conditions": "单一橡胶15N，恒速/恒加速度参考，500Hz；减小超调和误差，但20mm/s RMSE略高于AVAC。",
    "failure_conditions": "原文承认复杂人类轨迹与变刚度未充分验证；时延证明用一阶近似，不是任意通信延迟保证或未知患者运动恢复。",
    "source_ids": [
      "local-S04"
    ]
  }
]
```

## Consensus

```json
[
  "环境参数、机器人阻抗、参考和控制动作是不同更新对象。",
  "参考/阻抗协调、动态力补偿和联合接触估计均有具体先例，按条件比较。",
  "成像收益与接触误差有关但非同一指标；期望约束、被动性与逐时刻安全不能混用。",
  "显式力矩执行、位置伺服导纳、已知参考前馈和图像目标力调节是不同控制/决策层；运动任务或信息进入控制已有多种有效路线，不能以是否采用刚度/漂移分解划定现状。"
]
```

## Unresolved Contradictions

```json
[
  "Lin摘要对在线位置估计的概括需与首次接触备注区分。",
  "Cao名义视觉更新时间与人为延迟试验阈值不可直接等同。",
  "DOSMAC正文称已处理刚度与位置共同变化；所核实物实验是固定曲面沿扫描高度变化，不等于独立生理运动验证。二者保留其层次差异。",
  "ACT-US摘要强调视觉前馈，而方法主要写呼吸力带通/正弦预测，视觉和速度前馈的接口尚未完整给出；不得补造。",
  "DO-MPAC概括措辞称处理动态约束，但式25/Remark4为无约束解析预测控制，硬饱和在执行端。",
  "C11发表式11–12在正力、K=6时只保持或增加目标力，未体现低置信度纠正增力；原文图像反馈解释与公式方向疑点保留，未据此推断实验代码或整篇无效。"
]
```

## Coverage Record

```json
{
  "coverage_status": "SUFFICIENT",
  "research_effort_budget": {
    "queries": 32,
    "accepted_evidence_cards": 27
  },
  "stopping_reason": "已核读当前三篇人工补充的决定性直接近邻并修正路线边界；定向精确DOI复查无新增候选。结合此前机制骨干、近邻及引用链轮次，申请独立判断本范围的工作饱和，不声称无遗漏。",
  "coverage_gaps": [],
  "blind_spots": [
    "27卡为所注明选段，不是全文独立复证或数学定理认证。",
    "Jiang2024独立身份未闭合；C10已提供IAAC定义与调参比较，S03有复现比较；不据该待核引文作未检查能力或唯一性判断。",
    "少量旧身份/版本/公开日期仍未核准，仅保留发现背景，不承担决定性前沿主张。",
    "界定到既定法向接触控制问题，不覆盖纯成像/规划或整个机器人医学文献。"
  ]
}
```

## Unresolved Problem Leads

```json
[
  "本轮只补查上述能力和条件，不从综述未来工作直接推出新idea。"
]
```
