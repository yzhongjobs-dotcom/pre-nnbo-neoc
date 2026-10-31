# 2022–2026 直接相关文献补查（ARIS，balanced）

日期：2026-10-31。范围截至本次检索可见记录；原为会话临时输出，由Claude Code归档入项目。按指定 research-lit 的机制分析和证据边界执行，采用用户指定的 OpenAlex/Semantic Scholar 子技能及 ARIS helper；Crossref仅用于精确身份核对。本次是定向补查，未启动写入项目的 Controller 正式流程，也无独立覆盖审查。

**结论：有新的近期相关工作，不能说“近几年没有”；但本轮没有核实到同时满足“位置接口、仅力与编码器、运动中在线组织刚度、力换算运动并预测前馈、量测增益误差的反馈耦合分析”的新论文。后半句是本轮检索/阅读结果，不是领域不存在的证明。** 新候选中最需警惕的是R03的刚度—力联合估计与R04的接触参数可观性/MPC；最直接的超声任务补充是R01、R05。新检出的T-RO 2026系统实际是已知C03，未重复计新。

## A. 新候选表

类别①超声运动/接触控制（含明确标出的几何旁支）；②接触参数估计并服务控制；③从力估运动并补偿。**本轮没有新增已核实属于③且直接满足目标任务的论文**。摘要没有交代的刚度、运动或接口能力统一记“待核”，不推断其不存在。下表均有正式出版DOI；作者稿只是阅读版本，不重复计作预印本新论文。

|编号|题目|作者|年份|venue|DOI|类别|刚度/运动信息来源|接口|与本文关系|证据级别|来源链接|
|---|---|---|---|---|---|---|---|---|---|---|---|
|R01|Precise Repositioning of Robotic Ultrasound: Improving Registration-Based Motion Compensation Using Ultrasound Confidence Optimization|Zhongliang Jiang、Nehil Danis、Yuan Bi 等|2022|TIM|[10.1109/TIM.2022.3200360](https://doi.org/10.1109/TIM.2022.3200360)|①体动/重定位|运动：RGB-D分割、ICP与US图像微调；控制器刚度为设定值，非在线组织刚度|力矩/笛卡尔柔顺，式(3)|必须正面比较的近邻（体动层）|开放作者稿全文已读，v2；正式版未逐字比对|[OpenAlex](https://openalex.org/W4292793935)；[S2](https://www.semanticscholar.org/paper/7c09e52580cd354fefd736161372b5f2d1c17e1d)；[开放正文](https://arxiv.org/html/2208.05383v2)|
|R02|Towards Autonomous Atlas-Based Ultrasound Acquisitions in Presence of Articulated Motion|Zhongliang Jiang、Yuan Gao、Le Xie 等|2022|RA-L|[10.1109/LRA.2022.3180440](https://doi.org/10.1109/LRA.2022.3180440)|①关节姿态/路径适配|MRI图谱＋RGB-D/OpenPose＋非刚性配准；US分割居中；机器人刚度设定|关节力矩传感的柔顺执行；底层API未核|应引用的背景；非持续呼吸力控|开放作者稿§IV–VI已读（局部），v1|[OpenAlex](https://openalex.org/W4285194566)；[开放正文](https://arxiv.org/html/2208.05399v1)|
|R03|Composite Disturbance Filtering for Interaction Force Estimation With Online Environmental Stiffness Exploration|Yanran Wei、Jiayin Wang、Wenshuo Li 等|2025刊期；2024在线|TMECH|[10.1109/TMECH.2024.3443310](https://doi.org/10.1109/TMECH.2024.3443310)|②刚度/力联合估计|机器人动力学＋力生成模型；EM探索刚度并估力；组织运动如何进入模型待核|摘要未说明执行接口；有动力学不等于开放力矩指令|必须正面比较的近邻（辨识层，待全文）|仅摘要（OpenAlex实际摘要）|[OpenAlex](https://openalex.org/W4402041477)|
|R04|Differentiable Compliant Contact Primitives for Estimation and Model Predictive Control|Kevin Haninger、Kangwagye Samuel、Filippo Rozzi 等|2024|ICRA|[10.1109/ICRA57147.2024.10611406](https://doi.org/10.1109/ICRA57147.2024.10611406)|②参数/几何估计→MPC|关节位置＋电流/力矩与动力学；EKF估刚度、静止点或接触位置；MPC演示在线更新静止点、刚度固定|MPC输出阻抗参考位置；依赖底层笛卡尔阻抗/力矩模型，非裸位置接口|必须正面比较的近邻（联合辨识/控制层）|开放作者稿全文已读，arXiv v3；正式版未逐字比对|[OpenAlex](https://openalex.org/W4401417026)；[开放正文](https://arxiv.org/html/2303.17476v3)|
|R05|AI-Powered Robust Interaction Force Control of a Cardiac Ultrasound Robotic System|Ehsan Zakeri、Amanda Spilkin、Hanae Elmekki 等|2025刊期；2024在线|TIE|[10.1109/TIE.2024.3451138](https://doi.org/10.1109/TIE.2024.3451138)|①学习图像反馈＋鲁棒力控|US图像CNN确定目标力；力反馈＋自适应滤波准滑模；是否估组织刚度/运动未由摘要确定|混合位置/力控制；实际位置或力矩接口待核|必须正面比较的近邻（超声力控层）|仅摘要（S2/OpenAlex）|[OpenAlex](https://openalex.org/W4402569920)；[S2](https://www.semanticscholar.org/paper/735c71ca6b4ab7517d2e8b6a0ed73b3ccc825cc4)|
|R06|Ultrasound-Guided Robotic Autonomous Operation Based on Real-Time Deformation Tracking and Prediction|Teru Chen、Xingwei Zhao、Yuhao Zhang 等|2025刊期；2024在线|TII|[10.1109/TII.2024.3477549](https://doi.org/10.1109/TII.2024.3477549)|①图像运动预测（US引导手术）|US非刚性配准追踪组织/路径变形；遗忘因子RLS更新预测模型；非从接触力换算运动|动态向量场引导机器人；底层执行接口待核|必须正面比较的近邻（预测层；非表面恒力扫描）|仅摘要（S2/OpenAlex）|[OpenAlex](https://openalex.org/W4404132652)；[S2](https://www.semanticscholar.org/paper/1c59cf6493f38b8fa7d3c94143a5098cd72dabf3)|
|R07|Force-Based Viscosity and Elasticity Measurements for Material Biomechanical Characterization With a Collaborative Robotic Arm|Luca Beber、Edoardo Lamon、Giacomo Moretti 等|2025|TIM|[10.1109/TIM.2025.3581663](https://doi.org/10.1109/TIM.2025.3581663)|②黏弹辨识/触诊|协作臂力学接触测黏度、弹性；以高精度仪器作真值；环境运动处理与具体在线估计器待核|协作臂；实际执行接口待全文核对|应引用的背景；Beber系列后续|仅摘要（OpenAlex）|[OpenAlex](https://openalex.org/W4411484126)|
|R08|Tactile-Guided Robotic Ultrasound: Mapping Preplanned Scan Paths for Intercostal Imaging|Yifan Zhang、Dianye Huang、Nassir Navab 等|2025|IROS|[10.1109/IROS60139.2025.11246299](https://doi.org/10.1109/IROS60139.2025.11246299)|①/②触觉几何/肋间路径|触觉相关一维信号＋机器人跟踪生成肋骨点云，与图像骨表面配准；不是摘要已证实的刚度—运动联合估计|摘要未说明|应引用的背景（TUM后续，非呼吸补偿证据）|仅摘要（OpenAlex/S2）|[OpenAlex](https://openalex.org/W4416750255)|
|R09|Safe Contact Force Generation for Robotic Thyroid Ultrasound Imaging|Ryosuke Tsumura、Toshifumi Tomioka、Yoshihiko Koseki 等|2024|RA-L|[10.1109/LRA.2024.3349959](https://doi.org/10.1109/LRA.2024.3349959)|①超声接触力/机械顺应|串联弹性末端：线性执行器、力传感器、弹簧/恒力弹簧；靠弹簧压缩调力，不以环境刚度换算运动|末端线性执行器＋机械顺应；不需机械臂导纳反馈|应引用的背景（改变执行机构约束）|仅摘要（S2/OpenAlex）|[OpenAlex](https://openalex.org/W4390577974)|

预印本：没有为凑数单列新的直接相关预印本。L29–L30已知，未重列；R01/R02/R04的arXiv版本对应已正式发表论文。

## B. 优先正面比较（5篇新候选＋1篇已知升级）

1. **R04，Haninger等，ICRA 2024：辨识与控制层最重要。** 可微弹簧接触基元将刚度、接触点和静止点参数化，利用机器人动力学、位置与电流/力矩构造EKF，再更新MPC。作者稿§4.3、式(25)已有局部可观性条件；§6.3平面/铰链实验在线估静止点，刚度取固定值，目标力分别为3 N和15 N，MPC预测13步、步长0.03 s。它已证明“估接触位置并改善力控”，不能将联合参数/几何估计或可观性分析笼统写成本文空白。区别应写成本文运动下的投影辨识条件、以在线刚度换算运动量测及其耦合；该稿实验未验证这一完整链条。MPC输出阻抗平衡位置，底层仍使用笛卡尔阻抗及力矩动力学，不能简单标作与裸位置接口相同。[正文§4–6](https://arxiv.org/html/2303.17476v3)

2. **R03，Wei等，TMECH 2025（2024在线）：必须先读全文的辨识近邻。** 摘要明确把交互力视作机器人状态耦合扰动，联合机器人动力学和力生成模型，进行可分离性分析，以EM探索环境刚度并估计交互力，在手术平台验证。它估的是刚度与力；是否也估环境位移、怎样处理移动组织、是否将运动前馈到力环，**摘要不能回答**。因此既不能按“observer”归入固定刚度路线，也不能宣布其已解决本文问题。需逐式比较其可分离性条件与本文运动投影条件、参数误差进入控制环的路径。摘要无可比数值。[OpenAlex摘要](https://openalex.org/W4402041477)

3. **R01，Jiang等，TIM 2022：TUM motion-aware系列的直接体动先例。** RGB-D分割检测体动，停止扫描后ICP配准、重算余下路径，再用US置信图修正姿态与拼接。正文§III-C、§V明确面向可见刚体运动；快小生理运动仅由柔顺控制部分适应，未纳入其运动补偿模块。依赖RGB-D、图像与力矩柔顺执行；设定的机器人柔顺刚度不是在线组织刚度。约440 mm血管模体，随机移动范围约110 mm、旋转最高80°，两种遮挡条件各20次，配准补偿误差6.1±1.6/6.9±1.6 mm。它解决停扫后恢复和几何完整性，区别于本文持续法向恒力与时延预测；不能把其重规划误差与本文力RMSE直接对表。[正文§III-C、IV-D、V](https://arxiv.org/html/2208.05383v2)

4. **R05，Zakeri等，TIE 2025：超声学习/鲁棒控制路线。** 高层US图像CNN生成目标力，低层力反馈经自适应滤波准滑模形成混合位置/力控制。心脏超声模体实验中，作者摘要报告对FQSMC、SMC和PI的RMS/STD/最大值等指标改善21.87%–68.25%。这些百分数跨指标/对照，不是统一力RMSE；“cardiac”也不等于已验证呼吸或心脏表面运动。它把图像质量闭环引入力设定，与本文只从接触信号生成运动量测不同；是否在线估刚度、具体扰动和执行接口待正文核对。[OpenAlex摘要](https://openalex.org/W4402569920)

5. **R06，Chen等，TII 2025：图像预测补偿近邻。** 非刚性US图像配准追踪组织及手术路径变形，遗忘因子RLS更新预测模型，预测下一次变形，动态向量场引导机器人跟随变形路径；双臂离体试验报告动态路径误差小于1 mm。其运动信息来自US图像，任务是手术导航，不能冒称恒力扫描成果；但它表明近年US运动预测不限于固定周期模型。应与本文的力源量测、在线换算增益及延迟前馈比较，底层接口待核。[OpenAlex摘要](https://openalex.org/W4404132652)

6. **已知L15，Lei等，Sensors Journal 2025，应从候选升级为正面比较（不算新发现）。** ReME-Net用实时US帧差及时间/空间注意力估6-DOF肝脏呼吸运动，混合力/位置＋自适应变导纳维持接触，明确不依赖独立跟踪设备。摘要报告平移/旋转/力平均误差3.05 mm、2.52°、0.25 N，对应周期运动平移/旋转幅值19.72 mm、17.32°。它需要US图像，仍不满足本文“力＋编码器”的信息约束；但必须把“没有额外运动传感器”与“没有图像运动观测”区分。刚度估计、周期依赖与接口待全文核对。[DOI](https://doi.org/10.1109/JSEN.2025.3598958)；[S2实际摘要](https://www.semanticscholar.org/paper/e4dbaaaddf0b1ae45093efa1c182500f63beefdf)

**对三条路线及缺口的影响。** 可保留三条作为本文控制问题的组织视角，但不宜当作覆盖全领域的互斥分类：

- 路线一补充“参数与接触几何联合辨识”旁支（R04），明确已有可观性研究；R03的分离性分析也需核对后再定位。
- 路线二容纳图像决定目标力与鲁棒内环（R05），区分“学习期望力”和“学习/估计运动”；刚度在线估计已有F10、Beber系列，不能退回“既有方法都用固定刚度”的表述。
- 路线三若保留名称“测量预测再前馈”，须把R01/R02的重定位/路径适配列为相邻的几何补偿，而非称其已做连续预测前馈。若扩大为“运动观测驱动补偿”，则应细分重规划、图像反馈、预测前馈。R06/L15要求把“预测多依赖周期”限于所引具体方法，不能概括所有新近图像方法。

建议缺口只保留为可检验的窄表述：**在本文位置接口与仅力/编码器约束下，在线组织刚度同时参与力控缩放和受控探头的运动量测；需研究运动下辨识条件，以及换算误差使预测前馈经力量测形成反馈后的容许范围和逐窗口接受准则。** 本轮尚无文献证据消除该特定缺口，但R03与下列CEP待核文献可能改变边界，不能据此认证新颖性或写“首次”。

**其他已知但需升级处理。** C03是[10.1109/TRO.2026.3734315](https://doi.org/10.1109/TRO.2026.3734315)，摘要确有深呼吸人体与0.10 mm振荡结果；应作为机械顺应/专用末端路线比较，不能漏引后又声称新发现。C10、S01/S03/S04、L16、L28/L35、G13仍是已知核心。R07为Beber系列新增后续，但摘要只支持材料黏弹测量验证，不能将其升级成已验证运动补偿。

## C. 需用户获取全文的合并清单

优先级P1可能改变缺口/最近邻判断；P2补齐任务比较。均未自动下载或保存PDF。

|优先级|文献/DOI|要核对什么|
|---|---|---|
|P1|R03：10.1109/TMECH.2024.3443310|力生成模型是否包含环境位移/速度；EM估哪些状态；可分离性与激励条件；是否用于控制及执行接口；参数误差是否进入运动补偿回路。|
|P1|L15（已知）：10.1109/JSEN.2025.3598958|ReME-Net的运动观测、刚度与变导纳更新式；运动模式/频率变化实验；接口、延迟与接触控制指标。|
|P1|待判定新线索：10.1016/j.conengprac.2026.107163|Liu等，Composite admittance control framework for nonsmooth force tracking under unknown environmental dynamics。OpenAlex/S2/Crossref核实身份，但三库均未给实际摘要，出版页未读到正文；不得仅凭题名宣称相关能力。核对是否估刚度/运动、非光滑性的来源、接口与补偿结构。OpenAlex记2026-07-20，Crossref卷期2026-11，在线日期仍待出版社首页确认。[OpenAlex](https://openalex.org/W7169811421)；[S2](https://www.semanticscholar.org/paper/2da14d8039fb038f8ea47f6399ee3320f20e666f)|
|P2|R05：10.1109/TIE.2024.3451138|模型及AFQSMC估计量、人体/模体运动条件、执行接口；图像只调目标力还是也观测运动。|
|P2|R06：10.1109/TII.2024.3477549|图像形变与实际组织运动的对应、预测时域/延迟、是否有接触力约束；避免把导航性能当力控性能。|
|P2|R07：10.1109/TIM.2025.3581663|与ICRA/IROS 2024作者工作新增内容、压入深度/刚度联合估计、运动假设与控制用途。|

R01/R02/R04已有在线开放正文可读，本轮不要求为它们重复下载；如用来定稿引用或实现对照，再核对正式版差异。R08/R09目前为背景，不列首批下载任务。

## D. 实际检索记录与覆盖边界

**计数口径：一次提交给provider/helper的查询算1条逻辑查询；helper内部HTTP重试不另算。返回数不等于读过摘要或纳入数。** 检索原始JSON在本文件同目录，以下数量可从其中重算。

|检索环节|逻辑查询数|成功/失败|保存的返回记录|
|---|---:|---|---:|
|OpenAlex主题搜索，2022–2026，16组每组请求20条|16|0成功，16次429|0|
|OpenAlex题名过滤3组＋Jiang作者过滤1组|4|0成功，4次429|0|
|S2 relevance主题查询（已保存）|17|1成功，16失败429|20|
|S2 bulk首批3组|3|2成功，1失败429|1003|
|S2 bulk机制定向10组|10|8成功（含合法零结果），2失败429|393|
|OpenAlex五种子前向引文|5|4成功，Yuen失败429|37|
|OpenAlex重试Yuen引文|1|失败429|0|
|S2 Yuen前向引文|1|成功，返回全时期90条，筛得2022–2026年5条|5（近期）|

另有首次S2探测“robotic ultrasound motion”成功20条，其完整响应未落JSON，故单独记录，**不混入可复算去重数**。因此实际主题/作者检索为51条（含这1次探测），引文查询7次；不是“58条均完成”。在保存池中，主题返回1416条＋近期引文42条＝1458条，按DOI小写归一并合并规范化同题名，去重1431个发现记录；另外20条探测结果与该池重合量未完整计量。这个数是发现池规模，绝不是1431篇完成摘要筛选。

S2 bulk helper一次返回了1000条而非脚本意图的100条，且宽查询混入大量超声治疗、成像预测与微机器人条目；其余大集合也有明显噪声。此次只对标题/venue/机制线索命中者作选择性摘要筛查，再对A表9篇逐条核身份、读实际摘要（其中3篇还读开放正文的全部或指定章节）；**未逐条筛完1431条，更未完成bulk下一页**。A表9篇全部Crossref精确DOI返回同题名/作者/venue；R03/R05/R06明确区分在线年与刊期年，R03不是误写成另一篇2024论文。

**五种子实际前向结果：** Rossi 8条、Moreira 16条、DOSMAC 4条、Beber 9条（均为OpenAlex范围内全页，meta.count等于返回数）；Yuen由S2返回90条且无next，年筛后5条。这里“全页/无next”只表示该库当时索引的结果，没有保证出版方引文完整。

- Rossi链新增R04；TCST的Hybrid Adaptive Vision-Force Control等为相邻方向，未核成本文动态组织最近邻。
- Moreira链新增R03；TMECH 2022 Precision Force Tracking Control…的摘要是压电耳科装置跟踪S形/正弦目标力，没有据此认定它估生理运动。在线辨识、模型切换等其余线索未全部精读。
- DOSMAC链除已知S03外，有C节的CEP 2026待判定项；倾转旋翼文献明确排除，不能凑成机器人接触先例。
- Beber链新增R07/R08；L35已知去重；物体变形操作、周期交互学习等保留为外圈线索，未冒称运动组织扫描实证。
- Yuen近期5条中，G11、L17、C03已知；另两条为电生物阻抗感测的神经电极补偿（10.1109/JSEN.2023.3307489）和工具局部位置传感的移动物体打印（10.1016/j.jmapro.2023.01.020），不作为“接触力估组织运动”新发现。

**实际查询词及台账：** [OpenAlex 16组](aris_evidence/recent_2026-10-31/recent_openalex.json)、[S2 relevance 17组](aris_evidence/recent_2026-10-31/recent_s2.json)、[S2 bulk 3组](aris_evidence/recent_2026-10-31/recent_s2_bulk.json)、[S2机制定向10组](aris_evidence/recent_2026-10-31/recent_targeted.json)、[四条成功及Yuen失败引文](aris_evidence/recent_2026-10-31/recent_citations.json)、[题名/作者过滤及Yuen重试](aris_evidence/recent_2026-10-31/recent_oa_filters.json)、[Yuen S2原始引文](aris_evidence/recent_2026-10-31/recent_yuen_s2.json)、[精确DOI交叉核对](aris_evidence/recent_2026-10-31/recent_enrich.json)。主题覆盖体动/呼吸、Jiang、学习、在线刚度、ESO/扰动观测、MPC、非静止接触、黏弹/触诊、从力估运动、打磨；全部使用2022–2026过滤。精确身份查询和网页正文访问不计入主题查询数。

另试了6条限制在OpenAlex/S2页面的网页索引查询，实际主要返回S2托管PDF/参考文献片段，未将其计作数据库API覆盖，也未据这些片段建立A表科学结论。开放正文仅在线阅读，无PDF保存。没有在未获得回复时擅自扩展为一般网页检索；使用出版物/arXiv页面只核读数据库已经发现的身份。

**明确未完成：** OA主题/作者查询受限；S2 relevance绝大多数失败；未穷举各venue逐年逐卷、未完成宽bulk分页和全量摘要筛选、未做两库五条引文链的全交叉验证；TUM作者序列尚未穷尽；部分正文/接口未知；CEP待判定项可能改变最近邻判断。RA-L/TMECH/ICRA/IROS/TIE/TII已有新增，T-RO/IJRR/TMRB的本轮直接强命中主要为已知文献；不能把这当作TBME/MedIA/Automatica/TCST等没有相关工作。**不能宣称覆盖饱和，也不能把此次未发现完整同机制论文写成“确实没有”。**


**归档说明（Claude Code，2026-10-31）。** 查询脚本与小型台账存于 `aris_evidence/recent_2026-10-31/`；`recent_s2_bulk.json`（3.6 MB）与 `recent_targeted.json`（1.5 MB）为噪声较大的原始批量返回，未入库，可用同目录脚本重跑。R03 摘要经 OpenAlex 复核：只说联合刚度辨识与交互力估计、给出可分离性分析，未提组织运动或执行接口；R04 身份与“关节力矩/位置、EKF在线估静止点/刚度/接触位置、用于MPC”经 arXiv 摘要复核，可观性条件（§4.3、式(25)）为 Codex 读正文所得，未独立复核。
