# ARIS 第二轮：recent-related-work-20261031

## A. Controller 最终状态与用户来源

- 最终为 **HUMAN_SEARCH_REQUIRED**，子请求 **HUMAN_FULLTEXT_REQUIRED / fulltext_download_batch**，`stop: true`；已停止。
- 本轮17张新卡全部PASS并接受，累计 **12 → 29张**；Field Map通过结构验证，独立覆盖审查为 **CONTINUE**，未通过覆盖充分性判断。
- 已按Controller完成：finish-reading → 两综述相关章节的Initial Map → 27篇已有Primary证据的正式选择（复用合法卡，不重读）→ finish-reading → 修订Field Map → 实时重大分类审查请求 → 独立原生reviewer → 真实已完成子线程transcript attestation/提交 → QUERY_PLANNING → METADATA_RETRIEVAL → 补查 → PAPER_READING → 人工全文关口。Initial Map是本轮补成的认知交接，不能倒称它在上一轮阅读前已存在。
- 本轮补查13条查询，累计 **44条查询、396条候选/用户记录**；这些计数不是覆盖证明。旧筛查在新上下文显式沿用，未伪称补读了294条旧TITLE_ONLY排除的摘要。
- 四份PDF已正式登记为`USER_SUPPLIED_READ`，理由保留“用户在Claude Code中提供，为上轮临时检索发现、本轮正式检索未纳入的疑似遗漏”。结果如下：

|来源|正式paper_id|登记|本轮卡|
|---|---|---|---|
|Wei, TMECH2025, 10.1109/TMECH.2024.3443310|user_Wei_2025|成功|未出卡|
|Lei, Sensors J2025, 10.1109/JSEN.2025.3598958|user_Lei_2025|成功|未出卡|
|Zakeri, TIE2025, 10.1109/TIE.2024.3451138|user_Zakeri_2025|成功|未出卡|
|Chen, TII2025, 10.1109/TII.2024.3477549|user_Chen_2025|成功|未出卡|

**四来源未读原因是Controller接口阻塞：**登记后缺`identity_status`；`enrich-candidates`对`source_origin=user_supplied`直接跳过；`select-reading-subset`实际报`selected paper identity is not verified: user_Wei_2025`。未手改身份、状态或绕过准入，未重新注册别名规避。四份原PDF仍在`source-materials/user_{Wei,Lei,Zakeri,Chen}_2025.pdf`，无需重复下载；Wei/Lei的Claude笔记不等于本run接受卡。后续须由合法Controller入口解决身份/选择衔接，再出四卡。

## B. 本轮新读论文（PDF页序从1起；仅指定关键部分）

### 1. Optimizing Variable Admittance Control for Remote Ultrasound Scanning Under Uncertain Environment｜memfv6bMNssJ
- 问题：未知接触环境下远程超声的自由空间位置跟踪与接触恒力。
- 机制：线性弹簧模型；力与US置信度检测初始接触并记录位置，25Hz二阶Kalman在线估刚度；先调参考位置再调虚拟刚度，加补偿与能量罐。
- 数值/条件：PDF p1：体模/人体上肢误差范围0.2/0.4 N，SD0.02/0.2 N；p6刚度估计25Hz，初始0.3N/mm。
- 依赖：UR3e、ATI Nano25；US图像置信度；接触初始化和线性模型。
- 失效/边界：结构推断：p6记录pct后break，未验证持续更新移动环境位置；若位置变动，刚度估计可能混入运动。非作者明确静止假设或实测失效。
- 与本文：位置接口在线刚度+变导纳近邻；不能据此宣称已解决运动不变辨识和增益耦合。（PDF pp1,3–6）

### 2. A Passive Variable Impedance Control Strategy with Viscoelastic Parameters Estimation of Soft Tissues for Safe Ultrasonography｜ODwWwSVciuYJ
- 问题：异质胸腹组织安全远程超声。
- 机制：正弦触压离线最小二乘辨识HC黏弹参数，GP空间图；在线QP变刚度、力/穿透约束和能量罐。
- 数值/条件：p4：UR3e位置接口、BOTA SensONE与刚性球尖压头；beta1.35；有/无突然患者运动，VS-CF/VS-VF对CS/CF；未抄取精确误差数值。
- 依赖：离线参数图、穿透/力标定与空间对应。
- 失效/边界：p6作者局限：辨识要求患者久静止，需要在线估计；需要运动跟踪匹配模型变形，多患者验证尚缺。
- 与本文：直接近邻；纠正在线黏弹辨识说法，实际离线图+在线变刚度。（PDF pp2–6）

### 3. A shared-control ultrasound-guided robotic system for liver puncture under respiratory motion｜SS4c7arItfoJ
- 呼吸下共享控制肝穿刺，医生保留针朝向和插入。 位姿辅助+固定参数导纳的力调节接触稳定；集成探头和针操作机构。 p10：6N参考呼吸激活后0–12N描述在M/B/K调参之前，不能当控制后残余；p11称调节后稳定；p13犬在体为定性可行性。
- 依赖/边界：UR5e、US/CT、光学跟踪、人工操作与呼吸体模。 p13：样本有限，操作员依赖；犬在体针靶误差等未系统量化。 与本文：共享控制背景；不能用调参前波动否定固定导纳已实现的能力，也不证明运动/刚度分离。（PDF pp1–2,10–11,13）

### 4. Simultaneous Online Registration-Independent Stiffness Identification and Tip Localization of Surgical Instruments in Robot-Assisted Eye Surgery｜HokaNUBl3lAJ
- 眼科器械挠曲时在线定位尖端。 FBG巩膜力+双目挠度，AD/LSQ及免配准RIAD/RILSQ辨识器械3EI；力—变形模型构造KF位置量测。 p12：AD+KF挠曲区平均误差静态1.22mm、动态0.70mm；学习窗1/5s；pp11–12相机移动10mm及遮挡由既有数据合成。
- 依赖/边界：非零巩膜力/挠曲及视觉；辨识器械而非组织。 pp12–13：无挠曲不更新，适应增益速度/超调折中；颜色标记受光照影响，在体待验证。 与本文：外域的刚度增益—力换位置类比，非超声控制竞争者；所查内容未检验预测前馈耦合。（PDF pp3–9,11–13）

### 5. Toward Clinical Applications of Intelligent Robotic Ultrasound Systems｜AdPA6qQA2QMJ
- 问题：机器人超声的临床任务与工程转化障碍（Han等，RBME19,2026）。
- 机制：相关章节把接触力、示教/变阻抗、视觉与运动补偿、环境适应和软硬件共设计联系到扫描任务；不沿用单一综述分类。
- 数值/条件：PDF p9动态力控制长期成像；pp12–13呼吸/心跳、不可预测变形、图像退化和异构传感融合造成动态控制困难；综述不提供本轮独立实验。
- 依赖：综述层证据；仅审读接触力与运动/刚度变化相关章节及引用线索。
- 失效/边界：pp12–14：多数体模与少量健康志愿者，临床泛化/带宽与接口失配仍需验证；不是所有方法均失败。
- 与本文：交叉核对路线划分；提示外部/US运动观测、力反馈和硬件协同均需纳入。（PDF pp9,11–14,18）

### 6. Robot-Assisted Medical Imaging: A Review｜P4vOsQKTB1oJ
- 问题：成像设备精准位置与轨迹的机器人辅助（Salcudean等，Proc IEEE110,2022）。
- 机制：US章节并列力控制/被动顺应、感知刚度、图像视觉伺服及机器人弹性成像；明确体积US可用于器官运动补偿。
- 数值/条件：pp4–8：力反馈与深度信息、强度视觉伺服运动补偿、力/位置及刚度形变校正均已存在；无本轮独立实验数值。
- 依赖：仅相关US章节；跨成像综述，原始论文能力仍须核原文。
- 失效/边界：读者推断：图像质量/视野、接触力、外部感知和模型可靠性限制各分支；不能把综述引用直接变成原始控制证明。
- 与本文：三条路线是组织视角、相互重叠；不能遗漏早期图像运动补偿和机械安全。（PDF pp4–8,14–15）

### 7. Tactile-Guided Robotic Ultrasound: Mapping Preplanned Scan Paths for Intercostal Imaging｜txpI_b7Ta1IJ
- 肋间声窗狭窄下路径映射。 触觉信号学习提取肋骨点云、插值配准US模板，倾角补足骨下覆盖。 p1四种体模，路径MNND/HD3.41/3.65mm，重建相对CT0.69/2.2mm。
- 依赖/边界：p4 RGBD粗定位、人工选四角、预规划US模板、训练触觉骨识别。 推断：触觉几何映射不等同呼吸连续力补偿；近似平面与模板对应是边界。 与本文：弱相关；力觉能承载几何信息但未证明受控力/刚度运动观测。交付版本arxiv2507.20282v1，不能冒称正式终稿。（PDF pp1,4,7–8）

### 8. Ultrasound-Guided Assistive Robots for Scoliosis Assessment With Optimization-Based Control and Variable Impedance｜c2nzgpQiQ0IJ
- 问题：脊柱自动扫描多目标恒力/姿态协调。
- 机制：双环QP协调动力学与接触约束；专家示教GMM/GMR轨迹协方差逆映射平移和旋转阻抗。
- 数值/条件：pp6–7：体模及6名人，扫描3mm/s；体模7/15N、人体示例10N；未从图读出精确RMSE。
- 依赖：动力学模型、深度点云法向、US特征、专家示教分布。
- 失效/边界：p7：参数敏感性未研究；人体软组织与表面突起造成小振荡为作者解释。
- 与本文：真正的阻抗学习近年代表，但示教协方差不是在线环境刚度/运动分离。（PDF pp1,4–7）

### 9. AuSoScan: Automatic Scoliosis Assessment by Ultrasound Scanning With Soft Contact Control｜QPuQ4sfxSTIJ
- 问题：软组织扫描中的精确力/路径与脊柱重建。
- 机制：Hertz接触与速度相关前后接触半径模型，预先最小二乘识别k1/k2，再进入5步非线性MPC；US感知/重建。
- 数值/条件：p9：3/6/9/12N辨识数据；CasADi、5步预测，step/slope/wave对无软接触模型基线；表II为数值来源，未擅自OCR不可读值。
- 依赖：接触参数标定、软接触模型、机器人动力学、US/深度感知。
- 失效/边界：p10：严重脊柱侧弯/不平肩未验证，更强外部扰动反应留待未来。
- 与本文：模型+优化近邻；接触模型收益不能等同移动环境在线运动—刚度分离。（PDF pp1,5–6,9–10）

### 10. Dual-Agent-Based Robotic Ultrasound Path Planning and Interaction Control in External-Vision-Independent Environments｜D8AVZtwhD6gJ
- 问题：无外部视觉/几何先验的脊柱US路径与接触控制。
- 机制：DL定位脊突+Kalman跨椎间隙追踪；RL用6D力观测调方向并最小力矩；混合力到速度框架。
- 数值/条件：p1体模及人体：特征均误差<1mm，面外/内角误差1.61±1.1°/1.27±0.9°。
- 依赖：实时US、接触力、训练数据与解剖先验；p10需人工初始放置。
- 失效/边界：p10：交互控制仅用力观测、US可改善留未来；不能把整个系统称仅力（路径用US）。
- 与本文：路线2学习接触动作与图像路径混合；角度或特征误差不是组织运动估计误差。（PDF pp1,4–5,10）

### 11. Impedance Control on Arbitrary Surfaces for Ultrasound Scanning Using Discrete Differential Geometry｜i4wx6VF4kN0J
- 问题：曲面扫描中保持姿态/接触并限制穿透和失联冲击。
- 机制：以网格曲面有符号距离及参数化坐标构造被动阻抗；负目标距离产生接触力而非显式力前馈。
- 数值/条件：pp7–8 DLR MIRO硅胶体模：曲面参数RMSE0.917/1.465mm；空域距离0.871mm；负距离接触力SD0.498N，不是力RMSE。
- 依赖：已知网格/法向/主曲率及disk-topology参数化；选增益时利用已知硅胶刚度。
- 失效/边界：p8前馈力对照约7s失联并再冲击；推断：静态几何和负距离不能自动保证运动异质面精确恒力。
- 与本文：安全几何阻抗分支，不能将其作为已在线辨识组织运动的方案。（PDF pp1–2,6–8）

### 12. Optimization-Based Variable Impedance Control of Robotic Manipulator for Medical Contact Tasks｜GBHeodrYfZQJ
- 问题：静止/滑动医疗接触的力跟踪与被动性。
- 机制：N=1在线QP优化机器人刚度，力误差反馈与最小刚度正则，能量罐限制变刚度注能；低速忽略惯性。
- 数值/条件：p5：7/10/15N，8mm/s，软/硬/R2S/S2R；p6静态OS最大RMSE0.88N；p7扫描正文0.57N（摘要写0.5N，保留冲突）。
- 依赖：KUKA LWR4+、FT与FRI，机器人刚度100–500N/m；该范围不是本轮实测组织刚度。
- 失效/边界：作者验证含静态斜面/躯干不平；推断：不能当主动呼吸运动或运动不变刚度辨识证据。
- 与本文：优化变刚度直接近邻；输出机器人刚度不等同估计组织刚度，也不恢复组织运动。（PDF pp1,3–8）

### 13. Safe Contact Force Generation for Robotic Thyroid Ultrasound Imaging｜Ct8Urnw7rssJ
- 问题：颈部安全低力与诊断图像质量折中。
- 机制：串联弹性末端：线性伺服、弹簧/恒力弹簧提供低力；末端PID调压缩量，机械臂无需导纳反馈。
- 数值/条件：pp4–6：1–4N，线性伺服行程27mm；7人VAS；1–2N体模图像最接近手动，4N人体颈静脉完全压缩。
- 依赖：专用末端、力/位移测量、RGBD手工路径，受探头自重约束。
- 失效/边界：p7：尚需实时路径适应、考虑呼吸吞咽；速度对图像影响未评估。
- 与本文：机械顺应路线，机器人臂无导纳不等于全系统无反馈；低误差还须兼顾图像变形。（PDF pp1,4–7）

### 14. A Novel Ultrasound Robot With Force/Torque Measurement and Control for Safe and Efficient Scanning｜bGqa7m5UjxwJ
- 问题：安全恒力、快速调整成像深度并监测力矩。
- 机制：双并联弹簧/电机恒力机构，弹性变形光反射测力矩，临床经验力库调参考与深度。
- 数值/条件：p1：力/双矩误差3.53/6.68/6.11%FS，恒力误差<0.57N；p3力域4–12N；p10体模力误差最大0.48N。
- 依赖：定制末端和机械臂、标定、力库及医生图像判断。
- 失效/边界：pp10–11：力库8/10N为体模经验，实际临床库未建，部分在线流程使用假设数据。
- 与本文：机械顺应与主动调力代表；不是在线组织刚度辨识，受控探头运动传感也未证实。（PDF pp1,3,6,10–11）

### 15. A Unified Interaction Control Framework for Safe Robotic Ultrasound Scanning with Human-Intention-Aware Compliance｜AeZVQ15NBVsJ
- 问题：医生干预、患者移动和意外碰撞同时存在的扫描安全。
- 机制：分层阻抗+连续模式权重+能量罐；末端任务和零空间协同，按人意图区分跟随/恢复/扫描/避障。
- 数值/条件：pp6–7 Panda/ATI/AzureKinect颈动脉：患者推/转头、医生取凝胶与零空间碰撞演示；扫描约10N，跟踪误差范数<=0.04m为报告，不当精密恒力结论。
- 依赖：RGBD实时颈部跟踪与点云路径、Panda动力学和力矩控制、FT。
- 失效/边界：推断：运动适应是多模态分层控制，不能用此证明力编码器唯一传感或数十倍刚度鲁棒性。
- 与本文：安全/外部运动观察混合分支；患者运动已被覆盖，不能一概说先验全静止。（PDF pp1,6–7）

### 16. Diffusion-Based Impedance Learning for Contact-Rich Manipulation Tasks｜kexwZicrUfAJ
- 接触丰富操作中免逐任务手调阻抗。 力矩条件Transformer扩散重建sZFT平衡轨迹，能量估计器按方向调刚度阻尼；旋转使用SLERP噪声。 p14 TableIII三种插入各30次，全方法各100%；无方向适应各0%，无DM仅圆柱100%；pp1–2训练约1.6h。
- 依赖/边界：交付arxiv2509.19696v4(2026-08-11)预印本；iiwa实时力矩、示教与FT；未证明超声任务迁移。 p15：仅对角半正定刚度；更换示教基准阻抗待验证；降刚度被动，升刚度一般不保证被动。 与本文：外域学习阻抗代表；学习平衡轨迹+解析适应，非组织刚度/运动辨识。（PDF pp1–2,13–15）

### 17. Quasi-Static Elastography-Driven Automated Robotic Ultrasound Screening and Localization｜Vzv82t4dJgsJ
- 自动弹性成像筛查/异常定位所需的受控压缩。 顺应力控结合面外姿态/面内触压；B模无监督位移与多帧应变融合成3D应变图闭环定位。 p1两体模：腹部检出0.77/定位1.01±0.47mm，甲状腺0.73/3.44±0.84mm；人体可行性。
- 依赖/边界：图像灰度异质性和准静态压缩；B模替代无RF输入。 p13作者：不适用于移除强度对比的等回声模拟，设备处理后仍需稳定灰度形态；非动态呼吸刚度真值。 与本文：功能成像分支；应变对比与控制所需绝对刚度、运动分离应区分。（PDF pp1–2,11,13）

## C. Field Map 摘要与覆盖审查原文要点

**①激励/模型辨识刚度。** 接触模型与激励/标定得到材料信息，再用于控制适配，代表为Beber ICRA2024、Xiao Access2025、AuSoScan TMECH2025；Duan RA-L2022的示教协方差学习阻抗是相邻信息分支。依赖可用接触几何、材料模型、参数图或足够辨识信息。Beber明确局限于离线图及运动匹配不足，Xiao已检段落未验证接触点随运动更新。原因是力同时受组织运动、形变和材料影响，缺少分离信息时参数可能混淆；这是有条件推断，不能说所有估计器均已失败。

**②力反馈补偿运动。** 显式力、位置导纳、学习接触动作和机械顺应可直接抑制力偏差，代表为Piedra IJRR2026、Fu TIM2024、Wang TASE2025、Bao TIM2023、Tsumura RA-L2024、Wu TRO2026。依赖力信息和执行带宽，力矩控制需动力学/接口，机械方案需专用结构及工作域。快速运动、内环延迟、机械行程/速度边界或成像变形会限制效果。原因是这一路线改变扰动传递与反馈通道，无须先恢复组织位移，但仍受执行链和安全/图像质量折中约束。

**③观测运动再补偿。** 外部视觉、US图像、独立压头及受控探头力觉可提供不同运动线索，已有代表为Jiang TIM2022重定位、Qian ICRA2026高速US伺服、Yan IROS2024外部颈部跟踪。依赖观测可辨识性、标定、视野与时延，两个力觉子支在本run仍缺直接原始锚点。外部重定位限可见刚性变化，US方法可因旋转/剪切歧义失稳；Ebrahimi只能作器械刚度转位置的外域类比。原因是低时延不能消除观测歧义，受控探头同时测量运动时转换增益还参与闭环，需独立原始证据。

三路线允许交叉；学习动作、学习阻抗、辨识组织刚度不是同一事。未据本图裁决本文新颖性或宣布“无人做过”。

独立原生`coverage_reviewer`请求`930c89f6bafd4cc585debf8181cddcf4`，结论原文：**CONTINUE**。实际已完成子线程由`transcript_verifier`证明并被Controller接受，未声称Hook执行；审查绑定的是补查前的29卡/347候选快照。

> 六项绑定文件的SHA256均与Controller请求一致；31条查询按综述、骨干、前沿、定向补查顺序实际执行。

> 三路线允许交叉，区分环境刚度、机器人阻抗、运动观测和执行接口；现有证据支持局部机制比较。

> 347个去重候选中294项仅题名筛查、6项身份待核、5项缺真实摘要，尚未满足逐项题名加摘要筛查要求。

> 29张接受卡不能补足尚缺的历史机制和直接近邻；当前地图准确标为PARTIAL，不能判定覆盖充分。

其五项遗漏原文：

1. 恢复必要历史机制证据：独立压头与同一受控探头的力源运动观测均缺直接原始锚点，当前接受骨干未连接这些机制与近期动态异质接触任务；须正式定向检索、筛查并阅读必要历史论文。
2. 正式纳入并核对Wei、Lei、Zakeri、Chen近邻，以及Li H.-Y.移动环境变导纳、Jiang Actuators2024力跟踪、Tan TASE目标跟踪、Jiang RA-L2022关节运动；这些证据可能改变R1/R2/R3边界和当前前沿判断。
3. 补齐候选摘要筛查：294项TITLE_ONLY排除、6项HOLD_IDENTITY和5项缺摘要记录不能算完成题名加真实摘要审查；优先复查可能涉及运动辨识、接触控制及预测补偿的候选，并保留无法完成者的覆盖限制。
4. 核实六锚点的实际前向引文扩展：Yuen查询当前仅为身份检索，其余五条查询记录的cited_by为DOI；须用provider执行证据确认真正取得引用列表，并对新增分支继续必要分页，不能将精确题名检索当作前向扩展完成。
5. 补足互补综述与骨干覆盖的科学说明：当前两张接受综述卡均偏机器人医学成像，未充分展示非静止参数辨识与力源生理运动补偿的互补机制来源；应说明已撤回文献造成的具体覆盖影响，并通过适当替代证据或必要原始链补齐。

演进评估原文状态：`foundation_to_frontier: GAP`；`key_nodes_and_branches: GAP`；`transition_causality: PASS / DECLARED_TRACES_REVIEWED`；`explanatory_coherence: GAP`。其中PASS只针对已声明T1–T3，不是总体覆盖PASS。

补查尚未推翻上述结论：四来源只登记未读，六锚点本次仍是身份/题名核查，并未补成真正的provider引用列表；旧TITLE_ONLY、身份与摘要缺口仍保留。

## D. 被引线索核查

1. **Li H.-Y.等，Stable and Compliant Motion of Physical Human–Robot Interaction Coupled With a Moving Environment Using Variable Admittance and Adaptive Control。** 已经正式Scholar查询和Controller metadata enrichment核实：2018，IEEE Robotics and Automation Letters，DOI **10.1109/LRA.2018.2812916**；候选`xbXpC4gQiEAJ`已准入并选全文。真实摘要表明同时处理人、机器人及移动环境，变导纳+自适应产生时变交互刚度并讨论被动性，属直接相关前序；原始方法/数值/运动观测细节尚未读，不能先断言解决或未解决本文问题。

2. **Jiang J.等，Force tracking control method for robotic ultrasound scanning system under soft uncertain environment。** 本次正式完整题名查询0条对应命中。上一轮Piedra参考文献线索为Actuators 2024, 13:62，DOI 10.3390/act13020062；本轮未独立核实其身份与IAAC机制，不将引文线索冒充正式核验或Evidence。

3. **Tan等，Autonomous trajectory planning for ultrasound-guided real-time tracking of suspicious breast tumor targets。** 本次完整题名查询只回一条不相符记录，未核实。既有引用线索为TASE21(3):2478–2493，在线2023/卷期2024尚待核；DOI本轮未确认，连续US目标追踪能力待补查。

4. **Jiang Z.等，Towards autonomous atlas-based ultrasound acquisitions in presence of articulated motion。** 本次完整题名查询0条对应命中，未核实；两篇已读来源指向RA-L2022, 7(3):7423–7430，DOI未确认。需后续合法Query Plan用题名片段/作者/引用扩展继续核查，不能把未命中解释成论文不存在。

此外，Moreira2014、Haninger2024查询返回了Scholar引用cluster标识，但本次未执行真正的`cited_by`列表；不能将这些身份结果算作前向扩展。

## E. 当前唯一人工全文批

- **完整题名：** Stable and Compliant Motion of Physical Human–Robot Interaction Coupled With a Moving Environment Using Variable Admittance and Adaptive Control
- **作者：** Hsieh-Yu Li; Ishara Paranawithana; Liangjing Yang; Terence Sey Kiat Lim; Shaohui Foong; Foo Cheong Ng; U-Xuan Tan。
- **年份/出处：** 2018，IEEE Robotics and Automation Letters。
- **DOI：** 10.1109/LRA.2018.2812916。
- **目录：** `D:/pip_env/temp/scan_impedance_literature_20261031/aris_run/source-materials/`。
- **建议文件名：** `Li_2018_Stable_and_Compliant_Motion_Moving_Environment.pdf`；paper_id=`xbXpC4gQiEAJ`。
- 放置全文后应通过一次`submit-human-fulltext-batch`清单交付；尚未提交空清单或声称文件已到。该人工批只含Li，其他身份未核线索不伪造DOI或下载条目。

四用户来源已有文件，不在重下载清单。当前人工关口与四来源的Controller身份/选择接口阻塞是两个不同问题。

全部本轮草稿、绑定、查询结果、卡片与审查原文在`RUN_ROOT/read2_work/`；正式卡和Map由Controller写入规范目录。未修改ARIS源码或项目`E:/github_project/scan-impedance-update`文件，未执行Git pull/commit/push。
