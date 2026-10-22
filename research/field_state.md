# 运动与非均匀组织扫描中的接触估计与变阻抗控制研究现状

更新：2026-10-22。

机器人超声扫描需要在组织表面运动、局部力学性质变化及探头持续移动的条件下保持稳定接触。组织运动改变接触位置和相对速度，材料非均匀性则改变局部力—位移关系；二者均可能引起接触力波动。围绕这一问题，已有研究主要形成了自适应力跟踪、模型预测变阻抗、学习阻抗、生理运动补偿以及接触几何与力学参数联合估计等路线。当前讨论的研究方向是利用力与探头运动信息联合估计局部接触状态，并协调阻抗调整和参考补偿。以下从这些直接相关的原始研究展开，区分已有方法实际解决的问题与本项目尚未建立的技术进展。

## 1. 未知环境下的自适应与最优阻抗控制

针对环境刚度和接触位置不确定引起的力跟踪误差，早期研究通过自适应反馈或参考修正获得期望接触行为。Jung等TCST2004（E23）、Duan等RAS2018（F02）及Song等RA-L2025（E04）代表了这一类方法。它们并不都依赖完整恢复组织力学参数，而是直接利用交互误差修正控制输入。Roveda等RA-L2016（E44）进一步采用扩展Kalman滤波估计接触刚度，并更新外层导纳的LQR控制参数，将在线环境认识与控制器调整连接起来。[原文E23](../literature/papers/field_papers/E23.pdf)、[F02](../literature/papers/field_supplements/F02.pdf)、[E04](../literature/papers/field_papers/E04.pdf)、[E44](../literature/papers/field_papers/E44.pdf)。

最优阻抗研究也已经考虑未知环境动力学及接触位置。Ge等IJC2014（F03）研究未知环境下的最优交互；Liu等TCST2021（E45）将环境位置或轨迹与动力学参数一并纳入系统表示，通过输入输出历史重构及自适应动态规划求取阻抗。后者采用线性环境模型及相应学习条件，不能直接视为任意时变软组织问题的完整解，但它说明“现有最优阻抗统一依赖已知静态表面”的概括不准确。[原文F03](../literature/papers/field_supplements/F03.pdf)、[E45](../literature/papers/field_papers/E45.pdf)。

对非线性软接触，Kong等TSMC2025的NNBO采用Hunt–Crossley接触辨识、critic近似和Newton位置参考更新，并在软材料上验证最优交互控制。其材料预辨识、表面信息和位置参考接口构成了方法的具体条件。它与直接选择刚度、阻尼的算法不同，但仍是软接触优化控制的重要近邻。Pollayil等T-RO2023则从惯量、扰动及容许误差出发规划刚度与阻尼，并在特定结构下给出解析结果，表明阻抗优化的成本依赖于问题结构。[NNBO原文](../literature/papers/NNBO.pdf)、[Pollayil原文](../literature/papers/pollayil.pdf)。

## 2. 模型预测变阻抗与在线接触模型更新

Jin等TMECH2023（E41）将阻抗选择写成预测优化问题：单步形式采用QP，多步形式在预测时域内处理精度与柔顺之间的权衡，并结合增益和能量条件。Haninger等RAS2023利用Gaussian process描述交互力及不确定性，联合规划参考与阻抗。因此，Haninger原文的主要能力是概率模型支持的预测交互控制，不能用本项目采用的阻尼匹配规则概括整篇方法。[Jin原文](../literature/papers/field_papers/E41.pdf)、[Haninger原文](../literature/papers/haninger.pdf)。

Xue等RAS2025的MPVIC进一步结合在线HC参数估计、线性化接触预测、性能/保护模式和能量调节，并在非均匀软材料和人体交互中验证。其MPC控制量是误差运动的jerk，刚度和阻尼由模式与能量机制调节，并非任意任务代价下直接自由优化全部增益。Xue等IROS2025（E34）又引入学习调节和预定时间交互设计。两项工作已经把接触适应、响应速度和执行限制结合起来。[RAS原文](../literature/papers/xue.pdf)、[IROS原文](../literature/papers/field_papers/E34.pdf)。

本次补入的两篇直接相关工作扩展了在线预测适应的比较范围。Mu等CEP2025（G01）在自适应阻抗中引入动态观察器，通过历史误差滚动优化观察器权重，再由MPC生成速度补偿。Wang等TIE2026（G02）根据历史力预测误差在线调整阻尼和环境刚度模型参数，并用于下一控制周期的力预测和运动补偿。前者主要优化观察器参数，后者主要校正预测模型；二者均说明“历史接触信息—在线参数更新—预测控制”已有完整实现。[Mu原文](../literature/papers/field_supplements/G01.pdf)、[Wang原文](../literature/papers/field_supplements/G02.pdf)。

移动时域估计与变阻抗的结合也已有明确先例。Zhang等TIE2024提出机器人—环境观察器，利用历史信息同时估计接触环境阻抗参数与机器人内部扰动，再用于变阻抗加工控制。官方摘要及引言已核实这一关系，完整方程尚未取得。因此，可以确认MHE驱动VIC已有研究，但不能进一步断定该文是否覆盖了运动组织与材料非均匀性的具体分离机制。[官方来源](https://doi.org/10.1109/TIE.2023.3292856)。

在其他交互任务中，Minniti等ICRA2021（G25）将在线环境辨识或自适应控制与MPC结合，估计结果用于更新预测模型，并通过开门和负载操作验证。该文的任务及动作接口不同于组织扫描，但为在线辨识和预测执行的协调提供了直接方法依据。[原文](../literature/papers/field_supplements/G25.pdf)。

## 3. 学习阻抗与参考的联合调整

学习方法不统一等价于在线估计刚度后按规则设定增益。Li与Ge的TCST2014（G06）通过重复交互学习未知时变环境中的阻抗；Peng等TNNLS2022（G15）采用critic学习最优导纳；Roveda等AI2022（G08）通过actor–critic近似预测控制，在线调整参考与阻尼。Anand等RAS2023（E40）则学习机器人笛卡尔阻抗动力学，再以CEM进行预测刚度选择。这些方法分别在学习对象、训练方式及输出变量上作出不同设计。[G06](../literature/papers/field_supplements/G06.pdf)、[G15](../literature/papers/field_supplements/G15.pdf)、[G08](../literature/papers/field_supplements/G08.pdf)、[E40](../literature/papers/field_papers/E40.pdf)。

Li等T-RO2025（E27）针对未知地形上的持续接触力跟踪，结合随机地形预训练、变阻抗反馈与参考修正，包含静态和动态交互实验。Bogdanovic等RA-L2020（G23）及VICES IROS2019（G24）更早研究了同时输出期望位置与阻抗的学习接口；Zhang等RA-L2021（G07）将逆强化学习用于力相关任务的阻抗策略。这些研究说明，参考和阻抗的共同调整已有广泛基础，当前方法的差别需要落实到具体估计模型和控制机制。[E27](../literature/papers/field_papers/E27.pdf)、[G23](../literature/papers/field_supplements/G23.pdf)、[G24](../literature/papers/field_supplements/G24.pdf)、[G07](../literature/papers/field_supplements/G07.pdf)。

## 4. 生理运动补偿与机器人超声扫描

生理运动补偿已经直接进入接触力控制。Moreira等CMPB2014（F10）研究黏弹组织模型、主动观察器与生理运动补偿，并包含组织材料测试和呼吸运动模拟。Cao等TIM2025（G13）通过深度相机和AprilTag监测运动，用自校准模型预测下一周期接触位置，再将预测位置引入力传感导纳控制；其三维成像模块另外补偿运动造成的体积失真。该文已明确建立“运动预测—参考补偿—力反馈”的控制结构，所核控制器中的运动信息来自外部视觉。[Moreira原文](../literature/papers/field_supplements/F10.pdf)、[Cao原文](../literature/papers/field_supplements/G13.pdf)。

Velikova等ICRA2024（G14）处理的是另一层问题：在法向恒力采集和US图像引导下，利用隐式神经表示重建体积的强度和语义信息，降低呼吸造成的重建伪影。它与Cao的运动补偿控制相关，但不属于局部组织刚度估计或变阻抗决策算法。用户文字线索中的Li等RA-L2024则属于呼吸运动下的胸腹穿刺力—位置控制；目前仅完成文献身份核对，未据题名推断其组织参数估计方式。[Velikova原文](../literature/papers/field_supplements/G14.pdf)、[穿刺论文](https://doi.org/10.1109/LRA.2024.3386017)。

对于刚度与表面位置共同变化，DOSMAC CEP2026直接将二者建模为总扰动，通过扰动观察器和滑模导纳进行补偿。出版社摘要及贡献段明确报告了这一能力。这意味着“必须先分别恢复刚度和运动才能控制接触力”不是已有研究支持的普遍前提。其完整方程仍待取得，尚不能作为已复现算法进行精确数值比较。[官方来源](https://doi.org/10.1016/j.conengprac.2026.106764)。

任务层面的强相关研究还包括：AuSoScan TMECH2025（F08）的软接触MPC扫描；Wang等TASE2024（F09）的六维力/矩调节；Qian等TIM2026（G10）的历史表面几何更新、姿态调整及图像/力融合目标力；Ning等TBME2023（G17）和Wang等TASE2025（G16）的US图像与力信息引导。这些方法已经将接触调节用于实际扫描任务。无外部视觉的系统仍可使用US图像，不能等同于仅使用力与编码器的控制器；其成像结果也不能由本项目的力跟踪仿真替代。[F08](../literature/papers/field_supplements/F08.pdf)、[F09](../literature/papers/field_supplements/F09.pdf)、[G10](../literature/papers/field_supplements/G10.pdf)、[G17](../literature/papers/field_supplements/G17.pdf)、[G16](../literature/papers/field_supplements/G16.pdf)。

## 5. 接触几何、柔顺参数联合估计与稳定执行

Sanan等RSS2014的CARE（G21）通过接触力和探触位置联合估计组织柔顺图与配准关系，并分析可观测性。其模型采用先验表面、局部小变形和相应接触假设，实验利用不同力水平提供辨识信息。Chalasani等ICRA2016与RA-L2018进一步研究连续触诊中的组织几何和刚度在线更新，后两篇目前仍为机构摘要证据。由此可见，联合估计几何与力学参数并非新提出的思想；但静态配准或连续触诊的结果，也不能直接扩展为运动且非均匀组织上的闭环扫描保证。[CARE原文](../literature/papers/field_supplements/G21.pdf)、[ICRA2016](https://doi.org/10.1109/ICRA.2016.7487609)、[RA-L2018](https://doi.org/10.1109/LRA.2018.2801481)。

稳定执行方面，Beber等ICRA2024结合软组织黏弹参数估计、QP及能量池；Jin等TMECH2025的SNEUM（F06）采用神经势能表示状态相关刚度；Pan等RA-L2024（G05）通过复合学习处理机器人动力学不确定性并分析变阻抗稳定性和无源性。Pan的主要辨识对象是机器人动力学，不能解释为已经实现组织运动辨识。[Beber](../literature/papers/beber.pdf)、[SNEUM](../literature/papers/field_supplements/F06.pdf)、[Pan](../literature/papers/field_supplements/G05.pdf)。

Thelenberg与Ott的ECC2024（G04）在有限时域内优化刚度与阻尼，并利用终端条件处理时变阻抗；Poignonec等RAS2026（G03）将增益动态与无源约束同时纳入NMPC。后者采用给定目标阻抗轮廓、关节运动及增益跟随速率优化，所核预测模型保持时域内外力不变，并讨论使用环境模型改进预测的可能性。它说明预测与稳定执行已有专业完整的构造，不能把加入无源约束本身列为当前方法的独立创新。[G04](../literature/papers/field_supplements/G04.pdf)、[G03](../literature/papers/field_supplements/G03.pdf)。

原文中的计算成本也有明显差异：Mu报告所述MPC部分约13 ms；Thelenberg在七自由度仿真中每80 ms执行优化；Poignonec在所述UR5配置中以8 ms采样、80 ms预测窗口运行，包含NMPC的控制循环平均约5.6 ms。这些结果对应各自硬件、预测规模和实现，支持具体成本取舍的比较，不能直接证明本项目相对整个领域的速度优势。

## 6. 最接近研究与当前方法的具体关系

| 最接近工作 | 输入与估计内容 | 对控制的实际作用 | 与当前思路的具体关系 |
|---|---|---|---|
| CARE G21；Chalasani ICRA2016/RA-L2018 | 接触信息用于柔顺/几何或配准估计 | 接触探索、几何与力学认识 | 联合估计有直接先例；静态几何/配准与运动表面的区别须按具体模型说明 |
| Liu TCST2021 E45 | 未知环境动力学及位置/轨迹，历史重构 | 最优阻抗适应 | 不能把未知环境位置说成领域未处理的条件 |
| Zhang TIE2024 MHE-VIC | 历史信息、环境阻抗及机器人扰动 | 观察器驱动VIC | MHE与VIC结合已有先例；完整方程尚未核读 |
| Xue RAS2025；Mu CEP2025 G01；Wang TIE2026 G02 | 在线接触模型、误差历史或预测模型参数 | 预测运动控制及阻抗/补偿更新 | 是模型估计与预测控制路线的直接比较对象，不能只比较未更新模型 |
| Haninger RAS2023；Roveda AI2022 G08 | 概率力模型或学习交互模型 | 参考与阻抗/阻尼协调 | 联合调整参考与阻抗已被研究 |
| Cao TIM2025 G13 | 外部视觉运动预测与接触力 | 导纳参考补偿及力反馈 | 运动预测用于参考补偿已有实现；当前思路拟使用不同观测条件 |
| DOSMAC CEP2026 | 刚度/平衡位置变化的总扰动 | 扰动补偿与滑模导纳 | 不显式分离两种来源也能处理共同变化；全文和同任务比较仍待补足 |
| Beber ICRA2024；SNEUM F06；Pan G05；Poignonec G03 | 材料、状态或增益动态，视方法而异 | 稳定/无源的阻抗执行 | 为执行设计提供直接基础，不能将一般稳定性机制本身视为新的任务贡献 |

综上，现有研究已经能够在线适应接触模型、联合调整参考与阻抗，并通过视觉运动预测或总扰动补偿处理动态交互。当前提出的“从力与探头运动历史联合估计局部刚度和组织运动，分别用于阻抗调整与参考补偿”应在上述工作的基础上定位。其与视觉补偿的区别在于观测来源，与总扰动补偿的区别在于对接触响应的分解及其使用方式，与MHE-VIC和预测适应的区别则必须落实到实际估计模型和控制律。**目前文献依据不足以支持“首次联合估计运动与刚度并用于变阻抗”的结论，也不能仅凭这些模块的组合认定创新性成立。**

修复前C1 v2仿真可以作为方法线索，但不能用于定义领域空白。组合扫描中C1的RMSE约0.147 N，去掉漂移通道后约1.502 N，只保留增益变化约0.489 N，说明参考补偿及接触漂移建模对该实现有明显作用。与此同时，C1在硬组织推拉中的RMSE约1.372 N，仍高于估计＋Haninger阻尼规则适配的0.433 N。上述数值来自本项目的同任务适配，而不是相应原论文的能力上限，尤其“Haninger规则适配”不等于完整GP预测阻抗算法。

此外，C1代码中的漂移项是等效接触响应：材料变化、黏弹松弛和表面运动均可进入该项，力水平参数也吸收偏置。它没有单独估计绝对物理表面位置。因此，现有消融支持该通道的控制作用，尚不能证明已经准确分离了两个物理来源。新问题和新方法的贡献表述应与实际估计量保持一致。[实验来源](../experiments/cobot_scan_benchmark/completed_tasks/README.md)、[C1模型定义](../experiments/cobot_scan_benchmark/candidate1_model.py)。

当前控制实例已完成一次针对该缺陷的修复：先从刚度回归中消除局部三次漂移趋势，再更新刚度，并改用刚度相关的接触恢复。在原三个测试任务上，硬组织推拉RMSE降至0.291 N，20个硬组织场景均未记录到低于0.5 N的接触力；组合扫描为0.142 N，但静态刚度变化误差从0.060升至0.087 N。相同恢复与求解次数、改回原估计器后，硬组织误差为0.531 N，支持估计器替换的额外作用。该结果改进了原C1的实际能力，不改变上述已有工作的充分能力，也不建立真实组织运动恢复或完整新颖性结论。[修复方法](optimal_impedance_question_and_methods.md)、[原始结果](../experiments/cobot_scan_benchmark/contact_update_results/README.md)。

## 文献范围与来源

本次以用户指定的两个本地资料库为主要来源，结合出版社、作者稿和官方元数据补查。补入22份经过用途筛选的原文或作者稿及对应提取文本；其中直接相关工作用于上述现状，少量相邻任务、综述及方法资源单独标注用途。隐式策略适应G09、敏感度MPC G12和其他通用优化机制不被当作已经解决本扫描任务的直接证据。政策、斑块分类、无关注册和仅用于外域样本校准的论文未全部纳入。

原文版本、DOI、哈希、实际核读页码及尚未取得全文的条目见[本轮来源](../literature/motion_stiffness_review_sources.json)和[统一文献清单](../literature/local_materials_manifest.json)，完整文件入口见[文献索引](../literature/literature_index.md)。Cao与Velikova的方法页已视觉核对；选入和保存全文不等于全部通读或独立复证定理。Zhang MHE-VIC、DOSMAC和Chalasani等的全文缺口仍明确保留，未据此宣称其缺少某项能力。该现状是围绕当前问题的定向文献分析，不是全部领域的系统综述或新颖性不存在证明。
