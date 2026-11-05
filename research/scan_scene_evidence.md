# 超声恒力扫描：仿真环境与证据

2026-11-05。实现在 `experiments/rossi_failure_study/physics_plant3.py`、`scenes3.py`。验证：`check_curve3.py`、`check_plant3.py`、`smoke3.py`；对真实数据的核查：`real_data_audit.py` 和 `project_records/agent_runs/2026-11-05_codex_real_data_audit/final.md`。

研究范围由任务决定，不由数据来源决定。代表性难点是软硬突变（肋骨、硬夹杂）、组织运动、曲面和突然推拉。现有肋骨和骨脊工况能覆盖的，不另建场景。

## 仿真真值怎么来

| 部分 | 来源 | 性质 |
|---|---|---|
| 半球压头上的静态力—深度 | Remus S1 逐人逐部位拟合（`remus_s1.json`），r = 10 mm，约 5 mm/s | **人体在体实测** |
| 个体和部位差异 | 按部位（默认 R1–R5）抽受试者，用整条曲线，同一人放松和主动曲线配对 | **人体在体实测** |
| 从压头到探头 | 等应变 s = 深度/√面积 下，切线按面积^γ 换算，γ = .5（半空间）、1（薄层） | 物理推算，未被真实探头数据验证 |
| 探头几何 | 凸阵 69×17 mm、R60；线阵 48.43×11.56 mm 平底 | 尺寸来自规格；同一探头的完整组合是理想化 |
| 瞬态肌肉收缩 | 力在同一人的放松和主动曲线之间插值，活化水平为 a | 两条曲线实测；时程和强度是假设 |
| 动态响应 | 总相位 7、11、17、23°，由 HC 加记忆、弱幂律、Kelvin-Boltzmann 三种机制之一产生 | **假设**，人体腹部没有实测 |
| 骨 | 覆盖层 7、18、32 mm，肋间隙 6–19 mm | 覆盖和间隙来自尸体实测；肋宽、边缘形状是假设 |
| k_rib | 10–45 N/mm | 假设，整体胸廓约 9、单肋约 45 N/mm 只见检索摘要 |
| 扫描运动 | 速度 1–60 mm/s，加速度 .05–.5 m/s² | 速度有人工扫查实测；加速度是假设 |
| 呼吸 | 峰峰 3–12 mm，周期 3–6 s | 中心值来自人体在体测量，范围端点是假设 |
| 推拉、惊跳、咳嗽、心动、震颤 | 见 `scenes3.py` | 假设，作为任务压力情景 |
| 摩擦、切向柔顺、机械臂柔顺、身体支撑 | μ .03–.30，K_t 5–50 kN/m，kr 50–200 kN/m，m 15–40 kg、7–10 Hz | 假设；重复精度和可编程阻抗都不是刚度实测 |
| 测量 | 1 kHz 白噪声 .08 N，30 Hz 一阶滤波后输出约 .025 N | 假设。按 FT300 手册折算约 .18 N（输出约 .055 N），只作敏感性 |

未知量逐场景按范围抽样，可用 overrides 固定。结果按因素分开看，不把任一未证实的取值当真实场景。

## 对象给出的工作点（推算，不是探头实测）

全部个体曲线，8 N，P05–P95：

| 接触律 | 探头 | 放松：压深 mm，切线 N/mm | 主动：压深 mm，切线 N/mm |
|---|---|---|---|
| 局部支撑 | 凸阵 | 7.6–12.4，.84–1.59 | 3.1–8.5，1.42–3.88 |
| 局部支撑 | 线阵 | 9.3–18.3，.50–.79 | 2.4–10.5，.91–3.31 |
| 半空间 | 凸阵 | 11.8–21.2，.49–.81 | 4.0–13.4，.85–2.48 |
| 半空间 | 线阵 | 14.2–23.2，.43–.68 | 3.7–13.8，.71–2.15 |

全部工况和组织跑一遍（只用 L，不是比较）：324 次，数值失败 0。力峰值超过 15 N 45 次（肋骨工况 21、主动 12、混合 10、放松 2），超过 36.5 N 1 次（肋骨工况）。探测处增益、相位和 RMSE 的分布在 `results/smoke3.json`。

## 哪些得到了独立校核

| 检查 | 结果 | 说明 |
|---|---|---|
| 半空间换算对精确弹性解 | Hertz 球换算到圆形平压头，与精确解 k = 2E\*a 一致（误差 1e−16） | 数学上验证 γ = .5 的线弹性半空间极限，不验证组织是半空间 |
| 曲线隐含模量对 Remus 反演结果 | 5 mm 处隐含 E\* 为 12–19 kPa；Table 6 的 Ogden c = 5–8 kPa，按 FEBio 的 Ogden 形式小应变 μ = c/2，E\* ≈ 2c = 10–16 kPa | 量级一致。同一实验、不同分析路径，所以不是独立实验 |
| 曲线拟合质量 | S1 个体拟合 RMSE 中位 .40 N，P95 1.86 N | 真值曲线本身的不确定性 |

**没能验证的：**
- **压头到探头的绝对迁移。** Dahmani 是唯一的真实探头数据，但对象是明胶，缺模量和探头规格，0 mm 处约 6 N 的偏置原因不明，指数对零点极敏感（2–7）。它只说明真实探头记录有明显的非线性和零点问题。
- **人体腹部 1–50 Hz 的动态响应。** 见下。

## 数据各能说明什么

| 数据 | 能说明 | 不能说明 |
|---|---|---|
| Remus S1、S3 | 该压头、该加载速度下的曲线，个体和部位差异，主动变化 | 探头绝对曲线，平衡弹性与动态部分的比例（只有单一速度），空间相关场 |
| Dahmani | 见上 | 模量、接触零点、腹部参数 |
| Greiner（离体人脑） | 松弛跨多个时间量级、循环有明显耗能和预调节，作为机制形式 | 腹部参数，1–20 Hz 相位，速率指数 |
| 猪脊髓 Rheometer | 预调节趋势 | 处理后没有物理时间和闭合回线，不能算滞回和相位 |
| Callera（离体胰腺） | 汇总的 tanδ 中位 .11、范围 .01–.30 | 频率和振幅未给，不能放到 20.8 Hz |
| Boyer（皮肤，只有摘要） | 10–60 Hz 小幅响应的边界 | 腹部，逐频相位 |
| Beber（胸部仿体） | 力—位移形状，指数约 1.24 | 图中没有逐点时间，不能得出保持段松弛量和相位 |
| USForce（1 段约 4 s，30 FPS） | 这一段运动的量级 | 单位和力的方向只能推断；不足以决定默认范围；不能辨识压入、法向力或 20.8 Hz |

**动态的局限：** 现有三种机制在同一协议下覆盖不了脑组织的慢松弛和回线形状，时间常数都不超过 .8 s。这不说明腹部需要脑组织的常数，所以不移植。结论只限于这个任务的时间尺度，慢尾在 10 s 场景内近似一个缓慢扰动。7–23° 是探索档，不是腹部实测范围。现有数据既不能证明也不能排除 HC 速率项。

## 已知的简化和未解决项

- Remus 曲线是约 5 mm/s 加载下测的，这里当作弹性曲线再叠加记忆，可能部分重复计入。单一速度分不开，未做校正。
- Remus 含腹壁和身体的整体响应，再串联身体支撑和局部背衬，可能部分重复计入（支撑刚度远大于接触刚度，静态影响小）。
- 5 mm 下限是正则化假设：多项式在 5 mm 以下斜率不随接触面积归零。约 61% 的工作点落在这一段。
- 体表起伏改变压入，不改变无跟踪时的局部法向，姿态角只含宏观曲率和误差。
- 平底对齐初触时切线有阶跃，力连续，这是该理想几何的结果。
- 4 mm 覆盖只作端点，回放较慢（约 1.2 s/仿真秒）。
- 半空间闭式力在凸阵加覆盖层时与 ∫k 差约 13%，真值切线取力的实际斜率。

## 接触公式

局部压入：δ(ξ) = d − g(ξ) − ξ tanθ − ξ²/(2R_s)，线阵 g = 0，凸阵 g = ξ²/(2R)。θ > 0 表示表面前方下降：压力后移，摩擦使轴向力增大。

参考力：局部支撑 p = qδ/(1−δ/h)，覆盖层硬化到局部刚度达到肋骨和胸廓上限后线性延伸；半空间 k = 1.1E\*√A，F 为闭式。空间场只含相对变化，均值为 1。参考力经单调映射（PCHIP，切线严格一致）变成该受试者曲线给出的力。`q_scale` 作用在映射之后，是真正的整体倍乘。

动态：F = (f_e + Σq_i)(1 + λḋ)，Prony 状态因果更新；λ 在 8 N 均匀软组织上由 φ_HC 确定后固定。

观测：策略默认只拿到测量通道，真值只在 `episode(..., diagnostic=True)` 时提供，只用于诊断。Method 控制轴向力；评价分别报告轴向和法向。

## 来源与原文位置

| 标识 | URL或DOI、原文位置 | 原始数据与访问边界 |
|---|---|---|
| E1 | [GE探头表](https://www.gehealthcare.com/-/media/gehc/us/files/products/ultrasound/versana/versana-premier/transducer-guide-versana-premier-r3-july-2024-doc2993174.pdf)，p.1；[Mindray探头表](https://www.mindray.com/content/dam/xpace/en_us/products-solutions/products/ultrasound/point-of-care/tex/5016/TE_X_transducers_40630_Rev_3.0.pdf)，p.1 | 沿用 `final_web.md`；`final_web_verify.md`记有502/403，未完成独立复核；型号间数据不拼作同一探头规格 |
| E2 | AuSoScan [10.1109/TMECH.2025.3583041](https://doi.org/10.1109/TMECH.2025.3583041)，pp.7698–7701，式11–17及实验参数；T-RO与Moreira细节见 `2026-11-04_codex_contact_physics_audit/final.md` | 已有本地全文核读；仿体参数不能冒充人体耦合剂数据；Moreira [10.1016/j.cmpb.2013.11.005](https://doi.org/10.1016/j.cmpb.2013.11.005)为离体松弛 |
| E3 | [Remus全文](https://www.frontiersin.org/journals/bioengineering-and-biotechnology/articles/10.3389/fbioe.2024.1384062/full)，§2.3压头装置、§3.1/表4、§3.2/表5、§4 | r10 mm、5±1 mm/s；表5八个实测点用于双接触律回代，q/E*是迁移推算，未把Ogden参数误认杨氏模量 |
| E4 | [10.3233/WOR-172611](https://doi.org/10.3233/WOR-172611)，摘要Methods/Results | 沿用已查摘要；独立核对未能打开。数值是小样本分布，不是安全范围 |
| E5 | [Courrèges](https://link.springer.com/chapter/10.1007/978-3-319-22368-1_21)，Abstract及subscription preview | 本轮可见摘要，无腹部HC数值；不下载付费全文 |
| E6 | [人工扫描速度原文](https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2022.926006/full)，扫描速度结果段；具体数值沿用已有audit §B.5 | 扫描实验数据，非速度/加速度临床标准 |
| E7 | [10.1109/TIE.2020.3036215](https://doi.org/10.1109/TIE.2020.3036215)，摘要；[SAPM](https://pmc.ncbi.nlm.nih.gov/articles/PMC7614654/)，旋转调整能力结果 | 已有audit §B.3；平均绝对角差并非有符号误差分布 |
| E8 | [UR5e规格](https://www.universal-robots.com/manuals/EN/HTML/SW10_8/Content/prod-usr-man/complianceUR5e/H_g5_sections/appendix_g5/tech_spec_sheet.htm)，Pose repeatability；KUKA信息见已有final_web/verify | 重复精度、可设定阻抗不等于位置接口实测刚度 |
| E9 | [FT300手册](https://assets.robotiq.com/website-assets/support_documents/document/online/FT_Sensor_Instruction_Manual_Web_20181218.zip/FT_Sensor_Instruction_Manual_Web/Content/Specifications.htm)，§5.3；[Mini45规格](https://www.ati-ia.com/products/ft/ft_models.aspx?id=Mini45)，Ranges and Resolutions | 已有final_web/verify；FT300噪声与输出率已核实，Mini45独立核对超时；不得将分辨率当σ |
| E10 | [Virga机构摘要](https://portal.fis.tum.de/en/publications/use-the-force-deformation-correction-in-robotic-3d-ultrasound/)，Results | 5 N、55 mm成像深度、3.39±1.86 mm校正误差，不提供所需压入深度 |
| N1 | [Boyer 2007原始摘要](https://pubmed.ncbi.nlm.nih.gov/18003026/)，Abstract；[10.1109/IEMBS.2007.4353360](https://doi.org/10.1109/IEMBS.2007.4353360) | 4人，10–60 Hz、1–10 μm；k47.3–128.3 N/m、c.08–.121 N·s/m，复模量13.2–33.4 kPa。只见摘要，无逐频配对相位曲线。以未配对极值包络推算20.8Hz KV相位4.7–18.5°，非原文实测相位，更不是腹壁相位 |
| N2 | [Schimmoeller全文转录](https://www.researchgate.net/publication/329203670_Instrumentation_of_Off-The-Shelf_Ultrasound_System_for_Measurement_of_Probe_Forces_During_Freehand_Imaging)，§3.2–3.3、Fig.8–9、pp.122–123；[10.1016/j.jbiomech.2018.11.032](https://doi.org/10.1016/j.jbiomech.2018.11.032) | 可读原文转录；在体10.74 N/18%压缩与尸体7 mm分开报告；未查到同一腹部探头的可靠N—mm配对值 |
| N3 | [Patterson全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC7185178/)，Materials and methods/Tissue sample morphology、Fig.1、Table1 | 人体尸体、非活体；Table1 B行32 mm与汇总最大31 mm不一致，本文按逐行范围明确标出。rib spacing为文中影像肋间距离，不当作肋中心节距；精确边缘定义须随Fig.1使用 |
| N4 | [Quirk全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC5714366/)，§II及Table2、Fig.2 | 分析25名女性志愿者腹部AP体表位移，8±3 mm为峰峰、4.2±1.1 s；不是内脏SI运动；由周期幅值推算速度，不声称原文测得峰速 |
| N5（未采用数值） | [von Gierke 1952](https://doi.org/10.1152/jappl.1952.4.12.886) | 原文不可访问，2.5 kPa/15 Pa·s组合本轮**未查到原始证据**，不用于标定 |

N5未取得可用数值，不用于标定；其余条目保留已有原文核读及访问边界。对象实现未新增文献检索。
