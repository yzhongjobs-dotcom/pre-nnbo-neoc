# 超声恒力扫描：场景设置与对象边界

2026-11-05。设置实现在 `experiments/rossi_failure_study/physics_plant3.py`、`scenes3.py`；数值验收见 `results/plant3_checks.json`。旧对象及实验保留。不接 Method，不运行闭环或模型选型。以下推算、迁移近似及假设均不等同人体材料辨识；来源与原文访问边界保留于文末。

## 对象验收后的最终取法（优先于下表）

- **证据不足的量按范围逐场景抽样，不取单一代表值：** 总相位和生成机制分层（见下）；肋上覆盖 7、18、32 mm；接触律两种；探头、力档、滑动方向分层（在每种组织和工况内均衡）。μ .03–.30、k_d 0–5、K_t 5–50 kN/m、kr 50–200 kN/m、ρ ±.3、加速度 .05–.5 m/s²、呼吸 3–12 mm 和 3–6 s 连续抽样。用 overrides 可固定任一项，做敏感性检查。
- **找回的场景机制：** 惊跳和咳嗽时肌肉收缩使组织变硬（增益 2，由 Remus 主动与放松之比约 3 推算）；体表起伏 .3–2 mm、相关长度 3–6 cm（扫查由此产生法向运动）；心动微运动；姿态误差按在体 2.9±1.6° 取绝对值并随机取号。
- **肋骨不是刚性的：** 覆盖层硬化到局部刚度等于肋骨和胸廓刚度（10–45 N/mm 分摊到探头面积；整体胸廓约 9 N/mm、单根肋约 45 N/mm，均只见检索摘要，未核原文）后线性延伸。接触律因此没有奇点，探头压在骨上的刚度受机械臂和肋骨串联限制。
- **动态接触不只由 HC + Prony 生成：** 每个场景的总相位取 7、11、17、23°，由三种机制之一实现：HC 加标称记忆、三支路弱幂律（无速率项）、单支路 Kelvin-Boltzmann（无速率项）。三种分层均衡，选型结果不依赖某一候选的形式。
- **评价力：** Method 控制和测量的是轴向力，J 用轴向力计算；法向力单独报告。两者只差姿态和摩擦的投影 Fn = (Fa − Ft sinθ)/cosθ，主要在无姿态跟踪时显著（中位 RMSE .318 对 .113 N）。
- **姿态符号：** 表面前方下降时 θ > 0：压力后移，摩擦使轴向力增大，法向力把探头推向下坡。旧对象（plant2）这里的符号相反。
- **半空间接触律用闭式：** F = (1.1/√A)·(1.2·F_E) 增长期，占满后 (1.1/√A)·(F_E + .2·F_E(占满))。对齐均匀时与 ∫k 精确一致，凸阵加覆盖层时差 13%（切线 18%），属迁移假设的一部分；真值列用力的实际斜率。
- **肋骨边缘：** 覆盖层把骨头对表面的影响展开，边缘宽度取覆盖层厚度的 .5–1 倍（物理推理，非实测）。**混合组织：** 夹杂模量比 2–8（假设，参考仿体标称约 4 倍）；骨头只由肋骨工况表示。
- **已接受的例外：** 平底对齐初触时切线有阶跃（几何必然，力连续）；4 mm 覆盖端点回放较慢（约 4.9 s/仿真秒），只用于敏感性子集。

## 设置表

| 机制 | 设置与依据边界 |
|---|---|
| 探头 | 理想化凸阵长69、宽17、R60 mm；R43 mm为几何端点。线阵平底48.43×11.56 mm。尺寸依据E1，凸阵半径组合是假设，不冒称某款完整产品。抛物面近似与R60 mm圆弧在边缘高度约差10%。AuSoScan 80 mm与R30 mm不能当同一圆弧弦长和半径。 |
| 放松组织 | q=1.0–1.8 MN/m³、E*=15–31 kPa，按同一抽样分位对应；q=1.1来自R3（脐下腹直肌），q=1.6来自R1（上腹中线），均不是侧腹。对应E*约17/28 kPa，均为Remus压头数据迁移推算，非超声探头实测模量［E3］。 |
| 主动、混合组织 | 主动q=3.4–4.6 MN/m³、E*=54–72 kPa；q代表4、E*代表63。由主动满半球斜率1.08–1.43 N/mm分别除以πr²及2r推算。混合为放松软层中5–30倍硬夹杂的假设；均保留相关空间场。 |
| 力与方向 | 场景字段FD为4/8/12 N；8 N为中心档，非安全阈值［E4］。探头、力档、长/短轴滑动由(cell,i)均衡分层，不展开三者全组合矩阵。长轴滑动沿长轴平均空间场，短轴滑动沿宽度平均，法向刚度始终取完整足迹。 |
| 接触迁移 | 每个场景以contact_law选择foundation或halfspace；二者均为待验证的对象近似，不能从其中之一推出通用选型结论。公式见下。 |
| 肋与覆盖层 | 肋上覆盖18 mm代表，7/32 mm端点［N3，尸体数据］；另设4 mm薄层端点，代表胸骨/脊柱类浅骨的假设。肋带宽10–20 mm、平滑边缘1–4 mm是假设；肋间空隙6–19 mm，节距=本条肋宽+空隙。只在肋带局部增硬。 |
| HC | 20.8 Hz、8 N平衡点设单独HC相位0–30°，代表5°；45°只作应力检查。α=.35三支记忆下，裸接触总相位分别为6.4–22.8°、代表9.0°、应力32.8°；不能直接相加相位。λ=tan(φHC)k/(ω·8)，随后固定，不随FD或切向速度重设。腹部HC强度未辨识，以上是假设［N1/E5］。 |
| 记忆 | 因果Prony α=[.35,.35,.35]、τ=[.008,.08,.8] s；α整体×.5/1/1.5、τ×.3/1/3供条件检查，不为速率工况消除记忆［E2，记忆存在有依据，数值是假设］。 |
| 扫描轨迹 | 常规5–30、慢扫1–5、快速压力组40–60 mm/s；加速度代表.1、范围.05–.5 m/s²。启停暂停.2–.8 s；换向用五次平滑速度过渡。速度依据E6，暂停、加速度是假设。足迹及轨迹留在空间场内。 |
| 曲面与姿态 | R_s=.15–.5 m是假设；同一表面有/无姿态跟踪配对。跟踪残差±3°代表、±5°端点；E7的平均绝对角差2.9±1.6°不能当有符号高斯标准差。 |
| 呼吸 | 峰峰8 mm、周期4.2 s代表；3–12 mm、3–6 s通过覆盖参数检查［N4］。平滑cos²周期，代表峰速约6 mm/s为公式推算。当前代表形状吸呼对称，不冒称已拟合呼气不对称。 |
| 推拉、惊跳、咳嗽 | 慢推拉2–15 mm/.5–1 s，快推拉2–12 mm/.15–.4 s，五次位移过渡；峰速=1.875A/T。惊跳1–5、咳嗽.5–3 mm用起止零斜率sin²脉冲，持续时间πA/v_peak。常规总驱动峰速≤50 mm/s，strong组100–150 mm/s。以上为假设，快推拉属较强扰动。默认仅支撑位移；external单列外力驱动，取消同时叠加位移、外力、反射增硬。保留rest/breath/shift_slow/shift_fast/flinch/cough/tremor/rib/combo。 |
| 摩擦 | μk代表.10，范围.03–.30、端点.45；μs/μk=1–1.3，v_s=2–10 mm/s，平滑速度v_e=.5 mm/s；k_d=0–5、13 N·s/m端点。k1=.001、k2=.003 s仅取论文参考，不随机放大。AuSoScan式压力不对称默认开，可关闭；半径有探头边界，详见下述近似边界。人体耦合剂参数未测［E2］。 |
| 机器人与支撑 | kr代表100、范围50–200 kN/m；K_t代表20、范围5–50 kN/m。Cnn=1/kr，Ctt=1/K_t，Cnt=ρ/√(krK_t)，ρ=0代表、±.3端点，保证互易与正定。身体m=15–40 kg、7–10 Hz、ζ=.25–.5仍为支撑假设；重复精度不是刚度实测［E8］。 |
| 测量 | 默认沿用1 kHz、σ=.08 N、偏置游走2e−4 N/√tick及旧30 Hz一阶低通；可选sensor_lp=100 Hz、σ=.1 N。1 kHz及100 Hz低通是仿真链假设，不代表FT300输出率；FT300手册输出100 Hz［E9］。法向力、切向力、力矩独立输出；不将分辨率当σ。 |

## 接触几何与迁移公式

统一局部压入：

\[
\delta(\xi)=d-g(\xi)-\xi\tan\theta-\frac{\xi^2}{2R_s},\qquad
g(\xi)=0\;(\text{线阵}),\quad g(\xi)=\frac{\xi^2}{2R}\;(\text{凸阵}).
\]

θ为探头相对局部表面法向的倾角；部分接触由δ>0和有限足迹共同确定。平面R_s=∞。宏观表面高度、局部法向与曲率取同一抛物曲面；无跟踪时探头保持固定全局姿态。压力中心及法向力矩来自同一分布积分，力矩约定M_y=F_n·ξ_cp，工具切向力另加杠杆项。

foundation采用局部压力p=qδ/(1−δ/h_local)，无肋h_local=∞。肋带内部h_local=h；边缘用余弦窗w(x)平滑逆覆盖深度1/h_local=w/h，窗外为零。该肋边缘形状是假设，避免在无肋区域施加整片增硬。均匀无背衬凸阵：a=min(√(2Rd),L/2)，F=qW(2ad−a³/(3R))、k=2qWa；线阵F=qWLd。占满后两者继续增力。覆盖层必须先作用于每个局部δ再积分，不能用旧的整片F₀/(1−d/h)替代凸阵局部积分。

halfspace采用矩形半空间近似，β=1.1，均匀无背衬k=βE*√(W·2a(d))，F=∫k(s)ds。线阵F∝d，未占满凸阵F∝d^1.25；占满后F继续增长，F与k连续。空间场及肋上局部切线因子(1−δ/h_local)^−2作用于面积平均权重。压力中心由各深度增量载荷的一阶矩积分得到；这是对标量刚度规律的一个明确分布假设，不等同精确半空间压力解。

半球回代使用r=10 mm有限半球。foundation为q乘实际压入体积；halfspace使用圆形接触k=2E*a(d)，a=√(2rd−d²)至d=r，此后a=r。E*初值则按任务指定Hertz式F=(4/3)E*√r·d^1.5反推；Hertz不能无条件外推到d≥r。表5全部八点对两种迁移律的预测和误差保存在JSON，不通过调参消除偏差。常q在R1的5 mm点误差−27.8%，有限半球halfspace的最大误差约32.1%；均不是完整材料模型验证。若另假定变形深度H=20 mm，qH约22/32 kPa也只是等效换算，不是辨识出的杨氏模量。

## 动态、观测与验收边界

F=(f_e+Σq_i)(1+λḋ)，Prony状态按每步弹性力变化因果更新。平衡小信号为

\[
K_c^*(\omega)=k\left[1+\sum_i\alpha_i\frac{j\omega\tau_i}{1+j\omega\tau_i}\right]+j\omega\lambda F,
\qquad K_{series}^*=\frac{K_c^*}{1+K_c^*/k_r}\quad(\theta=\rho=0).
\]

日志0–19保留原列含义，E字段只给E*，不塞入q；追加20–25为Re(K*)、|K*|、相位rad、接触长、长轴压力中心、最小局部覆盖深度。复响应真值是冻结当前位置/速度、包含机器人柔顺与力投影的局部条件响应，不含身体支撑、空间移动、执行器及测量链传递函数。Obs仅公开测量及编码器/指令位置，真值列全部屏蔽。episode/replay入口保留，未接任何Method。

默认记录周期8 ms；开环幅相验收保持8 ms阶梯位置指令，额外以log_dt=.5 ms记录并分离阶梯输入的频谱镜像，避免125 Hz日志混叠改变20.8 Hz幅相。数值积分2 kHz，测量仍1 kHz。半空间在未触及边缘及占满区使用等价解析原函数，中间单边接触区保留深度求积；短轴相同覆盖深度的权重合并。计算时间只报告预热后的回放，不作实时认证；普通场景通过三倍限制，倾斜跨肋场景仍超限，具体比值见JSON。

理想平底线阵在平面上首次整面接触时，指定线性规律必有右侧刚度qWL或βE*√(WL)，脱离侧为零；力连续而切线不连续。它与“任何首次接触均切线连续”的要求数学上不相容，验收对此明确标失败；未暗加接触平滑层。凸阵占满、跨肋及其余脱离边界独立核验。

压力不对称采用有界前后权重重分配，保持δ确定的接触域，关闭时恢复对称积分；不把截断后AuSoScan的扩展半径原式当作真实超声探头压力解。其方向、幅值及半空间增量压力分配均需实物标定。μ=0时关闭该摩擦诱发的压力偏置。

## 复现与敏感性入口

`python -X utf8 experiments/rossi_failure_study/check_plant3.py`只执行静态、开环和场景验收。静态表涵盖双探头、双接触律、R3/R1/主动及4/7/18/32 mm整片背衬、4/8/12 N；局部跨肋另由密集积分核验。Remus重建只核对量级，不用于选型。

`cell_scene3(cond,tissue,regime,i,base,T,overrides)`确定可复现的场景。组织relaxed/active/rib/mixed与九个条件组成CELLS3；nominal/slow/stop_go/reversal/fast/curved/curved_untracked/strong/external分别保留。曲面配对仅改变姿态跟踪。

overrides支持contact_law、probe、force、slide、phi_hc、alpha_scale、tau_scale、q_scale、cover、rib_gap、mu、K_t、rho、motion_peak、sensor_lp；另可覆盖geometry、kr、pose_error、surface_radius、breath_pp、breath_period、accel等。代表值不自动展开为全范围矩阵；未知量条件检查仍应围绕几何迁移、动态响应、局部背衬、切向耦合和运动驱动五组组织。

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
