# ARIS 第一轮全文阅读（2026-10-31）

## A. Controller 状态与准入

- run_id：`recent-related-work-20261031`；RUN_ROOT：`D:/pip_env/temp/scan_impedance_literature_20261031/aris_run`。结束仍为 **PAPER_READING**；本轮 **7 张 Evidence Card 全部 PASS 并接受**，累计由 **5 增至 12 张**。未调用 finish-reading，未提交/进入 Field Map。
- 7 篇均经 `read-user-fulltext`、Controller 全文交付、独立原生 paper_reader 及主代理核查；只审读指定关键部分。卡片绑定原 PDF 哈希、read_event_id 和页码。卡片草稿/绑定见 RUN_ROOT 的 `read1_work/`；正式卡片由 Controller 写入 `.aris/canonical/recent-related-work-20261031/`。
- **4 篇新增来源未能正式登记/准入，不能记成已接受。** 已复制至 `source-materials/user_{Wei,Lei,Zakeri,Chen}_2025.pdf`，保留原文件；登记请求保存在 `read1_work/register_user_*_2025.json`。政策允许用户提供材料走 `USER_SUPPLIED_READ`，不受 venue/引用门槛限制；但当前 347 条候选中没有这 4 篇，PAPER_READING 的 allowed-actions 只有 `promote_user_source`，没有新增登记；现有 `register-user-source` 实现限定 METADATA_RETRIEVAL，而 `promote-user-source` 必须已有当前候选。因此未执行越阶段登记、未改状态，也未读取未准入全文。**这是 Controller 登记入口限制，不是相关性/政策排除。**

|用户来源|正式准入结果|本轮处理|
|---|---|---|
|Wei，TMECH 2025，10.1109/TMECH.2024.3443310|待登记；无正式 paper_id|复制及准备请求；未读，等 Claude Code 笔记|
|Lei，IEEE Sensors J 2025，10.1109/JSEN.2025.3598958|待登记；无正式 paper_id|复制及准备请求；未读，等 Claude Code 笔记|
|Zakeri，TIE 2025，10.1109/TIE.2024.3451138|待登记；无正式 paper_id|复制及准备请求；未读、无卡|
|Chen，TII 2025，10.1109/TII.2024.3477549|待登记；无正式 paper_id|复制及准备请求；未读、无卡|

四份请求均保留指定理由：**“用户在 Claude Code 中提供，为上轮临时检索发现、本轮正式检索未纳入的疑似遗漏”**。后续需先通过 Controller 合法入口接纳这 4 篇，再读 Zakeri/Chen；不应把恢复登记误当作要求用户重交全文。

未读 Wei、Lei、Xiao `memfv6bMNssJ`、Beber `ODwWwSVciuYJ`、Zhang Sci Rep `SS4c7arItfoJ`、Ebrahimi `HokaNUBl3lAJ`；其他未指定论文保留下一轮。未修改项目或 ARIS 源码，未执行 Git pull/commit/push。

## B. 本轮论文（页码均为 PDF 页序，从 1 起）

### 1. Qian 等，ICRA 2026｜P0yee_3nRaAJ

题名：*Breaking the Latency Barrier: Synergistic Perception and Control for High-Frequency 3D Ultrasound Servoing*；DOI：10.1109/ICRA57385.2026.11697234。
- **问题**：快速、大范围目标运动时，超声感知和动作推理延迟使参考切面跟踪掉队。
- **机制/为何有效**：双流分别估计平面内几何位移和平面外语义位移；A-UKF估位置/速度/偏置；单步 Flow 输出短时速度序列，避免迭代采样延迟。8帧输入、8步预测、执行前4步再规划（pp3–4）。
- **实验证据**：UR3e/CIRS体模、RTX4080，16.2 ms、约62 Hz；102.47 mm/s 时**全程平均误差6.124±0.386 mm，终端1.629 mm**；173.2 mm重定位终端1.5219 mm（pp5–6，Tables I–II）。不能写成动态误差1.5 mm。
- **依赖**：实时二维US及操作员选定参考图；输出三维期望平移速度。20000模拟图像对、1000模拟轨迹、50条实体专家轨迹；NDI用于真值，力传感器/力与速度限值用于接触安全（pp4–5、7）。
- **失效/原因**：DCL-Net替代感知虽73 Hz仍不收敛；Z轴旋转20°不稳、25°失败，作者归因剪切与平移歧义；X/Y旋转超出范围（pp6–7）。单志愿者结果是演示，非临床泛化。
- **关系/路线意义**：**图像运动感知与高速预测控制路线代表**。证明速度与观测可靠性需一起解决；不直接验证本文力/刚度换算观测及其增益耦合。A-UKF各组件作用主要为作者解释，不能当独立消融结论。

### 2. Jiang 等，TIM 2022｜vHYta8tPiJ8J

题名：*Precise Repositioning of Robotic Ultrasound: Improving Registration-Based Motion Compensation Using Ultrasound Confidence Optimization*；DOI：10.1109/TIM.2022.3200360。
- **问题**：患者改变位置使扫描轨迹失配、三维血管体积拼接中断。
- **机制/为何有效**：RGB分割掩膜Dice<0.95触发停扫；运动前后表面ICP更新轨迹、断点续扫；US置信图修正接触方向，再用机器人位姿及血管质心对齐改善拼接（pp4–6、9–10）。
- **实验证据**：手臂体模，无遮挡/机械臂遮挡各20次，重定位误差6.1±1.6/6.9±1.6 mm；运动范围110 mm、最大旋转80°；平均ICP迭代336 ms（pp8–9）。p=.14不等于两条件等效。
- **依赖**：CT血管/表面模板、Azure Kinect RGB-D及手眼标定；iiwa内置关节力矩、力矩柔顺执行，中心线恒力其余5DoF跟踪。柔顺刚度125–500 N/m按组织选用，非本轮证实的在线辨识（pp3–4）。
- **失效/原因**：ICP敏感初始化且仍有拼接间隙；作者限于可见刚性运动，震颤/血管运动未处理；称柔顺可部分适应小快运动，但无对应定量证明。精调改善可视化，作者明确不能提高三维几何准确性（pp9–10）。
- **关系/路线意义**：**外部几何感知、事件式重定位路线代表**；与本文连续恒力预测补偿任务相邻，但不能将“停扫—恢复”说成连续呼吸前馈能力。

### 3. Wang 等，TASE 2024｜f76lmG78c6cJ

题名：*Task Space Compliant Control and Six-Dimensional Force Regulation Toward Automated Robotic Ultrasound Imaging*；DOI：10.1109/TASE.2023.3282974。
- **问题**：位置接口机械臂如何同时实现6DoF柔顺、6D力/矩调节，并避免受限接触中的实现性失稳。
- **机制/为何有效**：体坐标几何导纳；接触时每周期用实测关节运动重置导纳误差状态，避免理论积分与实际伺服运动持续分离；自由空间保留理想积分以保轨迹精度，1 N阈值切换（pp5–7）。
- **实验证据**：UR5e低层滞后约20 ms，而控制周期2 ms；同参数指定仿真/手持轻摇中ideal失稳、practical/hybrid稳定。自由空间均路径误差<0.30 mm；3 N/.2 Nm调节收敛后平均相对力/矩误差<0.5%/<0.8%（pp5–9；Table IV目视核对）。
- **依赖**：UR5e位置/RTDE接口，ATI六轴F/T、编码器速度与数值加速度、负载重力/惯性补偿；2 ms外回路。人体脊柱/肾脏/颈动脉仅定性成像（pp7、9）。
- **失效/原因**：作者将ideal失稳归因实际跟踪滞后造成状态累积偏差；practical在空域精度较差。稳定性限所测接触条件，不是通用证明；低力残差接近传感器分辨能力，也不是运动组织动态误差（pp5–9）。
- **关系/路线意义**：**位置接口导纳执行的最近邻之一**；说明底层滞后及内部/实际状态一致性已有直接机制研究，但未在所阅核心中解决本文运动/刚度分离与力导出运动预测。

### 4. Zhang 等，RA-L 2025｜CkE-5mElf_IJ

题名：*Smooth Path Planning and Dynamic Contact Force Regulation for Robotic Ultrasound Scanning*；DOI：10.1109/LRA.2025.3604746。
- **问题**：乳房扫描路径不光滑，以及恒定较大力引起变形、较小力又不能充分贴合。
- **机制/为何有效**：路径双目标平滑；预触诊拟合三次力—应变，结合厚度/探头完整贴合深度预测下一点力；相邻点间用三分区US置信度更新目标力，底层有界力—法向速度映射（pp2、4–5）。
- **实验证据**：三静态硅胶/水凝胶体模、四路径、1 mm/s；变力以较低力获得接近恒力8/6/3 N的置信度均值/SD，未核得精确降力百分比。路径平均绝对偏差相对NFO/EBTSO降低21.3%/6.7%（pp6–8）。
- **依赖**：UR5速度命令、F/T、D435i点云、US置信图；预扫三位置触诊，假定相邻局部变形比γ=d/H连续；机器人125 Hz、力1 kHz（pp4–7）。
- **边界/疑点**：仅三体模，压缩导致路径变化待有限元在线修正。**原PDF p5式(11)–(12)已目视确认：ΔC≥0，正目标力更新因子≥1，仅能非减；与文字调力意图存在歧义**，不自行修符号，也不由此否定整体实验。
- **关系/路线意义**：**图像质量驱动力参考适配代表**。“动态”是期望力沿路径更新；不能推成动态组织辨识/呼吸预测，也不是在线刚度归一化阻抗。

### 5. Ning 等，TIE 2024（在线2023）｜YriQFb_sHqEJ

题名：*Inverse-Reinforcement-Learning-Based Robotic Ultrasound Active Compliance Control in Uncertain Environments*；DOI：10.1109/TIE.2023.3250767；首次在线2023-03-06。
- **问题**：未知形状、刚度及接触变形下仅用力反馈适应探头姿态和扫描接触。
- **机制/为何有效**：5帧六维力历史输入；学习期望法向力、roll/pitch矩及对应速度，再经**固定M/B/K**导纳变位置命令。先PPO，再用人工纠正的500回合示范做AIRL；纠正提供额外姿态/速度信息（pp3–7）。
- **实验证据**：静态姿态误差2.8°，未纠正IRL/RL为6.2°/6.1°；未见体模六路径各10次，相对静态模型姿态差2.3±1.3°/1.9±1.2°；平均压力9.0±3.6 N是变力统计，**非恒力误差**（pp7–9）。
- **依赖**：UR3位置接口、FT300六维F/T，动作10 Hz、转速≤0.1 rad/s；实体体模训练及人工修正；立体相机用于评价参考模型，US用于成像，非在线视觉控制（pp6–9）。
- **失效/原因**：作者明确快速目标运动导致力突变、可能误决策；仅静止/轻微呼吸，单名在体演示。静态3D法向及人工扫描是参考而非变形真值；方法明确单步反应而不预测未来（pp3、7、10）。
- **关系/路线意义**：**仅力反馈学习接触动作的最近邻之一**；已有未知表面适应能力，不能说仅力方法都依赖已知表面；但它并非学习阻抗矩阵，也未证明本文的运动/刚度分离及预测能力。

### 6. Wu 等，T-RO 2026｜wftPH5bwT0gJ

题名：*A Lightweight Robotic Ultrasound System With Hybrid Force Control and Modular Force-Range Switching*；DOI：10.1109/TRO.2026.3734315；正式在线2026-09-15。
- **问题**：呼吸低频变形与高频震颤同时存在时的轻量、可换力域恒力超声。
- **机制/为何有效**：被动恒扭矩机构缓冲快速扰动，主动PI维持机构变形补偿低频运动；重力/迟滞/摩擦补偿，换绕线轮改力域。**刚度匹配是机构刚度，不是组织刚度辨识**（pp2、6、8）。
- **实验证据**：0.5 Hz、3–7 mm叠加15 Hz台架扰动；角跟踪RMS由1.67°降至0.59°。力落在均值±5%/±10%区间占比作者报告提高34.98%/15.94%，不能当参考力RMSE；同一志愿者深呼吸图像质心平均振荡0.10 mm、灰度SD0.94（pp9–11、16）。
- **依赖**：610 g专用恒扭矩末端、三编码器、倾角计、伺服及补偿标定；量测机构形变而非估组织运动；可手持或接UR5遥操作。40 mm轮双补偿有效力域4.07–5.41 N（pp6、8、11）。
- **失效/原因**：恒扭矩窗口仅4.3°，电机速度约106 r/min，过快扫描失联；摩擦力峰、探头变重缩力域，残余高频振荡可模糊深部并恶化某些尺寸指标（pp8、10、13、16）。
- **关系/路线意义**：**机械顺应与主动补偿的直接相关代表**；改变扰动通道也能实现收益，不必先估组织运动。专用硬件不可作为排除理由；0.10 mm是图像指标，不是运动恢复或力跟踪精度。

### 7. Piedra 等，IJRR 2026｜BIP_InW4FSIJ

题名：*Explicit force-based control strategies for compliant robot-assisted ultrasound imaging*；DOI：10.1177/02783649261441635。
- **问题**：隐式导纳经过运动内环带来滞后，移动及软硬混合组织下需同时稳力和优化探头贴合。
- **机制/为何有效**：操作空间动力学及正交投影，轴向显式力PI+目标力前馈、两轴零矩、速度阻尼、动力学补偿；直接给关节力矩，绕过中间位置执行链。粒子滤波用于接触/约束方向，不是组织运动IMM（pp6–11）。
- **实验证据**：0–20 N、1.6 Hz目标下相滞5.76°，重调IAAC121°，均差1.01±0.61 N；13.4 N恒力、体模15 cm往返约5 cm/s，向/离均差1.25±0.37/1.43±0.45 N，重调IAAC5.61±1.21/5.22±1.71 N。混合软硬面3 cm/s扫查均差0.33 N（pp16–19，Tables2–4）。
- **依赖**：Panda力矩接口、动力学模型、ATI腕F/T及负载补偿；系统有D455世界相机/AprilTag追踪体模，逐实验视觉作用未完全隔离，不称全套仅力方案。100 Hz/5 Hz是机器人/US**记录率**；13.4 N来自专家记录（pp6、9、12–14）。
- **边界/原因**：全为体模；1.6 Hz峰误差仍3.37 N，法向/姿态扰动实验均差1.31 N。诊断最优力未解决；作者明确置信度只是声耦合。回放置信0.32对手动0.27，同时改变目标力和矩控制，不作单因素归因（pp16、20、24）。
- **关系/路线意义**：**动态恒力及异质组织控制的直接近邻**。执行接口本身可显著改善扰动抑制；本文位置接口与换算增益耦合问题须独立论证，不能借接口差异略过其已实现能力。

## C. 被引用但本轮候选库未收录的必要线索（未新增检索/准入）

以下已按题名片段检查当前候选库；仅为被引线索，不代表身份已核验，更不是“无人做过”的证据。部分可能已存在于此前项目背景记录。

1. **Force tracking control method for robotic ultrasound scanning system under soft uncertain environment** — Jiang J、Luo J、Wang H 等，Actuators 2024, 13:62，10.3390/act13020062；Piedra p26。其所对比IAAC的原始来源，需核清自适应组织刚度、接口及移动环境假设。
2. **Towards autonomous atlas-based ultrasound acquisitions in presence of articulated motion** — Jiang 等，RA-L 2022, 7(3):7423–7430；Zhang p8 ref9。补上刚性重定位之外的关节运动路线。
3. **Autonomous trajectory planning for ultrasound-guided real-time tracking of suspicious breast tumor targets** — Tan 等，TASE 21(3):2478–2493；Zhang p8 ref14印为2024卷期，Qian p8 ref12印为2023，在线/卷期身份待核。与Chen邻近的US目标实时跟踪路线，不能遗漏后只概括静态扫描。
4. **Full-coverage path planning and stable interaction control for automated robotic breast ultrasound scanning** — Wang Z 等，TIE 2023, 70(7):7051–7061；Zhang p8 ref11。Zhang力—速度控制/路径组织的直接前序，需恢复其完整接触能力。
5. **Deformation-aware robotic 3D ultrasound** — Jiang 等，RA-L 2021, 6(4):7675–7682；Jiang TIM p10 ref16。属于2022年前必要机制背景，不计2022–2026新发现；用于厘清按刚度修正压缩形变与本文估计运动的关系。

Qian 引用的 *Robust deep feature ultrasound image-based visual servoing: focus on cardiac examination* 已有候选 `2H7kUPP7rC8J`，因此未列为遗漏。Zakeri/Chen已在A列待登记，不在此重复计数。
