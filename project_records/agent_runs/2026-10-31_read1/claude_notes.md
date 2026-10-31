# Claude Code 精读笔记（ARIS run recent-related-work-20261031，第一轮）

按 skill 阶段一"五个问题"读到方法与实验；页码为 PDF 页。供 Codex 转成 Evidence Card，并用于改写研究现状。

## Wei et al., "Composite Disturbance Filtering for Interaction Force Estimation With Online Environmental Stiffness Exploration," IEEE/ASME TMECH 30(3), 2025（在线 2024-08），10.1109/TMECH.2024.3443310（用户来源）
- 问题：微创手术中无力传感器时估交互力，同时辨识环境刚度。
- 机制：把交互力看作机器人状态耦合扰动；广义动量动力学 + 线性各向同性力生成模型 F=Kδ 组成耦合系统，给出交互力可分离条件；EM（E 步滤波状态、M 步更新刚度）构成 EEFO（p.1984–1986，式(7)(9)–(13)）。
- 做到：仿真中常/变刚度辨识约 1.6 s 内收敛（p.1987 Table II）；Toumai 手术机器人+硅胶肾实验，对比 KF/GMO/DO/NDO/DKF 误差最小，恒力段稳定时间缩短 >63.4%（p.1989）。残余动力学用离线 GP 预训练补偿。
- 依赖：关节力矩与动力学模型（非位置接口）；**明确假设静止对象**："Considering only interactions with static objects, it follows that r = 0"（p.1985，式(9)下）——形变只由机器人运动产生。
- 停在哪/为什么：只估计不控制；作者把"环境动态特性（如跳动心脏）"和"用于动态交互中的阻抗控制"列为未来工作（p.1990 结论）。
- 与本文：路线一（激励/模型辨识刚度）2024–2025 最新代表；直接支持缺口——联合刚度辨识的最新顶刊工作以静止环境为前提，且未进入控制回路。

## Lei et al., "Robotic Respiratory Liver Motion Compensation for Stable Ultrasound Imaging of Small Lesions," IEEE Sensors J. 25(19), 2025, 10.1109/JSEN.2025.3598958（用户来源）
- 问题：呼吸使肝脏病灶偏离超声成像平面，需稳定成像以便穿刺。
- 机制：ReME-Net 以当前帧与目标帧超声图像差分特征 + 通道-空间注意力，端到端估计 6 自由度肝脏运动；机器人据此调整探头位姿；法向用按力误差在线调参数的自适应变导纳维持接触（式(5)附近，p.37002），不估计组织刚度。
- 做到：CIRS 腹部仿体装在六足台上做周期运动；平移幅值 19.72 mm、旋转 17.32° 时，平移/旋转误差 3.05 mm、2.52°，平均接触力误差 0.25 N（结论）。作者注意到补偿存在"轻微时延，主要来自补偿策略和执行器响应"。
- 依赖：实时超声图像与离线训练数据；周期运动；均质仿体。作者自述局限：仅仿体，真实人体伪影（声影等）会使运动估计不完整。
- 与本文：路线三（测量运动再补偿）2025 代表，运动来自超声图像而非受控探头的力；无法做仿真对比，但必须写入研究现状。

## Xiao et al., "Optimizing Variable Admittance Control for Remote Ultrasound Scanning Under Uncertain Environment," IEEE Access, 2025, 10.1109/ACCESS.2025.3568150（memfv6bMNssJ）
- 问题：位置控制的远程超声机器人在未知环境中精确力跟踪。
- 机制：线性弹簧 f=k_e(x_e−x)；融合力与超声置信度估计接触位置 x_e（接触时刻记录一次），另一线程在线估 k_e；由粗到细：先调参考位置、再调虚拟刚度并加误差补偿项；能量罐保无源（Sec. II–III）。
- 做到：UR3e + ATI Nano25；仿体与人上肢，稳定段力误差约 0.2 N 以内（Fig. 7, Table 2）。
- 依赖：位置接口（与本文相同）；环境在接触后视为静止（x_e 固定）；无运动补偿。
- 与本文：位置接口下"在线刚度+接触位置估计"的超声近邻，出处层级一般；可作路线一在超声中的近年实例。线索：其引用 [19] Li et al., "Stable and compliant motion of physical human–robot interaction coupled with a moving environment using variable admittance and adaptive control," IEEE RA-L（移动环境+变导纳，需核年份与机制）。

## Beber et al., "A Passive Variable Impedance Control Strategy with Viscoelastic Parameters Estimation of Soft Tissues for Safe Ultrasonography," ICRA 2024, 10.1109/ICRA57147.2024.10610167（ODwWwSVciuYJ）
- 问题：远程超声中安全的变阻抗控制。
- 机制：非线性黏弹（Hunt–Crossley 型）参数**离线**按网格触诊用最小二乘估计，GPR 建成弹性/黏性图；二次规划在线调刚度，加入由黏弹估计导出的力约束与能量罐（Sec. II）。
- 做到：位置控制的协作臂（Forward Dynamics Compliance Control）在胸腹仿体上，对比基线，含"突发患者运动"条件。
- 依赖：离线参数图；作者自述局限首条为离线黏弹估计。
- 与本文：路线一在超声中的近年代表，但刚度信息离线获得，运动只作为扰动由约束和能量罐兜底。**纠正**：此前记录中写作"在线黏弹参数"，应为离线图 + 在线变刚度。

## Zhang W. et al., "A shared-control ultrasound-guided robotic system for liver puncture under respiratory motion," Scientific Reports 16:26131, 2026, 10.1038/s41598-026-56017-3（SS4c7arItfoJ）
- 机制：共享控制穿刺系统；法向固定参数导纳维持 6 N 接触，切向最小化力；刚度只用于稳定性讨论（线性黏弹接触下闭环无源）。
- 做到：仿生呼吸装置上，呼吸开启后法向力在 0–12 N 波动（p.10）；犬体内定性可行性。
- 与本文：背景，说明固定导纳在呼吸下力波动大；不是最近邻。

## Ebrahimi et al., "Simultaneous Online Registration-Independent Stiffness Identification and Tip Localization of Surgical Instruments in Robot-Assisted Eye Surgery," T-RO 39(2), 2023, 10.1109/TRO.2022.3201393（HokaNUBl3lAJ）
- 机制：在线辨识**器械**弯曲刚度（最小二乘/自适应，视觉提供变形），用刚度把 FBG 测得的巩膜力换算成器械变形，作为 Kalman 滤波的针尖位置量测。
- 与本文：外域可迁移机制（以在线刚度为增益把力换算成位置量测），不作竞争者；它无闭环前馈，不存在本文的增益耦合。只可在方法背景中一句交代。
