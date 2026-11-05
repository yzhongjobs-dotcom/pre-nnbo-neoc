核查范围：仅读本地原文、代码和联网题录/全文；未修改项目文件、未提交、未运行仿真。联网文献的访问层级已注明。

## A. 接触模型原文核查

### A1. AuSoScan（Duan 等，TMECH 2025）

来源：[`F08.txt`](E:/github_project/scan-impedance-update/literature/extracted_text/field_supplements/F08.txt:316)（本地全文抽取，含行号）；[`F08.pdf`](E:/github_project/scan-impedance-update/literature/papers/field_supplements/F08.pdf)（本地全文 PDF，印刷页 7697–7699，即 PDF 第 5–7 页）。

| 内容 | 原文方程、符号与核查 |
|---|---|
| 状态、输入 | \(x=[q^T,\dot q^T,F^T]^T\)，\(u=\tau\)。\(q,\dot q\)：关节位置、速度；\(F\)：接触力；\(\tau\)：关节力矩。F08.txt L329–334，PDF p.7697。 |
| 机器人动力学 | \(\ddot q=M(q)^{-1}[\tau-C(q,\dot q)\dot q-G(q)-J(q)^TF]\)。\(M\)：质量矩阵；\(C\)：科氏矩阵；\(G\)：重力项；\(J\)：雅可比矩阵。F08.txt L338–340，PDF p.7697。 |
| 静态压力 | \(p(r)=p_s\sqrt{1-r^2/a^2}\)。\(p_s\)：最大压力；\(a\)：接触半径。 |
| Hertz 最大压力 | \(p_s=\frac{E^*}{2}\sqrt{d/R}\)。\(d\)：压入深度；\(R\)：等效半径；\(E^*\)：等效弹性模量。 |
| 等效模量、半径 | \(\frac1{E^*}=\frac{1-\nu_1^2}{E_1}+\frac{1-\nu_2^2}{E_2}\)；\(E^*=E_1/(1-\nu_1^2)\)（探头刚性，\(E_2\to\infty\)）；\(\frac1R=\frac1{R_1}+\frac1{R_2}\)，平面背部 \(R_1\to\infty\)，故 \(R=R_2\)。\(E_i\)：杨氏模量；\(\nu_i\)：泊松比。F08.txt L377–413，PDF p.7698。 |
| 滑动压力修正，式(11) | \[
p(r)=
\begin{cases}
p_m\sqrt{1-r^2/a_1^2},&-a_1\le r\le0\\
p_m\sqrt{1-r^2/a_2^2},&0<r\le a_2
\end{cases}
\] \(a_1\)：前方半径，\(a_2\)：后方半径，\(p_m\)：滑动时最大压力。作者假设滑动使前方接触区增大、后方减小。F08.txt L364–380，PDF p.7698。 |
| 前后半径，式(12) | \(a_1=\sqrt{Rd_0}-k_1v,\quad a_2=\sqrt{Rd_0}+k_2v\)。\(v\)：滑动速度；\(k_1,k_2\)：速度—半径线性系数。 |
| 静态压入，式(13) | \(d_0=4F/(\pi E^*L)\)。\(F\)：给定法向力；\(L\)：探头宽度。 |
| 法向力积分，式(14) | \[
F=\int_0^L\int_{-a_1}^{a_2}p(r)\,dr\,dl
=\frac{\pi Lp_m}{4}(a_1+a_2)
=\frac{\pi LE^*}{4}(2\sqrt{Rd_0}+\Delta kv)\sqrt{d/R}
\] \(\Delta k=k_2-k_1\)。F08.txt L387–415，PDF p.7698。 |
| 法向力，式(15) | \[
F=\frac{\pi LE^*}{8}(2\sqrt{Rd_0}+\Delta kv)\sqrt{d/R}.
\] |
| 摩擦 | \(f=f_d+f_b\)，\(f_d=\mu F\)。\(\mu\)：动摩擦系数；\(f_b\)：形变引起的内部摩擦。 |
| 摩擦积分，式(16) | \[
f_b=\int_0^L\int_{-a_1}^{a_2}
\frac{rP(r)}{\sqrt{R^2-r^2}}\,dr\,dl
\] \[
=\int_0^L\int_{-a_1}^{0}
\frac{p_mr}{a_1}\sqrt{\frac{a_1^2-r^2}{R^2-r^2}}\,dr\,dl+
\int_0^L\int_0^{a_2}
\frac{p_mr}{a_2}\sqrt{\frac{a_2^2-r^2}{R^2-r^2}}\,dr\,dl .
\] F08.txt L428–445，PDF p.7699。 |
| 法向动态，式(17) | 原文为 \[
\dot F=\frac{\pi LE^*\Delta k}{8}
\left(\dot v\sqrt{d/R}+
\frac{v\dot d}{2\sqrt{Rd}}\right).
\] |
| 摩擦动态 | \(\dot f=\mu\dot F+\dot f_b\)。F08.txt L447–460，PDF p.7699。 |
| 力矩 | 未给出接触压力积分得到的接触力矩方程；只在硬件中使用六轴 F/T 传感器。F08.txt L195–201。 |

**式(15)–(17)一致性：不一致。**

由式(15)直接求导应为：

\[
\dot F=\frac{\pi LE^*}{8}
\left[
\Delta k\dot v\sqrt{d/R}
+\frac{(2\sqrt{Rd_0}+\Delta kv)\dot d}{2\sqrt{Rd}}
\right].
\]

原式(17)遗漏了 \(2\sqrt{Rd_0}\) 项，并把 \(\Delta k\) 提取到整个括号外。特别地，\(v=0\) 时：

- 式(15)仍为 \(F=\frac{\pi LE^*}{4}\sqrt{Rd_0}\sqrt{d/R}\)，随 \(d\) 变化；
- 原式(17)给出 \(\dot F=0\)。

因此式(17)不是式(15)的一般导数；不能将其直接作为一致的动态方程。

**模型针对的问题与假设**

- 物理问题：机器人超声探头沿软背部表面滑动时的法向变形阻力和切向摩擦。
- 几何：圆柱形探头、人体背部或仿真体表面视为平面；探头刚性；小范围滑动修正。
- 接口：状态含 \(q,\dot q,F\)，输入为关节力矩 \(\tau\)；使用六轴 F/T 传感器，1 kHz。
- 未建模：探头—皮肤黏附、摩擦静态分支、接触力矩、组织空间各向异性和动态材料记忆。
- 作者承认：复杂场景、严重脊柱侧弯和突发扰动下的适应性仍不足，未来需要更柔顺、反应性更强的平台。F08.txt L659–676，PDF p.7702。

**辨识、数值与实验**

| 项目 | 原文数字 |
|---|---|
| 参数辨识 | 滑动速度 \(0.2\)–\(1\) m/s；法向力 3、6、9、12 N；最小二乘辨识 \(k_1,k_2\)。该速度仅用于辨识，不是扫描速度。F08.txt L539–575，PDF p.7701。 |
| 参数 | \(\mu=0.39\)，\(k_1=0.001\)，\(k_2=0.003\)，\(E^*=6.3\) kPa，\(R=0.03\) m。F08.pdf 本地全文 p.7701，表 I。 |
| 探头 | Scolioscan；孔径约 80 mm × 15 mm；F/T 传感器 FT300，1 kHz。F08.txt L190–201。 |
| 跟踪误差 | 无形变模型/有形变模型：Step 位置 0.97/0.53 cm、力 0.357/0.106 N；Slope 0.83/0.39 cm、力 0.393/0.185 N；Wave 1.51/1.02 cm、力 0.959/0.535 N。F08.pdf p.7701，表 II。 |
| 扫描实验 | 琼脂皮肤仿真体和一名人体对象；作者未给出实际扫描速度数字。F08.txt L603–614。 |

---

### A2. T-RO 2023：Wijayarathne 等

来源：[`wijayarathne.txt`](E:/github_project/scan-impedance-update/literature/extracted_text/wijayarathne.txt:197)（本地全文抽取）；[`wijayarathne.pdf`](E:/github_project/scan-impedance-update/literature/papers/wijayarathne.pdf)（本地全文 PDF，印刷页 3552–3554、3559–3562）；DOI：[10.1109/TRO.2023.3286070](https://doi.org/10.1109/TRO.2023.3286070)。

| 内容 | 原文方程、符号与核查 |
|---|---|
| Hertz 法向 | \[
d=\left(\frac{9F^2}{16E^2R}\right)^{1/3}.
\] \(d\)：中心压入；\(F\)：工具施加的总力；\(E\)：约化杨氏模量；\(R\)：球形工具半径。 |
| 约化模量 | \[
\frac1E=\frac{1-\nu_1^2}{E_1}+\frac{1-\nu_2^2}{E_2};
\quad E\approx E_2/(1-\nu_2^2)
\] 假设工具刚性。 |
| 应力张量 | \[
\sigma=
\begin{bmatrix}
\sigma_r&0&\sigma_{rz}\\
0&\sigma_\theta&0\\
\sigma_{zr}&0&\sigma_z
\end{bmatrix}.
\] |
| 法向应力 | \(\mathbf n=[\sin\theta,0,\cos\theta]^T\)，\(\sigma_n=\mathbf n^T\sigma_c\mathbf n\)。作者给出柱坐标到笛卡尔坐标变换后的展开式。 |
| 摩擦积分 | \[
df=\mu\sigma_n\,dS=\mu\sigma_n\frac{2\pi r\,dr}{\cos\theta},
\quad
f=\int df\cos\theta=2\pi\mu\int_0^a\sigma_n r\,dr.
\] \(a=\sqrt{Rd}\)，\(\mu\)：摩擦系数。 |
| 面积/压深修正 | \[
F_f=f=\mu F_z\left[1+(2\nu-1)\frac{3a^2}{10R^2}\right]\mathbf n_v+k_dv_e.
\] \(F_z=F\cos\theta_F\)：法向力；\(\mathbf n_v\)：切向速度单位向量；\(v_e=\|J\dot q\|\)：末端速度；\(k_d\)：切向阻尼系数。 |
| 法向动态 | \[
\dot z=-\dot d
=-\left(\frac1{6E^2RF_z}\right)^{1/3}\dot F_z.
\] |
| 摩擦动态 | \[
\dot F_f
=\left[
\mu\dot F_z+
\frac{3\mu(2\nu-1)}{10R}
(\dot F_zd-F_z\dot z)
\right]\mathbf n_v+k_d\dot v_e.
\] 等价地使用 \(\dot d=-\dot z\)。 |
| 合并接触动力学 | \[
\dot{\mathbf F}_e
=(6E^2RF_z)^{1/3}\dot d\,\mathbf n_z+
\left[
\mu\dot F_z+
\frac{3\mu(2\nu-1)}{10R}
(\dot F_zd+F_z\dot d)
\right]\mathbf n_v
+k_d\dot v_e.
\] |
| 机器人动力学 | \[
\ddot q=M(q)^{-1}
[\tau_u-C(q,\dot q)\dot q-G(q)-J^T\mathbf F_e].
\] \(M,C,G,J,\tau_u\) 含义同通常机器人动力学；文中为 7 自由度 KUKA。 |
| 轨迹—接触耦合 | 状态 \(x=[q^T,\dot q^T,\mathbf F_e^T]^T\)，输入 \(u=\tau_u\)，离散为 \(x[i+1]=f(x[i],u[i])\)。代价同时包括法向力误差、末端位姿误差和力矩正则项。 |
| 曲率约束 | \[
\frac{J^{-T}M(q)J^{-1}\|J\dot q\|^2}{\kappa}
\le \mu\,\mathbf N^T\mathbf F_e\,\mathbf N.
\] \(\kappa\)：轨迹曲率；\(\mathbf N\)：表面平均法向；实验中 \(\mathbf N=[0,0,1]^T\)。 |

方程和假设来源：wijayarathne.txt L197–340，PDF p.3552–3554；轨迹耦合与约束：wijayarathne.txt L295–380，PDF p.3553–3554。

**物理问题与假设**

- 物理问题：已知球形工具与软材料接触时，同时规划运动轨迹和法向/摩擦力。
- 工具：球形、刚性、明显比软材料硬。
- 运动：曲线路径、近似匀速；主要研究滑动摩擦。
- 忽略：黏附、由变形引起的滚动摩擦；作者明确说这些只在特定工具和应用中作为次要因素。
- 材料：已知或预先辨识的 \(E,\nu,\mu,k_d\)。
- 不等同于线阵探头：原实验是球形压头和 EcoFlex 软基底。

**辨识与实验数字**

| 项目 | 原文数字 |
|---|---|
| 参数辨识 | 沿材料表面预设运动，记录 ATI Mini45 的力/矩；杨氏模量由球形压头循环线性压入、非线性最小二乘得到；摩擦参数由切向运动数据拟合。wijayarathne.txt L674–688，PDF p.3559。 |
| 参数 | \(\mu=0.4512,\ k_d=13.1315,\ R^2=0.9103\)。wijayarathne.txt L699–706，PDF p.3560；图 8。 |
| 平台 | KUKA 机械臂、球形工具、ATI Mini45、EcoFlex 软材料；动态扰动平台 1 Hz。wijayarathne.txt L605–607、L659–683；PDF p.3558–3559。 |
| 计算 | 高层 MPC 5 Hz、1 s horizon、\(\Delta t=0.02\) s；低层力控制 100 Hz；1 s 轨迹平均计算约 150 ms，未使用 AutoDiff/解析模型约 1300 ms。wijayarathne.txt L483–506，PDF p.3556。 |
| 周期扰动误差 | Motion Path 1/2/3，MPC+FC 的路径/力 RMSE：0.0164 m/0.283 N、0.0191 m/0.319 N、0.0233 m/0.0084 N。PDF p.3562，表 III。 |
| 脉动平台 | MPC+FC 的力 RMSE：0.573、0.852、0.323 N；路径 RMSE：0.0258、0.0384、0.0153 m。PDF p.3561，表 II。 |

**作者承认局限**

- 接触模型只在局部范围有效；
- 大变形需要计算更重的 FEM；
- 高变形和润滑条件需要额外学习/建模变形摩擦；
- 当前机器人阻抗预设，尚未随任务调整；
- 未来需要更丰富材料模型、实时参数估计和三维曲面规划。

来源：wijayarathne.txt L786–826，PDF p.3562。

---

### A3. Hunt–Crossley（HC）

#### 基本方程

Moreira 等的本地全文表 1 给出：

\[
f(t)=
\begin{cases}
k\,x^n(t)+\eta\,x^n(t)\dot x(t),&x(t)\ge0\\
0,&x(t)<0.
\end{cases}
\]

其中：

- \(x\)：压入量；
- \(\dot x\)：压入速度；
- \(k\)：非线性弹性系数；
- \(\eta\)：速度耗散系数；
- \(n\)：压入非线性指数；
- \(f\)：法向接触力。

来源：[`F10.txt`](E:/github_project/scan-impedance-update/literature/extracted_text/field_supplements/F10.txt:178)，本地全文 PDF 印刷页 55。

HC 本身没有切向摩擦、接触力矩或机器人动力学方程。机器人动力学必须另行耦合。

#### Moreira 2014 的实际范围

来源：[`F10.txt`](E:/github_project/scan-impedance-update/literature/extracted_text/field_supplements/F10.txt:218)，本地全文 PDF 印刷页 55–57；DOI：[10.1016/j.cmpb.2013.11.005](https://doi.org/10.1016/j.cmpb.2013.11.005)。

| 项目 | 核查结果 |
|---|---|
| 针对问题 | 为软组织力控制选择适合实时应用的低复杂度组织模型。 |
| 实验 | 离体猪心和牛肉；每个试样 4 次 relaxation test；位置阶跃分别为 10 mm 和 8 mm。F10.txt L228–247、L289–310。 |
| 参数辨识 | 线性最小二乘；为统一辨识流程固定 \(n=1.26\)；该数值来自作者引用的既有结果。F10.txt L273–287。 |
| HC 数值 | 查到 \(n=1.26\)；HC 的 \(k,\eta\) 单独数值未在本地全文中查到，记为“未查到”。 |
| HC 误差 | 平均力估计误差 \(0.091\pm0.030\) N。F10.txt L258–271，表 2。 |
| 结论 | Kelvin–Boltzmann、HC 和 fractional 模型较好；HC 和 Kelvin–Boltzmann 对瞬态较好。 |
| 控制器 | 作者最终选 Kelvin–Boltzmann，而不是 HC，作为模型力控制器。静态组织参数 \(\alpha=34.3,\beta=200.0,\gamma=0.024\)；PVC 参数 \(\alpha=5.0,\beta=290.0,\gamma=0.02\)。F10.txt L803–851。 |
| 局限 | 试样为离体组织；relaxation test 对低频行为更敏感；非均质性导致参数随接触点变化；HC 的非线性和参数辨识复杂。F10.txt L218–247、L332–335。 |

#### Haddadi 2012

来源：DOI：[10.1109/TRO.2012.2183054](https://doi.org/10.1109/TRO.2012.2183054)；[题录/摘要页](https://www.researchgate.net/publication/254060226_Real-Time_Identification_of_Hunt-Crossley_Dynamic_Models_of_Contact_Environments)。当前可打开内容为题录和摘要，不是全文。

- 针对未知接触环境的实时动力学辨识；
- 比较 HC 单阶段、HC 双阶段和 Kelvin–Voigt RLS；
- 实验对象包括橡胶球、海绵、PVC 仿真体以及含硬夹杂的 PVC 仿真体；
- 摘要未给出可直接迁移的 \(k,\eta,n\) 数值，全文未查到。

#### Courrèges 2015/2016

来源：Springer 正式页面：[10.1007/978-3-319-22368-1_21](https://link.springer.com/chapter/10.1007/978-3-319-22368-1_21)。当前打开的是题录和摘要预览，不是全文。

- 问题：人体腹部接触阻抗的在体评估；
- 方法：比较原始非线性 HC 与 log-linearized HC；
- 摘要结论：log-linearized HC 应被视为独立接触模型，可同时用于辨识和控制；
- 在体实验确实存在，但摘要未给出 \(k,\eta,n\)、压入范围或速度数字；这些数值“未查到”。

---

## B. 项目场景中的物理过程

| 过程 | 判断 | 实测依据、大小与边界 |
|---|---|---|
| 1. 探头—皮肤滑动摩擦、静/动摩擦差、stick-slip | **真实存在；量级未确定** | AuSoScan 实测仿真体给 \(\mu=0.39\)，T-RO 的 EcoFlex 球形工具给 \(\mu=0.4512\)，但都不是“线阵探头+人体+耦合剂”的实测值。T-RO 明确将 stiction、润滑和材料不均匀性列为模型误差来源；液体或微粒会增加摩擦不均匀性。T-RO 全文 PDF p.3560。项目场景的 \(\mu_s,\mu_k\) 和 stick-slip 频率：**未查到**。 |
| 2. 压力不对称、前方堆积、速度相关阻力 | **真实存在或至少有直接物理模型支持；人体线阵探头的直接压力场测量未查到** | AuSoScan 明确假设滑动时前方接触区增大、后方减小，并给出式(11)–(16)；T-RO 的摩擦项随接触面积/压入量和速度变化。T-RO 实验观察到法向力增大时形变摩擦占主导。直接的探头下压力分布测量：**未查到**。 |
| 3. 探头轴线与表面法向不一致及轴向投影 | **真实存在；典型误差约数度** | Jiang 等在体组织上的探头法向估计平均绝对角差 \(2.9\pm1.6^\circ\)，仿真体 \(2.2\pm1.5^\circ\)；来源为 IEEE 题录/摘要页：[DOI 10.1109/TIE.2020.3036215](https://doi.org/10.1109/TIE.2020.3036215)。SAPM 临床数据要求约 \(5.08^\circ\) 的旋转调整；来源为 PMC 全文：[SAPM](https://pmc.ncbi.nlm.nih.gov/articles/PMC7614654/)。若偏角为 \(\theta\)，切向力轴向投影为 \(F_\mathrm{axis}=F_n\cos\theta+F_t\sin\theta\)；这是坐标分解，不是额外实测结论。 |
| 4. 机器人/工具切向和倾斜柔顺 | **真实存在；当前指定文献未给出协作臂末端六维刚度实测值** | T-RO 使用低阻抗和 100 Hz admittance 控制，并指出预设阻抗限制性能。AuSoScan 使用六轴 F/T，但未报告末端刚度。FARUS 报告探头弹簧连接 \(k=500\) N/m，但这是连接机构刚度，不是协作臂末端刚度：[正式全文/预印本](https://assets-eu.researchsquare.com/files/rs-1888871/v1_covered_93490491-f4e0-4f91-bdbf-ed31cb30b30e.pdf)。项目机器人的切向刚度、姿态刚度和切向载荷引起的位置误差：**未查到**。 |
| 5. 实际扫描速度及变化 | **真实存在，速度范围约 1–60 mm/s，启停和人为变化明显** | 人工扫描完整测量：中位速度 13.8、6.1、15.8 mm/s；最大速度 29.0、15.8、30.1 mm/s；在体纵向平均 \(38.4\pm4.8\) mm/s、最大 47.3 mm/s，径向平均 \(43.4\pm8.2\) mm/s、最大 58.9 mm/s。来源：[Frontiers 全文](https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2022.926006/full)。机器人例子：1 mm/s（PMC 全文，[thyroid scanning](https://pmc.ncbi.nlm.nih.gov/articles/PMC7392820/)）、5 mm/s（Virga 全文 PDF）、3.6 mm/s（JHU 正式学位论文链接）、AuSoScan 未报告扫描速度。 |
| 6. 法向力—压入非线性、toe 区、硬化、黏弹记忆、速率项 | **非线性、记忆和速率项真实存在；在体定量边界不足** | Moreira 的离体猪心/牛肉 relaxation tests 支持 HC、Kelvin–Boltzmann 和 fractional 模型；Courrèges 的 Springer 摘要支持人体腹部在体 HC/log-linearized HC 评估，但具体 \(k,\eta,n\) 未公开于可访问摘要。AuSoScan 和 T-RO 也使用 Hertz/HC 型速率项。指定材料中人体线阵探头的 toe 区和硬化曲线：**未查到**。 |
| 7. 六维力/矩传感器与可用量 | **研究型机器人超声中常见，但不是“标配”结论** | AuSoScan：FT300 六轴、1 kHz；T-RO：ATI Mini45；SAPM：六轴传感器并记录力矩；Virga：KUKA + F/T。SAPM 临床数据给出轴向力 8.01 N、绕 \(x/y\) 轴力矩 150 mN·m 的能力需求。六维传感器可提供切向力和力矩，但是否用于闭环取决于系统。项目当前仿真中切向/力矩通道没有真正进入控制。 |

---

## C. 当前仿真代码对照

来源：[`physics_plant.py`](E:/github_project/scan-impedance-update/experiments/rossi_failure_study/physics_plant.py:1)、[`scenes.py`](E:/github_project/scan-impedance-update/experiments/rossi_failure_study/scenes.py:1)、[`plant.py`](E:/github_project/scan-impedance-update/experiments/cobot_scan_benchmark/plant.py:1)。

| 过程 | 当前状态 |
|---|---|
| 法向有限接触面积 | **已包含**。`physics_plant.py` L45–57 的 `_fe()` 使用有限足迹积分；`scenes.py` L56–59 提供空间变化模量。 |
| AuSoScan 圆柱 Hertz 压力和前后不对称 | **缺失**。当前 `_fe()` 没有 \(a_1,a_2\)、\(\Delta k v\) 或式(11)压力分布。`law>0` 时 L50 使用统一 \(d^\text{law}\) 幂律；`law=1` 也不等于 AuSoScan 式(15)。 |
| T-RO 球形 Hertz 法向 | **未按原式实现**。当前默认是有限足迹积分；没有 \(d=(9F^2/16E^2R)^{1/3}\) 的球形接触接口。 |
| HC 法向速率项 | **部分包含但方程不同**。`physics_plant.py` L93、L101 使用 \(F=(f_e+\sum q_i)(1+\lambda\dot d)\)，不是 \(K d^n(1+B\dot d)\)。 |
| 黏弹记忆 | **已包含候选形式**。`physics_plant.py` L8–10、L89–101 使用三支 Prony 状态；`scenes.py` L138–140 给出 \(\alpha=[.35,.35,.35]\)、\(\tau=[.008,.08,.8]\)。这不是 Moreira 的 KB，也不是已辨识人体参数。 |
| 组织运动 | **已包含**。`physics_plant.py` L79–82 有身体动力学、表面 \(h(x)\)、编码器与接触深度关系；`scenes.py` L75–122 有心动、呼吸、推移、flinch、cough。 |
| 扫描速度 | **已包含固定速度**。`scenes.py` L54–55 为随机 0.01–0.03 m/s，即 10–30 mm/s。启停/换向速度轨迹没有独立建模；`physics_plant.py` L77 使用恒定 \(v_\mathrm{scan}\)。 |
| 压力不对称和前方堆积 | **缺失**。空间模量变化不是滑动造成的压力前后偏置。 |
| 摩擦系数、静/动摩擦、stick-slip | 当前 `physics_plant.py` **缺失**；其日志 L12–13 还注明 `Ft`、`My` 未使用。旧 `plant.py` L142–145 有 Stribeck 型切向摩擦和速度阻尼，但它属于旧对象，不是当前 Rossi plant。 |
| 切向力 | 当前 `physics_plant.py` **缺失**。`scenes.py` 的 `Ft` 没有物理计算。旧 `plant.py` L142–145 有 \(F_t\)，但未反作用到位置/法向力。 |
| 力矩 | 当前 `physics_plant.py` **缺失**；日志只保留占位列。旧 `plant.py` L144–145 有由足迹一阶矩和表面斜率构成的 \(M_y\)，但它没有机器人姿态动力学。 |
| 表面法向偏角/投影 | 当前 **缺失**。`physics_plant.py` 只使用标量法向接触；没有 \(F_t\sin\theta\)、姿态状态或三维法向。旧 `plant.py` L145 使用表面斜率 `hx` 对切向力和力矩作简化修正。 |
| 机器人切向/姿态柔顺 | 当前 **缺失**。`physics_plant.py` L94–96 只有法向串联刚度 \(k_r\)；没有切向刚度、姿态刚度或六维刚体动力学。 |
| 传感器六维可用量 | 当前 **法向力已使用**，见 `physics_plant.py` L107–110；切向力和力矩没有进入 `Obs`。旧 `plant.py` L216–218 暴露 `ft`、`my`，但当前 Rossi plant 没有对应通道。 |
| 代码中的参数真实性 | `scenes.py` L56、L134–140 的模量、身体质量、共振频率、阻尼、工具刚度和记忆参数均是仿真生成值；代码注释 L8–9 也明确部分量“没有直接测量，始终扫描”。 |

---

## D. 若加入缺失过程：方程与参数来源

不指定估计算法，也不偏向某一接触模型。

| 过程 | 建议加入的方程形式 | 参数及来源 |
|---|---|---|
| 法向非线性 | 并列保留三种可比较形式：\(F=K d^n\)；\(F=K d^n(1+B\dot d)\)；有限足迹积分 \(F=\int_Ap(r,x,t)\,dA\)。 | HC 指数 \(n=1.26\)：Moreira 本地全文 F10.txt L278–287；\(K,B\) 的项目可迁移数值：**未查到，需扫描**。 |
| AuSoScan 前后压力不对称 | 使用原文式(11)–(16)：\(p(r)\) 分段，\(a_1=\sqrt{Rd_0}-k_1v\)，\(a_2=\sqrt{Rd_0}+k_2v\)，\(F=\int p\,dA\)，\(f_b=\int rP(r)/\sqrt{R^2-r^2}\,dA\)。 | \(\mu=.39,k_1=.001,k_2=.003,E^*=6.3\) kPa，\(R=.03\) m：AuSoScan PDF p.7701；人体线阵探头适用性：**需验证，不能直接移植**。 |
| T-RO 接触面积/压深摩擦 | \[
F_f=\mu F_z\left[1+(2\nu-1)\frac{3a^2}{10R^2}\right]\mathbf n_v+k_dv_e,\quad a=\sqrt{Rd}.
\] | \(\mu=.4512,k_d=13.1315,R^2=.9103\)：T-RO PDF p.3559–3560；球形工具半径和材料参数对线阵探头：**未查到，需扫描**。 |
| 静摩擦、动摩擦、stick-slip | 可用 \[
F_t=-\left[\mu_k+(\mu_s-\mu_k)e^{-(v/v_s)^2}\right]F_n\tanh(v/v_s)-c_tv.
\] | AuSoScan/T-RO 只有动摩擦或等效摩擦系数；\(\mu_s,\mu_k,v_s,c_t\) 的“耦合剂+人体+线阵探头”实测：**未查到，需扫描**。 |
| 压力引起接触力矩 | \[
\mathbf M_c=\int_A\mathbf r\times[-p(\mathbf r)\mathbf n+\mathbf t(\mathbf r)]\,dA.
\] | 需要探头接触面几何、压力质心和 F/T 坐标变换；本地三篇文献未给线阵探头实测压力质心，**需扫描/标定**。 |
| 探头偏角和轴向投影 | \[
F_\mathrm{sensor,z}=F_n\cos\theta+F_t\sin\theta,\quad
\mathbf M=\mathbf r_c\times\mathbf F.
\] | 在体角度误差 \(2.9\pm1.6^\circ\)：Jiang DOI 摘要；SAPM 旋转需求约 \(5.08^\circ\)：PMC 全文。项目路径误差分布：**未查到，需扫描**。 |
| 机器人切向柔顺 | \[
M_t\ddot{\mathbf x}_t+D_t\dot{\mathbf x}_t+
K_t(\mathbf x_t-\mathbf x_{t,\mathrm{cmd}})
=\mathbf F_t.
\] | FARUS 弹簧连接 \(500\) N/m，但非末端刚度；实际协作臂 \(K_t,D_t,K_\theta,D_\theta\)：**未查到，需扫描**。 |
| 姿态柔顺 | \[
\mathbf I\dot{\boldsymbol\omega}
+D_\theta\boldsymbol\omega+
K_\theta(\boldsymbol\theta-\boldsymbol\theta_\mathrm{cmd})
=\mathbf M_c.
\] | 六轴传感器可测 \(M_x,M_y,M_z\)，但实际姿态刚度未查到；**需扫描**。 |
| 扫描启停/换向 | 采用梯形或 S 曲线：\[
\dot x=v(t),\quad
v(t)\in[0,v_\mathrm{max}],\quad
|\dot v|\le a_\mathrm{max}.
\] | 机器人文献：1、3.6、5 mm/s；人工中位 6.1–15.8 mm/s、最大约 30.1 mm/s；在体实验平均约 38–43 mm/s。\(a_\mathrm{max}\) 和换向曲线：**未查到，需扫描**。 |
| 组织记忆 | 可用 Prony/SLS：\[
\dot q_i=\alpha_i\dot F_e-q_i/\tau_i,\quad
F=F_e+\sum_iq_i;
\] 或 KB：\[
F+\gamma\dot F=K d+\eta\dot d.
\] | Moreira 的 KB：\((\alpha,\beta,\gamma)=(34.3,200,.024)\) 和 PVC \((5,290,.02)\)，但不是人体线阵探头参数；项目 \(\tau_i\) 需以在体/仿真体数据为依据，当前没有直接来源，**需扫描**。 |
| 速率效应 | 可分别比较 \(F=K d^n(1+B\dot d)\)、\(F+\gamma\dot F=Kd+\eta\dot d\)、Prony 记忆三种形式。 | T-RO \(k_d=13.1315\) 是 EcoFlex 球形工具结果；人体线阵探头速率系数：**未查到，需扫描**。 |
| 法向硬化/toe 区 | 分段形式：\[
F(d)=
\begin{cases}
K_1d^{n_1},&d<d_c\\
F(d_c)+K_2(d-d_c)^{n_2},&d\ge d_c.
\end{cases}
\] | 指定 AuSoScan、T-RO、Moreira 原文未提供人体在线阵探头的 toe 区 \(d_c,K_1,K_2,n_1,n_2\)，**未查到，需扫描**。 |
| 六维测量 | \[
\mathbf y_F=
[\mathbf F_c;\mathbf M_c]+\mathbf b+\boldsymbol\eta,
\quad
\mathbf F_c=[F_x,F_y,F_z]^T.
\] | AuSoScan FT300 1 kHz；T-RO ATI Mini45；SAPM 六轴传感器；具体噪声、偏置随机游走和坐标变换参数：**需由实际传感器数据给出**。 |

**结论性核查：** 当前 Rossi 仿真已经有法向有限足迹、空间模量变化、组织运动、Prony 记忆和法向速率项，但没有当前任务最关键的切向力—法向投影—接触力矩—机器人切向/姿态柔顺闭环。旧 `cobot_scan_benchmark/plant.py` 有简化 Stribeck 摩擦和力矩，但不能视为当前模型已包含这些过程。