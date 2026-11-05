# 场景证据：非学术渠道及补充

检索日期：2026-11-05。以下使用 15 个来源；同行论文仅用于针对性补充，不替代 ARIS 的系统检索。网页或 PDF 可读正文标“全文”，仅读到出版商摘要标“摘要”。未提取原始数据集，不能将数据集介绍当作实测统计。

## 证据表

| 问题 | 量 | 数值及条件 | 来源和位置 | 来源等级 | 访问层级 | 可靠性和局限 |
|---|---|---|---|---|---|---|
| A3 | 腹部扫查接触力 | 47 次检查、10 名健康志愿者、12 名操作者；均值 8.2±4.3 N，范围 1.2–36.5 N；BMI>25 的 4 人均值 10.5 N，BMI 18.5–25 的 6 人均值 7.9 N | [Dhyani 等](https://doi.org/10.3233/WOR-172611)，摘要 Methods/Results | 人体在体实测 | 摘要 | 小样本、腹部；不是各部位的分布或安全阈值，摘要中的范围不能解释为人群百分位。 |
| A3、B2 | 加压和姿势 | 避免过大的探头压力及握力，调整患者、设备以减少不良姿势 | [2017 行业人体工学标准](https://doi.org/10.1177/8756479317725531)，Sonographer/user 部分 | 规范指南 | 全文 | 职业健康建议，不提供定量接触力设定值。 |
| A4、B1 | 接触面几何 | GE 凸阵及线阵的 footprint 见下表 | [GE Versana Premier 探头指南](https://www.gehealthcare.com/-/media/gehc/us/files/products/ultrasound/versana/versana-premier/transducer-guide-versana-premier-r3-july-2024-doc2993174.pdf)，第 1 页表 | 厂商规格 | 全文 | footprint 为厂商标称外形接触尺寸，非受压实际接触斑；不能从声学 FOV 当成曲率半径。 |
| A4 | 凸阵曲率 | Mindray SC6-1s 半径 60 mm；SC10-2RCs 半径 43 mm | [Mindray TE X 探头指南](https://www.mindray.com/content/dam/xpace/en_us/products-solutions/products/ultrasound/point-of-care/tex/5016/TE_X_transducers_40630_Rev_3.0.pdf)，第 1 页相应型号 | 厂商规格 | 全文 | 未给出两型号 footprint；FOV 是成像视野，并非接触面长度。 |
| B2、B1 | 胸部扫描方式 | 相对缓慢移动、肋间暂停，探头垂直胸廓；患者条件许可时扫前、侧、后胸 | [BMUS 肺超声指南](https://www.bmus.org/static/uploads/resources/COVID19__Lung_Ultrasound_BMUS.pdf)，第 4 页 Standard Technique | 规范指南 | 全文 | 定性建议；没有速度、加速度、倾角速率，也不是所有腹部扫查。 |
| A2 | 压痕刚度 | 10 名男性、仰卧、半球压头；放松肌肉小位移下前侧至侧面线性段 0.24–0.34 N/mm，后侧 0.72 N/mm；较大变形分别约 0.41–0.62 和 2.31 N/mm | [Remus 等](https://doi.org/10.3389/fbioe.2024.1384062)，§3.2、§4 | 人体在体实测 | 全文 | 为指定压头与线性分段斜率，不是超声探头的局部切线刚度或通用力—深度关系；身体区域与压入深度重要。 |
| B3 | 呼吸运动 | 5 名志愿者胸腹表面多维跟踪；腹部平均跟踪误差 1.09 mm | [Schaerer 等](https://doi.org/10.1088/0031-9155/57/2/357)，摘要 | 人体在体实测 | 摘要 | 1.09 mm 是跟踪误差，不是呼吸位移幅值；摘要未给所需幅值或速度。 |
| B2、B3 | 公开同步记录 | US-PPNR 声称提供超声、机器人位姿、力矩与光学跟踪同步轨迹，涵盖探头扫查和调整 | [US-PPNR README](https://github.com/PieceZhang/US-PPNR-Dataset/blob/main/README.md)，At a Glance / Dataset Overview | 开源数据集 | 全文（说明） | 仿体/组织及机器人操作，不等同自然人体呼吸；尚未下载验证采样率、标定或轨迹。 |
| B2、A3 | 公开力数据 | ForceU-VLA README 称含肝、肾扫查的力、视觉和动作；目前仅部分公开 | [ForceU-VLA README](https://github.com/VMVLab/ForceU-VLA/blob/main/README.md)，Dataset | 开源数据集 | 全文（说明） | 尚未核实原始数据的受试者、力定义及可用位姿，不直接作为临床分布。 |
| C1 | FT 300 噪声和输出率 | 2017 年 10 月之后的版本 Fx/Fy/Fz 信号噪声 0.1 N；输出 100 Hz；旧版 Fx/Fy 1.2 N、Fz 0.5 N | [Robotiq 手册](https://assets.robotiq.com/website-assets/support_documents/document/online/FT_Sensor_Instruction_Manual_Web_20181218.zip/FT_Sensor_Instruction_Manual_Web/Content/Specifications.htm)，§5.3 | 厂商规格 | 全文 | 噪声定义为 1 秒标准差；输出率不是传感器幅相带宽，不能由此证明 20.8 Hz 测量保真。 |
| C1 | Mini45 分辨率与刚度 | SI-145-5、16-bit DAQ/ECAT/Net 等系统，Fx/Fy/Fz 有效分辨率 1/16 N；计算刚度 Kx/Ky 7.4×10⁷、Kz 9.8×10⁷ N/m | [ATI Mini45](https://www.ati-ia.com/products/ft/ft_models.aspx?id=Mini45)，Ranges and Resolutions / Specifications | 厂商规格 | 全文 | 分辨率按去除四个噪声计数定义，非噪声标准差；传感器本体计算刚度不等于整机末端刚度；未给完整测量链带宽。 |
| C1 | 机器人重复精度 | UR5e ±0.03 mm，ISO 9283 | [UR5e 技术规格](https://www.universal-robots.com/manuals/EN/HTML/SW10_8/Content/prod-usr-man/complianceUR5e/H_g5_sections/appendix_g5/tech_spec_sheet.htm)，Pose repeatability | 厂商规格 | 全文 | 重复精度不是精确度、动态误差或刚度。 |
| C1 | 机器人重复精度及控制刚度 | iiwa 7 R800 / 14 R820 的位置重复精度分别 ±0.1 / ±0.15 mm；另一版手册列 X/Y/Z 可编程笛卡尔刚度 0–5000 N/m | [KUKA 产品表](https://www.kuka.com/-/media/0993c763460846859aff78300f15570f.pdf)，技术表；[KUKA 手册](https://www.kuka.com/-/media/kuka-downloads/imported/87f2706ce77c4318877932fb36f6002d/kuka_lbr_iiwa_ja.pdf)，Technical data | 厂商规格 | 全文 | 可编程阻抗设定范围不是位置接口固有或实测末端法/切向刚度；不同手册对 14 型重复精度标法不同，按具体版本核实。 |
| A2、B1 | 仿体模量 | CIRS 049/049A 背景 18 kPa，病灶 6、9、36、70 kPa；批样为准静态加载测量 | [CIRS 产品规格](https://www.cirsinc.com/wp-content/uploads/2020/03/049-049-DS-030920.pdf)，Specifications 及脚注 | 厂商规格 | 全文 | 杨氏模量不是探头接触刚度；仿体不同于活体，需结合几何和边界条件。 |

## 探头几何表

长度为表中 footprint 的较长尺寸、宽度为较短尺寸；厂家未说明凸面弧长或弦长，勿自行视作精确接触弧长。线阵平面记作不适用；未列出者为未查到。

| 型号 | 类型 | 长度 mm | 宽度 mm | 曲率半径 mm | 来源 |
|---|---|---:|---:|---:|---|
| GE C1-5-RS | 凸阵 | 69.3 | 17.2 | 未查到 | [GE 指南](https://www.gehealthcare.com/-/media/gehc/us/files/products/ultrasound/versana/versana-premier/transducer-guide-versana-premier-r3-july-2024-doc2993174.pdf)，第 1 页 |
| GE 4C-RS | 凸阵 | 66.2 | 18.3 | 未查到 | 同上 |
| GE L4-20t-RS | 线阵 | 48.43 | 11.56 | 不适用 | 同上 |
| GE 12L-RS | 线阵 | 47.1 | 12.7 | 不适用 | 同上 |
| Mindray SC6-1s | 凸阵 | 未查到 | 未查到 | 60 | [Mindray 指南](https://www.mindray.com/content/dam/xpace/en_us/products-solutions/products/ultrasound/point-of-care/tex/5016/TE_X_transducers_40630_Rev_3.0.pdf)，第 1 页 |
| Mindray SC10-2RCs | 凸阵 | 未查到 | 未查到 | 43 | 同上 |

## 各问题结论

- **A2：** 腹部压头实测表明部位、肌肉状态和压入深度均影响斜率；未查到同一临床探头在 4 N 下对应的可靠压入深度和切线刚度。不能将压头斜率直接移植，需结合几何推算并做敏感性分析。
- **A3：** 腹部小样本有实际接触力和 BMI 分组，但未查到跨部位、体型的完整分布；4 N 可以作为被研究的设定，不能称临床标准或典型均值。
- **A4：** 两种凸阵和两种线阵 footprint 及另外两种凸阵半径已有厂商直接规格；同型号同时具有长、宽、半径的资料未查到。不同型号的数据不能拼成一个真实探头。
- **B1：** 胸部肋间暂停有指南支持；未查到可直接使用的肋上局部接触刚度、软组织厚度及肋间距联合数据。跨肋几何和刚度须另行依据解剖资料或推算。
- **B2：** 指南支持慢扫、暂停和调整方向；未查到临床滑动速度、加速度、启停、法向速度及倾角变化的可靠数值分布。公开轨迹值得分析，但现未核验；研究中的 1–11 mm/s 是实验范围，不是行业标准。
- **B3：** 查到人体体表呼吸跟踪方法，但未取得可用的幅值与速度统计；咳嗽、惊跳也未查到直接数据。跟踪误差不能当作位移幅值。
- **B4：** 未查到针对超声探头—耦合剂—人体皮肤、且说明载荷和滑动条件的可信摩擦系数；宜作条件化敏感性分析，不用干皮肤系数冒充。
- **C1：** 可直接使用上述型号的噪声或分辨率规格作为设备参考，但不得混同噪声、分辨率、输出率和带宽；KUKA 的可设定阻抗不代表真实位置控制刚度。UR、Franka 和 iiwa 的同姿态法向/切向实测刚度，以及 FT 300、Mini45 完整链路带宽尚未查到；Franka 和 Kyoto Kagaku 的适用规格本次未核得。

## 建议 ARIS 核读全文

1. Robin Remus, Christian Sure, Sascha Selkmann, Eike Uttich, Beate Bender. **Soft tissue material properties based on human abdominal in vivo macro-indenter measurements** (2024). [DOI: 10.3389/fbioe.2024.1384062](https://doi.org/10.3389/fbioe.2024.1384062)。重点核对力—位移曲线、压头尺寸、各区域和分段定义及临床力区间。
2. Manish Dhyani, Shawn C. Roll, Matthew W. Gilbertson, Melanie Orlowski, Arash Anvari, Qian Li, Brian Anthony, Anthony E. Samir. **A pilot study to precisely quantify forces applied by sonographers while scanning: A step toward reducing ergonomic injury** (2017). [DOI: 10.3233/WOR-172611](https://doi.org/10.3233/WOR-172611)。本次只读摘要；核对力分量、范围统计单位和扫描流程。

