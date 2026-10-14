"""One-time, lossless organization of the R15 handoff package."""
from pathlib import Path
import csv
import hashlib
import json
import re
import os

ROOT = Path(__file__).resolve().parents[1]
RECORD = ROOT / '90_整理记录'
METHOD = '01_当前方法与证明'
FIELD = '02_研究现状与文献'
EXPERIMENT = '03_实验与后续工作'
RESOURCE = '04_规范与外部资源'
SOURCE = '05_原文摘录与阅读材料'
HISTORY = '99_历史归档'
MAIN = 'R15_RESEARCH_PROBLEM_SPATIAL_PREVIEW_CONTACT_MODEL_2026-09-12.md'
APP = 'R15_RESEARCH_PROBLEM_SPATIAL_PREVIEW_CONTACT_MODEL_APPENDIX_2026-09-12.txt'
PROOF = 'R15_ANEOC_CERTIFIED_DEPLOYMENT_PROOFS_2026-09-13.txt'
if (RECORD / '文件迁移清单.csv').exists():
    raise SystemExit('Already organized; refusing to run again.')

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def checked(rel):
    p = (ROOT / rel).resolve()
    if not p.is_relative_to(ROOT) or p == ROOT:
        raise ValueError(f'Path outside intended workspace: {p}')
    return p

def destination(rel):
    if rel.startswith('tmp/'):
        return f'{SOURCE}/{rel}'
    if rel.startswith('补充规则/'):
        return f'{RESOURCE}/{rel}'
    if rel.startswith('历史版本/'):
        return f'{HISTORY}/05_原有历史版本/{rel[len("历史版本/"):]}'
    if rel.startswith('__pycache__/'):
        return f'{HISTORY}/06_旧路线代码/{rel}'
    name = Path(rel).name
    if name in (MAIN, APP, PROOF):
        return f'{METHOD}/{name}'
    if name in ('02_FIELD_STATE_ANALYSIS_AND_FINAL_RESEARCH_PROBLEM.txt',
                '05_EVIDENCE_AND_COMPARATOR_INDEX.txt',
                'R15_STRONGEST_PRIOR_FINAL_AUDIT_2026-09-07(2).txt',
                'R15_CONTACT_WRENCH_PREVIEW_LITERATURE_AND_CONTROL_FUSION_2026-09-08.txt'):
        return f'{FIELD}/{name}'
    if name == '06_E2E_HARNESS_REUSE_GUIDE.txt':
        return f'{RESOURCE}/{name}'
    if name == 'verify_ds_aneoc_nnbo_derivatives.py':
        return f'{HISTORY}/06_旧路线代码/{name}'
    if name.startswith(('00_', '01_', '03_', '04_', '07_', '08_')) or 'HANDOFF' in name:
        return f'{HISTORY}/01_历史交接与状态/{name}'
    if name in ('如何理解新method.md', 'R15_METHOD_REPORT_EXPLANATION_2026-09-11.md'):
        return f'{HISTORY}/04_历史解释与汇报/{name}'
    category = '03_历次审计' if ('AUDIT' in name or 'GATE' in name) else '02_历次方法与阶段稿'
    match = re.search(r'2026-09-\d{2}', name)
    date = match.group() if match else '未标日期'
    return f'{HISTORY}/{category}/{date}/{name}'

originals = sorted(p for p in ROOT.rglob('*') if p.is_file() and RECORD not in p.parents)
rows = []
for p in originals:
    old = p.relative_to(ROOT).as_posix()
    new = destination(old)
    target = checked(new)
    if target.exists():
        raise FileExistsError(target)
    rows.append({'原路径': old, '新路径': new, '字节数': p.stat().st_size, 'SHA256': digest(p)})
if len({r['新路径'] for r in rows}) != len(rows):
    raise ValueError('Duplicate destination.')

# Directory renames retain their internal structure, including empty directories.
directory_moves = [
    ('tmp', f'{SOURCE}/tmp'),
    ('补充规则', f'{RESOURCE}/补充规则'),
    ('历史版本', f'{HISTORY}/05_原有历史版本'),
    ('__pycache__', f'{HISTORY}/06_旧路线代码/__pycache__'),
    ('R15_STAGE1_CONTACT_CLOSURE_2026-09-12', f'{HISTORY}/07_原有空目录/R15_STAGE1_CONTACT_CLOSURE_2026-09-12'),
]
for old, new in directory_moves:
    checked(old)
    if checked(new).exists():
        raise FileExistsError(new)

with (RECORD / '文件迁移清单.csv').open('w', encoding='utf-8-sig', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['原路径', '新路径', '字节数', 'SHA256'])
    writer.writeheader()
    writer.writerows(rows)

for old, new in directory_moves:
    src, dst = checked(old), checked(new)
    if src.exists():
        dst.parent.mkdir(parents=True, exist_ok=True)
        src.rename(dst)
for row in rows:
    src = checked(row['原路径'])
    if src.parent == ROOT:
        dst = checked(row['新路径'])
        dst.parent.mkdir(parents=True, exist_ok=True)
        src.rename(dst)

def write(rel, text):
    p = checked(rel)
    if p.exists():
        raise FileExistsError(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text.strip() + '\n', encoding='utf-8')

def link(label, rel):
    target = rel.replace(' ', '%20').replace('(', '%28').replace(')', '%29')
    return f'[{label}]({target})'

write('README.md', f'''
# R15 项目入口

整理日期：2026-10-09。现有研究记录的最新节点为 2026-09-13；本次仅整理文件与已有结论，不重新开展文献检索、方法推导或实验。

**当前主线：空间接触预瞄 → 柱面 Hertz 模型与 MHE → NNBO/single-critic HJB → ANEOC 候选直接认证发布 → 实际阻抗。Newton 在需要时校正。**

## 从这里恢复进度

1. {link('当前进度与版本判断', METHOD + '/00_当前进度与阅读顺序.md')}：五分钟恢复项目状态。
2. {link('完整方法正文', METHOD + '/' + MAIN)}：第 1—2 节总问题，第 3 节接触与估计，第 4 节最优控制，第 5 节直接策略更新，第 6 节实际执行。
3. {link('推导附录', METHOD + '/' + APP)}：A—G 接触与估计，H—M 控制、预测、执行及数值核对。
4. {link('八项认证证明', METHOD + '/' + PROOF)}：从根误差一直到可行性、局部近最优性和条件性计算优势。
5. {link('研究现状阅读入口', FIELD + '/00_研究现状阅读说明.md')}：已有文献及审计的适用边界。
6. {link('实验与后续工作清单', EXPERIMENT + '/00_实验与后续工作清单.md')}：已做核对、待验证项和公平比较要求。

## 文件结构

| 目录 | 内容 |
|---|---|
| `01_当前方法与证明` | 当前有效正文、附录、八项证明与进度说明 |
| `02_研究现状与文献` | 领域现状、证据索引、最强先验审计、空间预瞄文献资料 |
| `03_实验与后续工作` | 从现有正文、附录和证明提取的待办及实验要求 |
| `04_规范与外部资源` | 方法设计规范、本地 E2E/ARIS 与外部资产定位 |
| `05_原文摘录与阅读材料` | 原 `tmp` 中的论文提取文本和公式核对截图 |
| `99_历史归档` | 历史交接、方法、审计、讲解与旧路线代码；保留原文件名 |
| `90_整理记录` | 逐文件原路径、新路径、大小与 SHA256；完整性验证结果 |

**目前完成的是方法正文与条件性证明文本，不是具体平台认证或实验验证。** 证明文件记录第 1—7 项已获同意继续，第 8 项待审核。

## 回溯与引用

- {link('历史版本导航', HISTORY + '/README.md')}：旧文件的角色、时间线和易混淆点。
- {link('完整文件索引与迁移关系', '90_整理记录/文件索引.md')}：按原文件名查找新位置。
- {link('迁移清单 CSV', '90_整理记录/文件迁移清单.csv')}：全部原文件的哈希与路径。
- 原文件内容保持不变；旧文中的相对路径和“当前/FINAL”字样保留其历史语义，请以迁移索引和当前入口解析。正文、附录和认证证明仍同目录，三者互相引用的文件名保持有效。
- 历史材料的指令仅供回溯；其时点、旧公式和旧验收状态不能自动覆盖当前主线。
''')

write(f'{METHOD}/00_当前进度与阅读顺序.md', f'''
# 当前进度与阅读顺序

整理于 2026-10-09；依据现有截至 2026-09-13 的文件内容。

## 当前方法

有限接触足迹提前包含前方材料信息。由位姿和力/力矩观测，经柱面 Hertz 有效接触模型和 MHE 得到当前/前方接触参数、边界、表面状态及不确定性，再构造最优控制描述符 η。

NNBO 在本项目中迁移的是 single-critic HJB 机制；实际动作、接触动力学和约束以当前正文为准，不直接照搬原论文的 force-as-action。固定基、节点和约束布局，建立有限维根 F(Z,η)=0。由 F_Z C_η=-F_η 求根灵敏度，从认证锚点预测新候选。

**候选满足误差、闭环可行性和局部性能预算即直接发布；Newton 只在有必要且条件允许时校正。** 候选 critic 恢复增益率动作，经积分形成实际阻抗；权重本身不是阻抗。

## 最近一次关键推进

早期路线是 NEOC 预测后再 Newton 校正。随后发现 cached-J chord 与预测器有相同的一阶方向，一次 Newton 也可能具有二阶根误差，不能据此宣称 NEOC 独有优势。最新路线建立预测候选的直接发布证书，把更新目标从“每次精确回根”推进为“误差预算内的可认证近最优策略”。

计算优势仍是条件性的；需要把灵敏度准备、认证、恢复和实际查询次数一起计入。

## 当前文件与状态

| 文件 | 阅读位置 | 状态 |
|---|---|---|
| {link('方法正文', MAIN)} | §1—6，尤其 §5.2—5.3 | 当前完整主线；文件名日期不代表其内容停留在第一阶段 |
| {link('技术附录', APP)} | A—M | 接触、控制、导数、误差、执行和补充数值核对 |
| {link('认证证明', PROOF)} | 第 1—8 项 | 八项文本已写出；记录第 1—7 项已获同意继续，第 8 项待审核 |

## 已有与待做

- 已有：完整数学问题与执行接口，条件性认证证明，若干合成代数/导数和低维根核对（附录 G.2、M.1—M.4）。
- 待做：第 8 项审核；指定平台全维根与认证条件实例化；域内逼近质量；带噪 MHE—控制闭环；真实组织精度与运行时间；与强基线公平比较。
- 尚不能声称：全局 HJB 最优、普遍加速、真实平台已认证、真实组织实验有效。
- “任务中心越界前完成全部调整”已不是方法成立条件；如实验实现，可作为更强结果。

## 阅读优先级

当前三份核心文件优先于历史方案和交接；原论文用于核对继承机制。发生科学冲突时应核实原始证据，不以文件名里的 FINAL 或历史验收标签解决。

历史的 `R15_CONTINUATION_HANDOFF_AFTER_STAGE1_2026-09-13.md` 仍称正文止于 §3.7、附录止于 G；当前正文已到 §6、附录已到 M，因此该交接记录是早期快照。旧“如何理解新method”和汇报版也保留在历史区，不作为最新方法摘要。

后续待办见 {link('实验与后续工作', '../' + EXPERIMENT + '/00_实验与后续工作清单.md')}。
''')

write(f'{FIELD}/00_研究现状阅读说明.md', f'''
# 研究现状与文献阅读说明

这里收集已有研究现状和文献证据，不代表截至 2026-10-09 的新检索。历史问题表述与当前方法可能不同，应结合 {link('当前方法', '../' + METHOD + '/' + MAIN)} 阅读。

| 文件 | 用途 | 使用边界 |
|---|---|---|
| {link('领域现状与研究问题', '02_FIELD_STATE_ANALYSIS_AND_FINAL_RESEARCH_PROBLEM.txt')} | §3 的已有能力和主要对照；问题演化背景 | 2026-09-08 快照，其中“最终研究问题”不是当前版本的替代 |
| {link('证据与对照索引', '05_EVIDENCE_AND_COMPARATOR_INDEX.txt')} | 核心论文、DOI、venue 标尺与最强先验的用途 | 旧索引部分文件仅在外部 E2E/ARIS 中，当前包未必包含全文 |
| {link('最强先验审计', 'R15_STRONGEST_PRIOR_FINAL_AUDIT_2026-09-07(2).txt')} | 有界搜索中的关系级对照与证据链 | 审计完成不等于当前方法或实验已验收 |
| {link('空间预瞄文献与控制融合', 'R15_CONTACT_WRENCH_PREVIEW_LITERATURE_AND_CONTROL_FUSION_2026-09-08.txt')} | 力/力矩空间观测的文献依据；§11 历史实验候选 | 文献材料仍可回查，旧 governor、deadline 和控制结构不能整体迁入当前主线 |

## 恢复研究背景时抓住什么

已有研究分别覆盖组织参数表征、可变阻抗、接触感知 MPC 和 NNBO/HJB。项目要检验的增量是：扫描中前方接触信息能否支持同一参数化控制问题的快速更新，并通过误差、可行性和性能证书落到实际阻抗。

与直接力/力矩控制、反应式变量阻抗等更简单方案的对照仍然重要。不能把单个已有模块、应用迁移或二阶预测误差本身当成已成立的新颖性。

NNBO、NEOC、隐式微分和柱面接触的当前引用见正文参考文献 [1]—[7]。论文提取文本及截图见 {link('阅读材料', '../' + SOURCE + '/README.md')}；原始 PDF 定位见 {link('外部资源说明', '../' + RESOURCE + '/README.md')}。
''')

write(f'{EXPERIMENT}/00_实验与后续工作清单.md', f'''
# 实验与后续工作清单

本清单是对现有要求的归集，不是新设计或已经执行的实验协议。平台参数、样本量、阈值和具体实验日程尚未确定，不填假数值。当前方法以 {link('正文', '../' + METHOD + '/' + MAIN)} 为准。

## 1. 已完成的核对不等于实验结果

| 已有材料 | 内容 | 来源与限制 |
|---|---|---|
| 接触代数与导数 | 合力/力矩、Jacobian、时间链式关系 | 附录 G.2、M.1；模型内部合成核对 |
| 能量与动作接口 | 能量恒等式、四维动作约束、单节点导数排列 | 附录 M.2；不是机器人闭环验证 |
| 低维有限根 | 活动/非活动分支的灵敏度、预测阶数、Newton 校正 | 附录 M.3—M.4；没有全维组织动力学和计时 |
| 八项证明 | 条件下的误差、可行性、近最优性及计算优势 | 认证证明第 1—8 项；不等于具体数值配置已通过条件 |

## 2. 后续推进顺序

以下顺序是便于执行的依赖整理；要求来自现有正文、附录与认证证明。

| 顺序 | 工作 | 必须留存的结果 | 出处 |
|---|---|---|---|
| 1 | 完成第 8 项审核，确认条件性优势的表述 | 同阶/等价边界及公平成本口径 | 认证证明 §8 |
| 2 | 指定平台与固定控制配置 | 机器人接口、周期、时域、权重、增益/速率/接触/力矩限制、尺度、critic 基与固定节点 | 附录 L.5 |
| 3 | 建立全维 HJB–KKT 根实例并核验认证条件 | 根残差、正则域、灵敏度误差、接触/活动裕量及可靠界来源 | 认证证明 §1—3；附录 I—J |
| 4 | 验证域内逼近、轨迹与性能证书 | Bellman 残差覆盖、轨迹管、约束裕量、预测与配置误差分别记账、允许性能预算 | 认证证明 §4—6 |
| 5 | 比较更新与认证计算成本 | 同一查询序列的候选质量、发布成功/失败、准备与在线延迟、恢复成本 | 正文 §5.3；附录 M.5；认证证明 §8 |
| 6 | 接入带噪 MHE 和实际增益动态 | 从材料变化到估计、根更新、发布及实际阻抗响应的全过程 | 附录 F、L、M.5 |
| 7 | 验证真实接触模型与实际平台运行 | 组织/仿体预测误差、闭环性能、真实运行时间和约束表现 | 正文范围声明；附录 G.1、L.5、M.5 |

## 3. 当前已明确的计算对照

- ANEOC：已维护灵敏度的预测候选，通过认证可直接发布。
- cached-J chord：允许复用旧 Jacobian 分解的强基线。
- one-step Newton：新参数下的一次 Newton；不自动代表所有具体 RTI 实现。
- 完整 cold solve、普通 warm solve：补充对照。

固定同一 F、基、节点、约束、时域、尺度和局部根分支，用相同发布精度、可行性和性能预算。基线获得相同传感信息、提前准备机会与结构化求解优化。不同根、不同任务质量或未完成发布不能混作同等速度。

记录 Jacobian 评估与分解、灵敏度多右端求解、预测、校正、发布认证、等待/数据搬运、通信与实际施加。后台准备计入总成本，事件后关键路径另报；记录每锚点有效查询数及失败恢复。报告平均值、分位数、最坏观测值与未完成事件。

不能仅凭 O(‖Δη‖²) 证明更快；按认证证明 §8 的关键路径和摊销条件判断。若基线缓存相同参数灵敏度运算，其与 ANEOC 的计算等价情形也必须保留。

## 4. 接触感知与闭环实验候选：保留但待适配

历史文献融合文件 §11 包含可回收的候选：同一原始日志下比较强 F-history 估计与 F+M 估计、边界/参数真值 oracle；控制对照包括固定阻抗、反应式变量增益、显式力/力矩控制。指标涉及边界/参数误差、误检漏检、接触力峰值与接触丢失、轨迹/姿态误差、图像质量、实际阻抗响应。

这些是旧路线的候选库，尚未全部确认适合当前直接认证发布方法。旧 governor、强制越界前完成、deadline 或 passivity projection 的整套消融不能照搬。适配后应另行形成正式协议。

候选来源：{link('历史文献融合 §11', '../' + FIELD + '/R15_CONTACT_WRENCH_PREVIEW_LITERATURE_AND_CONTROL_FUSION_2026-09-08.txt')}。

## 5. 解释结果时的边界

- 全维 HJB 近似精度不足，Newton 收敛也不能消除配置误差。
- 接触模型失配，不能用小拟合残差代替真实组织验证。
- 候选未通过发布证书，不能当成成功发布事件计时。
- 只有候选生成快而认证成本主导时，整体加速可能很小。
- 频繁分支失效、重锚或高维参数可能抵消灵敏度准备收益。
- 中心越界前完成调整是可能的更强结果，不是方法成立的硬截止条件。

直接来源：{link('附录 G、L、M', '../' + METHOD + '/' + APP)}；{link('认证证明第 7—8 项', '../' + METHOD + '/' + PROOF)}。
''')

write(f'{RESOURCE}/README.md', f'''
# 规范与外部资源

`补充规则` 保留四份原方法设计规范。{link('E2E/Harness 资产复用说明', '06_E2E_HARNESS_REUSE_GUIDE.txt')} 保留原路径和使用边界；其中进度与旧路线描述属于历史时点。

原始论文 PDF 不在本交接包内。历史记录列出的主要位置如下，本次整理不迁移这些外部目录，也未验证其全部文件是否仍存在：

- NNBO：`D:/桌面/科研Agent Harness设计/impedance-control-e2e/source-materials/user-batch-mpvic/Neural-Network-Based_Optimal_Impedance_Control_for_Robots_in_Physical_Interaction_With_Soft_Environments.pdf`
- NEOC：`D:/桌面/科研Agent Harness设计/impedance-control-e2e/source-materials/user-batch-topvenue-calibration/Neighboring Extremal Optimal Control Theory for Parameter-Dependent Closed-loop Laws.pdf`
- 隐式微分：`D:/桌面/科研Agent Harness设计/impedance-control-e2e/source-materials/user-batch-mpvic/Revisiting Implicit Differentiation for Learning Problems in Optimal Control.pdf`
- 其他资产：`D:/桌面/科研Agent Harness设计/ARIS`、`impedance-control-e2e`；完整历史定位见复用说明。

提取文本和截图已保留在 {link('阅读材料', '../' + SOURCE + '/README.md')}。正式公式引用仍应回到原论文。历史规范与交接中的命令不自动成为新任务指令。
''')

write(f'{SOURCE}/README.md', '''
# 原文摘录与阅读材料

原 `tmp` 目录原样迁入本目录：论文提取文本、公式页截图和多轮渲染全部保留，包括看起来重复的截图。

- [NNBO 提取文本](tmp/pdfs/nnbo.txt)
- [NEOC 提取文本](tmp/pdfs/neoc.txt)
- [隐式微分提取文本](tmp/pdfs/idoc.txt)
- 其余材料按原子目录保留在 `tmp/pdfs` 与 `tmp/pdf_review_nnbo_v2`。

历史文件里的 `tmp/...` 应解析到本目录下的 `tmp/...`。这些是辅助定位材料，不等价于出版原文全文。
''')

write(f'{HISTORY}/README.md', '''
# 历史归档导航

全部历史文件保留内容与原文件名。此处的 FINAL、CURRENT、accepted 和冻结状态描述的是各自写作时点。

| 子目录 | 内容 |
|---|---|
| `01_历史交接与状态` | 最初 00—08 交接包中的历史状态、蓝图与规则，以及 9 月 10 日、13 日的旧交接 |
| `02_历次方法与阶段稿` | 按文件名日期保留旧方法、研究问题与阶段路线 |
| `03_历次审计` | 按日期保留独立审计、必要性、验收及理论门记录 |
| `04_历史解释与汇报` | 如何理解新 method、9 月 11 日汇报解释版 |
| `05_原有历史版本` | 原 `历史版本` 目录内容，与同名其他版本分别保留 |
| `06_旧路线代码` | DS-ANEOC-NNBO 的 HC 导数核对脚本及原缓存，不作为当前柱面模型实现 |
| `07_原有空目录` | 原空阶段目录，保留其目录身份 |

## 方法演化的阅读线索

1. 9 月 7—8 日：最强先验、问题与接触预瞄审计；相关仍可用文献在 `02_研究现状与文献`。
2. 9 月 9 日：PIB/decision-focused/DS-ANEOC-NNBO 等候选路线、理论审计与收敛。
3. 9 月 10—11 日：空间预瞄、柱面接触、最优阻抗与 ANEOC 的阶段设计；汇报版包含旧 deadline 和额外结构。
4. 9 月 12 日文件名的正文与附录：实际内容已续写完整主线，见当前成果区。
5. 9 月 13 日认证证明：预测候选直接认证发布；承认与 chord/Newton 的同阶和等价边界。

这里的时间线根据文件名和内容恢复，不声称知道每次写入的精确顺序。特别是 9 月 13 日的 AFTER_STAGE1 交接仍称只完成第一阶段，已落后于当前正文，保留供回溯。

按原路径寻找文件请查 [完整迁移索引](../90_整理记录/文件索引.md)。本次未删除或去重任何历史文件。
''')

index = ['# 完整文件索引与迁移记录', '',
         '整理日期：2026-10-09。仅移动原文件，原始字节内容不改写。SHA256 见 CSV。', '',
         f'原有文件共 {len(rows)} 个；同名历史文件分别保留。旧相对引用可按本表解析。', '',
         '| 原路径 | 新位置 |', '|---|---|']
for row in rows:
    index.append(f'| `{row["原路径"]}` | {link(row["新路径"], "../" + row["新路径"])} |')
write('90_整理记录/文件索引.md', '\n'.join(index))

failures = []
for row in rows:
    p = checked(row['新路径'])
    if not p.is_file() or digest(p) != row['SHA256'] or p.stat().st_size != row['字节数']:
        failures.append(row['新路径'])
if failures:
    raise RuntimeError(f'Integrity failures: {failures}')
result = {
    '整理日期': '2026-10-09', '原文件数': len(rows), '迁移后校验通过数': len(rows),
    '原文件内容变化数': 0, '原文件删除数': 0,
    '保留的原空目录': f'{HISTORY}/07_原有空目录/R15_STAGE1_CONTACT_CLOSURE_2026-09-12',
    '说明': '原始文件逐项按字节数和 SHA256 核验；新增导航与待办仅用于整理现有内容。',
}
write('90_整理记录/完整性校验.json', json.dumps(result, ensure_ascii=False, indent=2))
print(json.dumps(result, ensure_ascii=True, indent=2))
