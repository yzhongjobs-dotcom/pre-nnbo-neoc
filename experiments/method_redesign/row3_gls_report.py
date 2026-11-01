"""Render only completed row3_gls JSON summaries; no fitting or selection."""
import json
from pathlib import Path


def f(x):
    return '—' if x is None else f'{x:.3f}'


def pct(x):
    return '—' if x is None else f'{100*x:.1f}%'


def markdown(r):
    cfg = r['selected']['params']
    lines = ['# E4：刚度 GLS 漂移先验重新表述', '',
        '## 交付与评测边界', '',
        '主目录拒绝写入，本次在可写副本运行。补丁仅新增 experiments/method_redesign/ 下四个 Python 文件和本报告/JSON；没有改动原文件，没有提交或推送。既有 frozen snapshot 与 130000+ final 确认集均未访问。',
        '调参种子190000+1000×族序号、每族10场景；12组事先固定候选（与第2步E1/E2/E3相同评估次数、相同轨迹及目标），目标为逐场景 median(|ln(k̂/k)|)+0.25 P95，先场景再族等权。长度尺度∈{0.008,0.02,0.05,0.12} s，力标准差比例∈{0.03,0.3,3}。未按覆盖率、危险方向或闭环误差选参。E0固定，E2直接用第2步选定参数，不重调。',
        '开环：同一最终Proposed轨迹，诊断集110000+、每族20场景。闭环：同一诊断集每族50场景，两种方法共享场景和传感器噪声种子；sweep沿用两组各6档、每档20个场景（160000+推拉斜坡时长、170000+刚度比），每档两种方法。已用诊断集上的先验失效分析构建本候选，因此这些结果是开发诊断，不能改称独立确认。',
        '只用当前及过去的法向测力、编码器和已发指令；E4没有使用Ft/My、表面真值或刚度真值。评分才使用真值。开环清空真值列并逐拍截断；闭环继承FINAL，只覆盖刚度估计与量测通道接受条件。', '',
        '## 模型、消元与方差', '',
        r'窗口模型为 $z=kY+B\beta+r+\varepsilon$。$B$ 为 common._build 生成的、与Y和z同处理域的三次时间基；$\beta$ 取平坦先验。令 $Q^TQ=I,\ Q^TB=0$，$y=Q^TY,\ u=Q^Tz$。力域核 $K_{ij}=v\exp[-(t_i-t_j)^2/(2\ell^2)]$，$v=(a F_w)^2$，$F_w=\max(0.5,\operatorname{median}|F_j|)$ N，使用该窗口保留区间的测得端点力。这里v是方差，a是标准差比例；超参实际报告为a及a²。',
        r'$C=\sigma_\varepsilon^2 I+Q^TKQ$，$\hat k=(y^TC^{-1}u)/(y^TC^{-1}y)$，$V_k=(y^TC^{-1}y)^{-1}$。默认 $\sigma_\varepsilon=0.04$ N，窗口16拍×8 ms，三次基与编码器插值沿用FINAL。长度尺度以秒计，核作用于控制区间方程误差的力域，不把物理表面位移GP直接当成力GP。',
        r'令 $P=QQ^T$，则原坐标写法为 $Y_\perp^T\Omega^{-1}z_\perp/(Y_\perp^T\Omega^{-1}Y_\perp)$，$\Omega=\sigma^2I+PKP$。P和PKP秩亏，但正的σ²I使Ω本身正定；只在正交补上计算最直接，也避免对Pσ²IP的秩亏逆。完整联合GLS对β做Schur消元与此式相同。',
        r'固定Y、B和C时，$Q^T(z+B\gamma)=Q^Tz$，故任何三次漂移精确消失；v→0恢复等权投影。代码用完整QR正交补与投影核特征分解，按相同保留区间模式缓存几何，每拍随力尺度更新谱权重。缓存采用8 ms均匀网格；累计时间的浮点误差仅为舍入量。原实现的1e−12 Gram正则与这里的QR精确投影存在微小数值差异，验证中最大刚度差见JSON。',
        '重要边界：力尺度Fw是经验贝叶斯式的可观测协变量，代数不变性指给定该尺度/协方差后对Bβ的消去；如果人为改变原始测力同时重算Fw，谱权重也可能改变，不能声称整个数据依赖算法对此变换无条件不变。既没有给多项式系数加收缩先验，也没有把历次重叠窗口当独立数据递归累积。',
        r'方差公式以给定Y、Fw时，r+ε条件零均值且协方差正好为C为前提。确定性残余漂移下，条件偏差仍是 $(y^TC^{-1}Q^Tr)/(y^TC^{-1}y)$；GLS重新加权而不保证消除这个偏差。闭环Y内生、非平稳运动、滤波后相关噪声、窗口内刚度变化、估计力尺度和模型失配都不由该名义方差完整覆盖。', '',
        '## 与 Rossi 带通的谱加权关系（推导，不是原算法等价复现）', '',
        r'若均匀采样、忽略有限窗边界，并把残差协方差近似为平稳循环矩阵，则傅里叶基将其对角化，$w(\omega)=1/[\sigma^2+S_r(\omega)]$，$\hat k=\sum_\omega w\,\Re(\overline{\tilde Y}\tilde z)/\sum_\omega w|\tilde Y|^2$。平方指数核的谱满足 $S_r(\omega)\propto v\ell\exp(-\ell^2\omega^2/2)$，故较低频的残余漂移受到较小权重，硬带通对应通带有限权、阻带协方差趋于无穷（或将权置零）的极限。',
        r'在同一滤波器H作用于Y和z后做最小二乘，其谱权为$|H|^2$；选择$|H|^2\propto1/(\sigma^2+S_r)$即得到这种谱GLS。项目Rossi型“同频带通后回归”在此意义下可视为谱加权特例。Rossi的EWRLS遗忘、滤波瞬态、有限窗边界，以及三次投影与平稳核一般不对易，均不由上述理想化推导保证等价。SE+白噪声自然给出软高通倾向，不等于特定20.8 Hz窄带滤波器；还需带外噪声假设才会得到真正带通权。', '',
        '## 接受规则与闭环接口', '',
        r'预先固定 $2\sigma_{GLS}/\hat k\le\ln1.2$，只接受正且未越过既定刚度范围[200,300000] N/m的原始估计。delta法给出 $\operatorname{sd}(\ln k)\simeq\sigma/k$，因而名义约95.45%的±2σ对数区间对应$k/\hat k\in[1/1.2,1.2]\subset[0.75,1.2]$。即使直接用线性高斯±2σ区间，它落在[1−ln1.2,1+ln1.2]≈[0.818,1.182]，也在既有冻结分析容许区间内。对数区间解释是一阶近似，冻结区间和GP后验均不构成实际变增益闭环保证。',
        '该规则替换原b风险评分；FINAL原有precision与consistency前提保留。consistency改用GLS白化残差χ²、自由度为有效区间数−rank(B)−1，阈值仍为FINAL的2.6797；不在诊断集调任何阈值。反馈仍按FINAL的相对精度门控、cσ上裕量、rho_d下降限速、无效窗口仅上调与初始化规则，只是k̂、σ、χ²由E4给出；反馈不另套更严的量测门。IMM、激励、前馈和接触监督不变。',
        '开环E4全时段报告与E0同样的范围截断估计（启动/无效窗口保持）；accepted行是通过上述量测门的当拍子集，并非保持后的量测增益误差。实际保持与选择反馈的任务影响由闭环评价。方差门可能因先验大而长期不更新，或者因错配而自信地接受偏差窗口；须结合接受率、覆盖率和闭环结果判断。', '',
        '## 选定超参', '',
        f'ℓ={cfg[0]} s；a={cfg[1]}；方差比例a²={cfg[1]**2:g}；v=(a Fw)² N²。调参目标={r["selected"]["objective"]:.6f}。', '',
        '| ℓ s | a | 调参目标 |', '|---:|---:|---:|']
    for row in r['tuning']:
        lines.append(f'| {row["params"][0]} | {row["params"][1]} | {row["objective"]:.6f} |')
    if 'prior_example' in r:
        a=r['prior_example']
        lines += ['',f'幅度解释：在满16区间、测力{a["force_N"]:g} N的规则窗口上，先验经三次投影后的正交补均方标准差为{a["projected_residual_rms_N"]:.5f} N（不是3×4=12 N）；大部分长尺度GP分量已被三次基消去。投影核谱见JSON prior_example。该例只描述先验，不是真实运动误差测量。']
    lines += ['', '## 开环：同轨迹估计', '',
        '误差为|ln(k̂/k)|；危险偏低为k̂/k<1/1.2，危险偏高为k̂/k>1/0.75（与第2步一致，并非诊断中更紧的0.83界）。1σ/2σ指真值落入k̂±σ/±2σ的时间比例。全部统计从1 s起、k真值>0；不因失联或估计器状态筛掉困难样本。表为有效拍汇总，不把拍当独立重复；JSON另保留逐场景摘要。', '',
        '| 族 | 方法 | 中位 | P95 | 偏低 | 偏高 | 1σ覆盖 | 2σ覆盖 | P50/P99 ms |',
        '|---|---|---:|---:|---:|---:|---:|---:|---:|']
    for fam, methods in r['open_loop']['summary'].items():
        for method, values in methods.items():
            a = values['all']
            lines.append('| '+' | '.join([fam,method,f(a['median']),f(a['p95']),pct(a['under']),pct(a['over']),pct(a['cover1']),pct(a['cover2']),f(a['p50_ms'])+'/'+f(a['p99_ms'])])+' |')
    for segment, title in [('transition','过渡段'),('ramp','推拉斜坡段'),('accepted','E4通过量测门的当拍子集')]:
        lines += ['', '## '+title, '',
            '过渡段定义沿用第2步：|d ln k/dx|>10 m⁻¹；也含运动引起的工作点变化。推拉斜坡按场景事件的进入和退出时间区间标记，只用于评分，不传入估计器；不是整个推拉事件或保持段。', '',
            '| 族 | 方法 | 拍数 | 中位 | P95 | 偏低 | 偏高 | 1σ覆盖 | 2σ覆盖 |',
            '|---|---|---:|---:|---:|---:|---:|---:|---:|']
        for fam, methods in r['open_loop']['summary'].items():
            for method, values in methods.items():
                if segment not in values: continue
                a = values[segment]
                lines.append('| '+' | '.join([fam,method,str(a['ticks']),f(a['median']),f(a['p95']),pct(a['under']),pct(a['over']),pct(a['cover1']),pct(a['cover2'])])+' |')
    lines += ['', '| 族 | E4接受拍比例（共同评分域） |', '|---|---:|']
    for fam, methods in r['open_loop']['summary'].items():
        lines.append(f'| {fam} | {pct(methods["E4"]["acceptance_rate"])} |')
    lines += ['', '## 闭环：诊断集每族50', '',
        'RMSE单位N；失控为真实力峰值>50 N；长失联为最长单次失联>0.5 s。失败单列，RMSE与接触指标仅统计未发散场景，配对仅含两方法均成功的场景。bootstrap为按场景成对重采样10000次、随机种子291101，对均值差E4−Proposed给百分位95%区间；逐项探索性区间，未作多重比较校正。', '',
        '| 族 | Proposed中位/均值 | E4中位/均值 | 失控P/E | 长失联P/E | 失败P/E | 配对均值差 [95%区间] |',
        '|---|---:|---:|---:|---:|---:|---:|']
    def closed_line(name, a):
        p,e,b = a['Proposed'],a['E4'],a['paired']
        return f'| {name} | {p["median"]:.3f}/{p["mean"]:.3f} | {e["median"]:.3f}/{e["mean"]:.3f} | {p["runaway"]}/{e["runaway"]} | {p["long_loss"]}/{e["long_loss"]} | {p["failed"]}/{e["failed"]} | {b["mean_delta"]:+.3f} [{b["ci95"][0]:+.3f}, {b["ci95"][1]:+.3f}] (n={b["n"]}) |'
    for fam in S_FAMS:
        lines.append(closed_line(fam,r['diagnostic']['summary'][fam]))
    lines += ['', '| 族 | E4整段从未接受的场景 | 接受/有效拟合拍（全时段） |', '|---|---:|---:|']
    for fam in S_FAMS:
        rr=[a for a in r['diagnostic']['episodes'] if a['group']==fam and a['method']=='E4']
        n0=sum(a['accepted']==0 for a in rr)
        rate=sum(a['accepted'] for a in rr)/max(sum(a['attempted'] for a in rr),1)
        lines.append(f'| {fam} | {n0}/{len(rr)} | {pct(rate)} |')
    lines += ['', '## Sweep：两组各6档、每档20', '',
        '| 剂量 | Proposed中位/均值 | E4中位/均值 | 失控P/E | 长失联P/E | 失败P/E | 配对均值差 [95%区间] |',
        '|---|---:|---:|---:|---:|---:|---:|']
    for kind, doses in [('ramp',(.6,.45,.3,.2,.15,.1)),('ratio',(3.,10.,20.,30.,45.,60.))]:
        for dose in doses:
            key = f'{kind}:{dose}'
            lines.append(closed_line(key,r['sweep']['summary'][key]))
    lines += ['', '## 计算时间', '',
        '开环表为每拍设计矩阵、QR/谱缓存访问、标量拟合和方差输出时间；不含场景生成与评分。闭环表为各场景控制器总调用时间P50/P99再取场景中位，剔除首1 s，含刚度、IMM、导纳与监督，排除对象积分与文件写入；二者不是同一计时口径。4线程计时包含GIL与CPU竞争，不是8 ms硬实时认证。', '',
        '| 族 | Proposed总调用P50/P99 ms | E4总调用P50/P99 ms |', '|---|---:|---:|']
    for fam in S_FAMS:
        a=r['diagnostic']['summary'][fam]
        p,e=a['Proposed'],a['E4']
        lines.append(f'| {fam} | {p["p50_ms"]:.3f}/{p["p99_ms"]:.3f} | {e["p50_ms"]:.3f}/{e["p99_ms"]:.3f} |')
    if 'serial_timing' in r:
        lines += ['', '另按预先指定的每族首个诊断场景串行复跑，避免把本任务线程竞争误认为算法代价；共15个开环重放、10个闭环运行，全部科学指标与主评测逐项完全一致。未控制机器上的外部负载；这是每族1场景的时间样本，不是50场景的总计时结论。', '',
            '| 族 | E0估计P50/P99 ms | E4估计P50/P99 ms | Proposed总调用P50/P99 ms | E4总调用P50/P99 ms |',
            '|---|---:|---:|---:|---:|']
        for fam,a in r['serial_timing']['summary'].items():
            fields=[]
            for key in ('E0_estimator','E4_estimator','Proposed_controller','E4_controller'):
                fields.append(f'{a[key]["p50_ms"]:.3f}/{a[key]["p99_ms"]:.3f}')
            lines.append('| '+fam+' | '+' | '.join(fields)+' |')
    lines += ['', '## 结果判断与局限', '']
    for fam in S_FAMS:
        o=r['open_loop']['summary'][fam]
        p,e=o['E0']['all'],o['E4']['all']
        b=r['diagnostic']['summary'][fam]['paired']
        direction='更低' if b['ci95'][1]<0 else ('更高' if b['ci95'][0]>0 else '区间跨0')
        lines.append(f'{fam}：开环中位/P95由{p["median"]:.3f}/{p["p95"]:.3f}变为{e["median"]:.3f}/{e["p95"]:.3f}；闭环均值差{b["mean_delta"]:+.3f} N，95%区间{direction}。')
    favorable = [fam for fam in S_FAMS if r['diagnostic']['summary'][fam]['paired']['ci95'][1]<0]
    adverse = [fam for fam in S_FAMS if r['diagnostic']['summary'][fam]['paired']['ci95'][0]>0]
    lines += ['', f'诊断族中，配对均值区间完全低于0：{favorable}；完全高于0：{adverse}。这种逐族比较不自动证明族间差异。',
        f'整体上，E4的2σ覆盖为{pct(r["open_loop"]["summary"]["overall"]["E4"]["all"]["cover2"])}，但通过量测门的子集仅为{pct(r["open_loop"]["summary"]["overall"]["E4"]["accepted"]["cover2"])}；后者的危险偏低比例仍为{pct(r["open_loop"]["summary"]["overall"]["E4"]["accepted"]["under"])}。这直接反驳将名义后验门当作真实95%可靠性保证。软组织推拉在固定Proposed轨迹上的接受率为0，说明更严格门会损失量测通道的可用性。',
        '是否替换：本轮不建议直接替换既有投影。GLS提供了自洽的条件模型与谱权重，但必须同时满足开环误差、方差覆盖、更新可用性、闭环尾部和算力的要求；本轮是12候选、诊断集上的一个实例，不能用部分族的中位数改善替代完整证据。具体正负方向以上述逐族与剂量表为准。',
        '解释边界：GLS降低先验认为漂移占优方向的权重，同时损失辨识信息并改变接受窗口；它能缓解低频泄漏，却无法识别与Y完全共线的残余运动，也没有建模时变k。闭环收益或恶化包含点估计、方差门、保持旧量测增益和保守反馈刚度的共同作用，本次没有以额外消融把它们分别归因。用GLS白化残差一致性替换白噪声一致性也可能放宽部分失配窗口。开环0接受率不意味着E4闭环从未接受：闭环反馈、激励幅值和信噪比随估计改变，分布已不同，应以上表整段接受统计为准。',
        '主要局限：GP是工程先验而非被独立验证的真实运动分布；SE核很平滑，与推拉分段运动及非平稳接触可能失配；力域尺度为经验选择，未包含其超参不确定度；窗内k变化和力滤波相关噪声没有完整建模；线性冻结容许带不证明时变/切换稳定。诊断集与已查看过的sweep都不能充当新确认集；没有机器人通信或真机证据。', '',
        '## 验证与复现', '',
        '已通过：任意三次漂移不变性（给定协方差）、v=0退化、正交补与稠密投影求解一致、联合GLS与方差一致、2万组合成高斯样本方差检查、真值列污染/未来前缀不变；E0/E2诊断统计逐项复现第2步（详见JSON verification）。原文件SHA256见JSON，补丁只含新增代码与汇总，不含轨迹缓存或二进制。', '',
        '`python -X utf8 experiments/method_redesign/row3_gls.py --workers 4`',
        '`python -X utf8 experiments/method_redesign/row3_gls.py --verify-only`',
        '`python -X utf8 experiments/method_redesign/row3_gls.py --workers 4 --resume`',
        '`python -X utf8 experiments/method_redesign/audit_row3_gls.py`',
        '第一条从头调参、生成轨迹和全部评价；第二条做有科学意义的代数/泄漏检查；第三条只复用同名JSON中已完成阶段，不重新选择配置。重跑计时会变化。', '']
    return '\n'.join(lines)


S_FAMS = ('S2','T3m','T2m-soft','T2m-stiff','T2f-stiff')


if __name__ == '__main__':
    path=Path(__file__).resolve().parent/'results/row3_gls.json'
    result=json.loads(path.read_text(encoding='utf-8'))
    assert result['status']=='complete'
    path.with_suffix('.md').write_text(markdown(result), encoding='utf-8')
