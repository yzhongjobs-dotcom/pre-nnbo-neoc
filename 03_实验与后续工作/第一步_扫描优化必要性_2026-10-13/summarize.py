"""Regenerate the grounded Stage-A report, diagnostics and artifact hashes."""
import hashlib,json,platform
import numpy as np
import scipy
from scipy.optimize import brentq
import scan_validation as v

P=v.HERE
def load(name):return json.loads((P/name).read_text(encoding='utf8'))
def write(name,value):(P/name).write_text(json.dumps(value,indent=2,ensure_ascii=False),encoding='utf8')

def main():
    r=load('results.json');repro=load('reproducibility.json');integration=load('integration_audit.json')
    assert repro['trajectory_exact']
    source=hashlib.sha256((P/'scan_validation.py').read_bytes()).hexdigest()
    assert source==load('frozen_config.json')['source_sha256']
    diagnosis=[]
    for E in [30000.,45000.,60000.]:
        d=brentq(lambda d:v.contact(d,E)-v.CFG['force_target'],.00005,.002)
        ke=float((v.contact(d+1e-7,E)-v.contact(d-1e-7,E))/2e-7)
        diagnosis.append(dict(E=E,force_target_indentation_m=d,contact_tangent_stiffness_N_m=ke,
            local_static_optimal_Kn_N_m=.2*(v.CFG['force_scale']/v.CFG['normal_scale'])**2/ke))
    write('quasistatic_diagnosis.json',diagnosis)
    rng=np.random.default_rng(42);energy=[]
    for i in range(40):
        t=float(rng.uniform(0,7));speed=.003;g=rng.uniform(v.LOW,v.HIGH)
        y=v.initial(speed,g)[None,:];y[0,0]=speed*t+rng.uniform(-.0004,.0004)
        y[0,1]=v.surface(speed*t)[0]+v.DNOM+rng.uniform(-.00015,.00015)
        y[0,2:4]+=rng.uniform(-.002,.002,2);target=rng.uniform(v.LOW,v.HIGH)[None,:]
        rhs,_,_,out=v.dynamics(y,t,speed,target,None,'ramp')
        fn,fr,e,ev,n,tan,delta,vt=v.geometry(y,t,speed,None,'ramp');u=out[4];vel=y[:,2:4];refv=vel-ev
        diss=np.sum(y[:,6:8]*ev**2)+v.CFG['robot_damping']*np.sum(vel**2)
        lhs=np.sum(v.CFG['mass']*vel*rhs[:,2:4])+np.sum(e*y[:,4:6]*ev)+.5*np.sum(e*e*u[:,:2])+rhs[0,8]
        fc=-fn[:,None]*n-fr[:,None]*tan
        expected=np.sum(vel*(v.CFG['force_target']*n+fc))-.7*diss-np.sum(refv*(y[:,4:6]*e+y[:,6:8]*ev))
        energy.append(float(abs(lhs-expected)))
    write('energy_audit.json',dict(samples=40,max_energy_identity_error=max(energy)))
    assert max(energy)<1e-10
    gate={}
    for name,rows in r.items():
        simple=min(['fixed','reactive','force'],key=lambda m:rows[m]['cost']);b=rows[simple];o=rows['mpc']
        cost_gain=1-o['cost']/b['cost']
        fg=1-o['force_rmse_N']/b['force_rmse_N'];tg=1-o['track_rmse_m']/b['track_rmse_m']
        valid=all(o[k]==0 for k in ['contact_loss_samples','force_violation_samples','gain_violation_samples','rate_violation_samples']) and o['min_energy_J']>=v.CFG['energy_min']
        passed=valid and ((cost_gain>=.1 and max(fg,tg)>=.1) or (o['task_tolerances_pass'] and not b['task_tolerances_pass']))
        gate[name]=dict(strongest_simple_by_cost=simple,cost_improvement_fraction=cost_gain,
            force_improvement_fraction=fg,track_improvement_fraction=tg,screen_pass=passed,
            optimized_vs_fixed_cost_fraction=o['cost']/rows['fixed']['cost']-1)
    write('screen_decision.json',gate)
    write('environment.json',dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))
    lines=['# 第一步：扫描任务的优化必要性验证','',
    '2026-10-13。**仿真已跑通、数值核对通过；本轮没有建立持续优化阻抗的必要性，未通过进入第二步的任务收益筛查。** 新增切向滑动、曲面几何和摩擦后，调优固定阻抗仍足够，显式力反馈更好。结果涉及这个任务及控制族，尚不是整个 R15 或 NEOC 原理的否定。','',
    '## 1. 做了什么','',
    '2D 质点探头沿已知平滑表面滑动；柱面 Hertz 接触使用当前材料模量、局部表面曲率和滑动摩擦/黏性。法向与切向力沿局部法线/切线耦合进入 Cartesian 动力学。实际 Kt、Kn、Dt、Dn 作为状态积分，动作受增益率限制，保留能量账本与共同保护。','',
    '扫描长度 24 mm，速度 3/6 mm/s，模量 30→60 kPa 单调或往复变化，目标力 0.14 N。表面幅度 0.5 mm、波长 48 mm。物理积分 RK4 2 ms，控制查询 100 ms，预测时域 0.4 s。材料按实际位置变化；控制器只获得当前准确 E、μ、kd，预测冻结材料，参考表面几何为各方法共有。','',
    '开发路径为慢速单调扫描；其余三条为未参与调参的验证路径。27 组固定增益候选加连续优化，18 组材料/力误差响应调度和 5 组显式力反馈系数只在开发路径调参，选择后冻结。基线均使用相同物理增益与执行限制；显式力反馈额外调整法向力命令，是检验任务必要性的替代控制架构。','',
    '有限时域在线优化采用两个阻抗目标块，经受限增益率接口执行；粗候选排序加两次 L-BFGS-B 局部求解。它是直接轨迹优化参照，尚不是原 R15 NN critic/HJB 的完整实现或全局最优证明。另做两条整段离线优化：既有目标接口，也有原始增益率动作；未来材料只用于这些诊断，不混入因果对照。','',
    '代价积分由力误差、切向/法向位置误差、速度误差、命令消耗和增益率构成；系数与配置见 [源码](scan_validation.py) 和 [冻结配置](frozen_config.json)。预设力 RMSE 容差 0.02 N、切向 RMSE 容差 0.5 mm；实质收益筛查为相对最强简单控制代价改善≥10%且力/轨迹之一改善≥10%，或新增容差通过且不新增违反。这是本轮筛查配置。','',
    '## 2. 主要结果','',
    '| 路径 | 固定代价 | 在线优化代价 | 优化相对固定 | 显式力反馈代价 | 筛查 |',
    '|---|---:|---:|---:|---:|---|']
    for name,rows in r.items():
        lines.append(f"| {name} | {rows['fixed']['cost']:.6f} | {rows['mpc']['cost']:.6f} | {100*gate[name]['optimized_vs_fixed_cost_fraction']:+.2f}% | {rows['force']['cost']:.6f} | 未通过 |")
    lines+=['','正百分比表示在线优化代价更高。四组各方法均保持接触，并通过所检查的增益/增益率、力上限和能量约束采样检查；全部满足本轮力/切向轨迹容差。','',
    '| 路径 | 固定力 RMSE / N | 优化力 RMSE / N | 力反馈 RMSE / N | 固定切向 RMSE / mm | 优化切向 RMSE / mm |',
    '|---|---:|---:|---:|---:|---:|']
    for name,rows in r.items():
        lines.append(f"| {name} | {rows['fixed']['force_rmse_N']:.6f} | {rows['mpc']['force_rmse_N']:.6f} | {rows['force']['force_rmse_N']:.6f} | {1000*rows['fixed']['track_rmse_m']:.4f} | {1000*rows['mpc']['track_rmse_m']:.4f} |")
    lines+=['','两条未来全知的离线诊断结果如下。它们检验潜在任务收益，不用于声称因果控制或预瞄优势。','',
    '| 路径 | 目标接口离线代价 | 相对固定改善 | 直接增益率离线代价 | 相对固定改善 |',
    '|---|---:|---:|---:|---:|']
    oracle=load('offline_oracle_audit.json');rate=load('direct_rate_oracle_audit.json')
    for name,a in oracle.items():
        c=a['metrics']['cost'];cr=rate[name]['metrics']['cost'];fixed=r[name]['fixed']['cost']
        lines.append(f'| {name} | {c:.6f} | {100*(1-c/fixed):.2f}% | {cr:.6f} | {100*(1-cr/fixed):.2f}% |')
    lines+=['','## 3. 根因定位','',
    '先修复求解诊断：单初值受增益率饱和平台影响，增加固定/当前增益与粗网格候选、两个初值局部求解。原单初值结果保存在 [历史诊断](00_单初值求解诊断/results.json)，新版使用冻结基线重新执行。','',
    '当前任务的主要局限有解析解释。准静态、局部线性接触中，令环境切向刚度为 kc，法向位置误差为 e，F−Fd≈−Kn e，e≈(Fd−Fref)/(kc+Kn)。忽略瞬态的力—法向位置代价给出 Kn*≈α(Fscale/nscale)²/kc。此处 α=0.2，得到 Kn*≈3.52、2.53、2.00 N/m，均低于允许下限 5 N/m，见 [诊断数值](quasistatic_diagnosis.json)。因此不同材料下都偏向低法向刚度。','',
    '切向准静态误差约为 (μFn+kd v+Br v)/Kt；当前目标与约束下，高 Kt 可直接压低扫描误差。开发路径调参确实选出 Kt=1600、Kn=5 N/m。材料缓变没有自动产生必须连续优化的内部折中：同一对边界增益就能覆盖整段扫描。此解释属于局部诊断，未替代非线性系统的全局证明。','',
    '短时域冻结材料的在线最优解，执行到真正时变环境时还可能略差于固定方案；整段全知优化只获得小幅收益，支持任务收益空间有限这一诊断。两条离线诊断也仅是有限参数化的收敛局部解，不能据此断言所有可行控制器都不可能更好。','',
    '## 4. 核对与复现','',
    f"- 独立标量根检查接触力插值：最大相对误差 {load('geometry_audit.json')['max_relative_force_table_error']:.3g}；表面相对曲率保持正值。局部抛物面接触是模型假设，尚无真实组织验证。",
    f"- 40 个非零误差/速度/增益率状态的机械—弹簧—能量账本恒等式：最大误差 {max(energy):.3g}。",
    '- 批量中心差分梯度通过独立方向核对；六个优化快照用更细预测积分和独立 SLSQP 核对，未找到进一步代价下降。所有主对照选中的局部求解均报告成功。',
    f"- 四组优化策略固定日程重放，物理步长 2→1 ms：累计代价最大相对变化 {max(a['cost_relative_difference'] for a in integration.values()):.3g}，原步长重放状态完全一致。该核对隔离积分误差，不代替重新优化。",
    '- 留出路径 ramp_fast 的完整因果优化闭环独立复跑：全部轨迹数值逐项一致；耗时重复测量，未要求耗时相等。其余三组不声称完整闭环重复验证。','',
    '## 5. 当前决定','',
    '**第一步已执行，必要性证据尚未建立，第二步 NEOC 对照未启动。** 下一次讨论应明确实际扫描中简单控制难以兼顾的物理需求，例如接触力/姿态与扫描精度之间受约束的折中，再据此修订任务或研究主张；不能仅把基线调差或放大权重来制造优势。现阶段值得保留缓变解跟踪机制，但“扫描必需持续优化阻抗”的故事仍需真实任务支持。','',
    '## 6. 文件与运行','',
    '- [主程序](scan_validation.py)：模型、调参、四组因果闭环；[配置](frozen_config.json)、[基线选择](tuning.json)、[结果](results.json)、[筛查](screen_decision.json)。',
    '- [独立核对](verify_scan.py)：求解器/积分、整段目标接口诊断与复跑；[直接增益率核对](direct_rate_audit.py)：另一种动作参数化的整段离线求解；[报告生成](summarize.py)。',
    '- 各路径 CSV 保存完整状态、接触力、误差、命令、增益率、能量与阶段代价；对应 solves.json 保存每次初值、状态、环境、解、收敛标志与时间。',
    '- [复现记录](reproducibility.json)、[独立求解核对](independent_solver_audit.json)、[积分核对](integration_audit.json)、[离线诊断](offline_oracle_audit.json)、[增益率诊断](direct_rate_oracle_audit.json)、[依赖版本](environment.json)、[文件哈希](artifact_manifest.json)。','',
    '在本目录执行：','',
    '```powershell','python -X utf8 scan_validation.py','python -X utf8 verify_scan.py','python -X utf8 direct_rate_audit.py','python -X utf8 summarize.py','```','',
    '仅复用冻结基线重跑对照时使用 `scan_validation.py --mode verify`。依赖 NumPy、SciPy；各命令顺序执行，避免核对时读取正在写入的结果。']
    (P/'README.md').write_text('\n'.join(lines)+'\n',encoding='utf8')
    entries={str(p.relative_to(P)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest()
             for p in sorted(P.rglob('*')) if p.is_file() and p.name!='artifact_manifest.json' and '__pycache__' not in p.parts}
    write('artifact_manifest.json',entries)
    print(json.dumps(dict(artifacts=len(entries),screen=gate,report=str(P/'README.md')),indent=2))

if __name__=='__main__':main()
