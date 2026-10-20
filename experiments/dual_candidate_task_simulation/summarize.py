"""Report shared physical errors separately from gain-preference OCP costs."""
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent

def load(name):
    data=json.loads((ROOT/(name+'.json')).read_text());arrays=np.load(ROOT/(name+'_traces.npz'))
    rows=data['results'];scores=[]
    for r in rows:
        a=arrays[str(r['case'])];dt=float(a[1,0]-a[0,0])
        stage=(a[:,1]/.002)**2+(a[:,2]/.001)**2+.04*np.sum((a[:,3:5]/.015)**2,axis=1)+((a[:,10]-.6)/.1)**2
        terminal=.1*((a[-1,1]/.002)**2+(a[-1,2]/.001)**2+((a[-1,10]-.6)/.1)**2)
        scores.append(float(dt*stage.sum()+terminal))
    return data,rows,np.asarray(scores)

def stats(name):
    data,rows,score=load(name)
    keys=['force_rmse_N','tangent_rmse_mm','normal_rmse_mm','cost','decision_median_ms','decision_p95_ms','mhe_total_seconds']
    out={k:float(np.mean([a[k] for a in rows])) for k in keys}
    out['common_sampled_physical_score']=float(score.mean());out['n_cases']=len(rows)
    out['violating_samples']=sum(a['constraint_violating_samples'] for a in rows)
    out['skipped_updates']=sum(a.get('skipped_busy_updates',0) for a in rows)
    out['expired_policy_queries']=sum(a.get('expired_policy_queries',0) for a in rows)
    out['case_ids']=[a['case'] for a in rows]
    out['strata']={str(j):{k:float(np.mean([r[k] for r in rows if r['kind']==j])) for k in keys[:4]} for j in range(4)}
    return out

def paired(a,b,metric='cost'):
    da,ra,sa=load(a);db,rb,sb=load(b)
    assert [x['case'] for x in ra]==[x['case'] for x in rb]
    assert da['cases']==db['cases']
    va=sa if metric=='physical' else np.array([r[metric] for r in ra]);vb=sb if metric=='physical' else np.array([r[metric] for r in rb])
    d=100*(vb-va)/np.maximum(vb,1e-10);rng=np.random.default_rng(16901)
    boot=np.mean(d[rng.integers(0,len(d),(10000,len(d)))],axis=1)
    return {'metric':metric,'mean_paired_improvement_percent':float(np.mean(d)),'paired_bootstrap_95_percent':np.percentile(boot,[2.5,97.5]).tolist(),
            'ratio_of_means_improvement_percent':float(100*(vb.mean()-va.mean())/vb.mean()),'n':len(d)}

def main():
    names=['candidate1_test','candidate1_without_neighbor_test','candidate1_full_value_reference_test',
           'candidate2_test_seed16111','candidate2_test_seed16112','candidate2_test_seed16113','jin_ms_test',
           'candidate1_preview_audit_full','candidate1_preview_audit_removed','candidate2_preview_audit_full','candidate2_preview_audit_removed',
           'candidate1_shifted_test','candidate2_shifted_test','candidate1_without_moment_test',
           'compiled_timed_candidate1','compiled_timed_candidate1_without_neighbor','compiled_timed_full_value','compiled_timed_candidate2','compiled_timed_jin_ms']
    out={'scope':'Synthetic planar model. Primary nominal results use zero computation wait. Timing rows are a separate isolated 8-case publication experiment. No real robot or journal acceptance claim.',
         'shared_metric':'sampled physical tracking/velocity/force plus endpoint term; no gain preference or gain-rate penalty',
         'confidence_scope':'paired case bootstrap for fixed fitted methods; no simultaneous-inference correction; neural seed summaries are reported separately',
         'methods':{name:stats(name) for name in names},'comparisons':{}}
    for name,a,b,metric in [
        ('neighbor_effect','candidate1_test','candidate1_without_neighbor_test','cost'),
        ('value_quality_reference','candidate1_test','candidate1_full_value_reference_test','cost'),
        ('preview_candidate1','candidate1_preview_audit_full','candidate1_preview_audit_removed','cost'),
        ('preview_candidate2_seed16111','candidate2_preview_audit_full','candidate2_preview_audit_removed','cost'),
        ('candidate1_vs_jin_shared_physical','candidate1_test','jin_ms_test','physical'),
        ('candidate2_vs_jin_shared_physical','candidate2_test_seed16111','jin_ms_test','physical'),
        ('timed_candidate1_vs_jin','compiled_timed_candidate1','compiled_timed_jin_ms','physical'),
        ('timed_candidate2_vs_jin','compiled_timed_candidate2','compiled_timed_jin_ms','physical')]:out['comparisons'][name]=paired(a,b,metric)
    seed_methods=[out['methods'][f'candidate2_test_seed{s}'] for s in [16111,16112,16113]]
    out['withdrawn_ablation']={'files':['candidate1_without_preview_test.json','candidate2_without_preview_test.json'],
        'reason':'old-region-only masking also removed updated current material after crossing; not a valid spatial-preview ablation',
        'correction':'uniform current-contact equivalent retains estimated current Fn/Ft; separate no-preview actor retrained; new seed16501 cohort n24',
        'old_improvement_numbers_not_used':[26.86276488049843,6.172869256143549]}
    out['candidate2_seed_mean_and_range']={k:{'mean':float(np.mean([s[k] for s in seed_methods])),'range':[min(s[k] for s in seed_methods),max(s[k] for s in seed_methods)]}
        for k in ['force_rmse_N','tangent_rmse_mm','normal_rmse_mm','cost','common_sampled_physical_score']}
    (ROOT/'summary.json').write_text(json.dumps(out,ensure_ascii=False,indent=2))
    print(json.dumps({'comparisons':out['comparisons'],'neural_seeds':out['candidate2_seed_mean_and_range']},indent=2))
    lines=['# 两候选各自仿真与比较结果','','24个物理场景，每类6例。三个学习种子复用同一批场景，不能称72个独立物理试验。零计算等待主表如下。',
           '','| 方法 | 力RMSE/N | 切向RMSE/mm | 法向RMSE/mm | 共同物理评分（越低越好） | 采样约束违反 |','|---|---:|---:|---:|---:|---:|']
    for label,key in [('候选1：局部价值反馈＋邻近更新','candidate1_test'),('候选1关闭邻近更新','candidate1_without_neighbor_test'),('较充分同族价值求解','candidate1_full_value_reference_test'),('候选2：种子16111','candidate2_test_seed16111'),('候选2：种子16112','candidate2_test_seed16112'),('候选2：种子16113','candidate2_test_seed16113'),('Jin多步方法扫描适配','jin_ms_test')]:
        s=out['methods'][key];lines.append(f"| {label} | {s['force_rmse_N']:.6f} | {s['tangent_rmse_mm']:.4f} | {s['normal_rmse_mm']:.4f} | {s['common_sampled_physical_score']:.5f} | {s['violating_samples']} |")
    lines+=['','这里的共同物理评分只包含运动、速度和法向力，按相同采样积分；不含某控制器族的增益偏好/率正则。不能用Jin在本候选完整OCP代价中的大数值声称胜出。单项误差显示法向位置与力/切向运动的取舍。',
            '','| 实测等待发布：8场景 | 力RMSE/N | 切向RMSE/mm | 法向RMSE/mm | 每场景决策中位耗时的均值/ms | 跳过的更新 |','|---|---:|---:|---:|---:|---:|']
    for label,key in [('候选1','compiled_timed_candidate1'),('关闭邻近预测','compiled_timed_candidate1_without_neighbor'),('较充分价值求解','compiled_timed_full_value'),('候选2','compiled_timed_candidate2'),('Jin扫描适配','compiled_timed_jin_ms')]:
        s=out['methods'][key];lines.append(f"| {label} | {s['force_rmse_N']:.6f} | {s['tangent_rmse_mm']:.4f} | {s['normal_rmse_mm']:.4f} | {s['decision_median_ms']:.3f} | {s['skipped_updates']} |")
    lines+=['','计时是所有方法共用编译标量接触核后的CPU原型，2ms发布网格，非机器人实时认证。旧timed_*文件是NumPy逐点版本诊断，不作为当前速度结论。模型更新时间与网络逐次推理的粒度不同，不直接拿两个中位时间作全链加速比；估计开销另计且已进入发布延迟。',
            '','机制、分层数据、置信区间、分布变化及训练种子范围见[summary.json](summary.json)。小区间只反映固定实现的场景差异，不是顶刊达标结论。']
    lines+=['','逐例分项与示例轨迹：[PNG图](results.png)、[PDF图](results.pdf)。']
    def preview_text(key):
        d=out['comparisons'][key];ci=d['paired_bootstrap_95_percent']
        return f"{d['mean_paired_improvement_percent']:.3f}%，配对95%区间[{ci[0]:.3f}%,{ci[1]:.3f}%]"
    lines+=['','## 这些结果具体说明什么','',
        '预瞄消融经核查重做：旧版本把新材料进入当前区后的信息也删掉了，原26.86%/6.17%不作为贡献证据。现保留与当前接触Fn/Ft匹配的均匀模型，只去掉空间变化；学习版重新独立训练。在新seed16501的24场景比较，未用原测试挑参数。',
        '该对照匹配当前Fn/Ft，不匹配全部力矩和局部导数；因此检验空间参数化信息的整体作用，不能把收益单独归因于时间提前量。',
        '',
        '候选1的修正预瞄消融：完整任务代价逐例平均改善'+preview_text('preview_candidate1')+'。NEOC开/关仍只改善约0.007%，区间[-0.260%,0.207%]；当前带邻近预测还比同预算关闭它更慢。信息作用与邻近更新价值必须分别判断。较充分同族求解也不是已知真值的全局最优控制，不能将微小实测胜幅解释成超过最优。',
        '',
        '候选2的修正预瞄消融（同训练种子16111、各自训练）：逐例完整代价改善'+preview_text('preview_candidate2_seed16111')+'；这是固定种子的条件结果。三个主种子的共同物理评分均值为0.06029，Jin适配为0.05646，约高6.78%。不能拿单次网络推理快或一个种子的消融当作总体优于领域方法。',
        '',
        '两候选法向位移误差略小，力与切向运动误差仍更大；这是真实取舍。目前没有顶刊级方法优势结论。变软和更快表面8场景中，候选1/2力RMSE均值约0.01792/0.02482N；这个额外队列没有外部基线，不作为超越已有工作的证据。约束检查均为当前采样轨迹，主要场景能量约束不活跃，不能据此声称证明了主动约束场景或真实接触安全。',
        '',
        '下一轮应解释并改进：候选1为何额外预测/校正未换来比普通反馈更新更好的成本；候选2为何利用了预瞄仍未达到强方法的力/切向效果，并检查训练目标与实际任务取舍。可沿已有强阻抗结构改进或学习修正，但必须明确改什么并用新场景验证；不通过改分数、虚构计算瓶颈或新增模块数量宣布成功。实际平台/标定仍待核实。']
    (ROOT/'results.md').write_text('\n'.join(lines)+'\n')

if __name__=='__main__':main()
