"""Independent algebra and observation-interface checks for the contact model pilot.

No simulation run: reuse measured diagnosis inputs. This also quantifies two
identifiability issues instead of mistaking regularised Hessians for evidence.
"""
import json
import math
from pathlib import Path
import numpy as np

from finite_patch_model import FinitePatch, Geometry, orthogonal_residual
from line_patch_model import LinePatch
from run_finite_patch_pilot import DepthFit, observed_window, verify

HERE = Path(__file__).resolve().parent


def main():
    report = {'physics': verify()}
    lm = LinePatch(spatial=True); th = np.array([math.log(3e4), math.log(3), .001, .013])
    xx = np.linspace(-.015, .025, 7); dd = np.linspace(.0003, .003, 7)
    ff, _, kk = lm.response(xx, dd, th); h = 1e-7
    fd = (lm.response(xx, dd+h, th)[0]-lm.response(xx, dd-h, th)[0])/(2*h)
    report['line_physics'] = dict(derivative_relative_error=float(np.max(np.abs(fd-kk)/kk)),
                                 inverse_relative_error=float(np.max(np.abs(lm.inverse(xx, ff, th)-dd)/dd)))
    assert report['line_physics']['derivative_relative_error'] < .03
    assert report['line_physics']['inverse_relative_error'] < 1e-5
    a = np.load(HERE/'results/finite_patch_input_S2_110000.npz')
    z = a['measured'].copy(); end = 480
    obs = observed_window(z, end, 16, tau=a['baseline_tau'][end-1])
    zz = z.copy(); zz[end:] = 1e6
    other = observed_window(zz, end, 16, tau=a['baseline_tau'][end-1])
    for key in obs:
        np.testing.assert_array_equal(obs[key], other[key])
    report['causal_prefix'] = True

    # Same candidate, same causal observations, but realistic systematic wrench
    # offsets and force-to-moment coupling. These are not white-noise proxies.
    corrupted = z.copy(); corrupted[:, 6] += .035+.005*z[:, 5]+.003*z[:, 3]
    obs2 = observed_window(corrupted, end, 16, tau=a['baseline_tau'][end-1])
    N = np.c_[np.ones(len(obs['f'])), obs['ft'], obs['f']]
    nuisance_error = float(np.max(np.abs(orthogonal_residual(N, obs['my'])-
                                         orthogonal_residual(N, obs2['my']))))
    assert nuisance_error < 1e-10
    k0 = float(a['baseline_kraw'][end-1])
    r1 = DepthFit(spatial=True, moment=True).fit(obs, k0)
    r2 = DepthFit(spatial=True, moment=True).fit(obs2, k0)
    change = abs(math.log(r1['k']/r2['k']))
    assert change < .005, change
    report['bias_and_lever'] = dict(bias_Nm=.035, tangent_lever_m=.005, normal_lever_m=.003,
                                   projected_difference_Nm=nuisance_error, fitted_logk_change=change,
                                   limitation='constant bias and lever only; drift, scale/cross-axis errors untested')

    # Uniform p=1 must reproduce independently computed FORWARD regression.
    yp, fp = (orthogonal_residual(obs['B'], obs[k]) for k in ('y', 'f'))
    direct = float(yp@fp/(yp@yp))
    fit = DepthFit().fit(obs, k0)
    relative = abs(fit['k']/direct-1)
    assert relative < 2e-4, relative
    report['uniform_forward_limit'] = dict(independent_k=direct, fitted_k=fit['k'], relative_error=relative)

    # At a single force and tangent, infinitely many power-law contacts agree.
    F, K = 4., 1e4
    alternatives = []
    for p in (1., 1.2, 1.5):
        d = p*F/K; C = F/d**p
        alternatives.append(dict(p=p, depth=d, coefficient=C, force=C*d**p, stiffness=p*C*d**(p-1)))
    report['operating_point_ambiguity'] = alternatives
    # For y=a sin(wt), second-harmonic force amplitude from curvature is F'' a²/4.
    # q1 is MEASURED first-harmonic force ripple, not the command-domain amplitude.
    report['curvature_second_harmonic_N'] = {
        str(q1): {str(p): (p-1)/(4*p)*q1*q1/F for p in (1., 1.5)}
        for q1 in (.03, .1, .167)}

    # Nuisance parameter rank: at constant load and sliding speed, Ft=-mu F.
    N0 = np.c_[np.ones(16), np.ones(16)*(-1.6)]
    report['constant_sliding_wrench_nuisance_rank'] = int(np.linalg.matrix_rank(N0))
    assert report['constant_sliding_wrench_nuisance_rank'] == 1
    results = HERE/'results/finite_patch_ablation_v4.json'
    if results.exists():
        report['task_data_validation'] = score_results(results)
    # A local implementation failure is retested with a larger solver budget,
    # holding the old v3 model/interface fixed. Selection uses solver status,
    # never the reference stiffness error.
    old = json.loads((HERE/'results/finite_patch_continuous_v3.json').read_text(encoding='utf-8'))
    budget = []
    for fam in ('S2', 'T3m'):
        r = next(v for v in old['rows'] if v.get('family')==fam and v.get('model')=='patch_wrench'
                 and 'k' in v and not v['solver_success'])
        a = np.load(HERE/'results'/f'finite_patch_input_{fam}_{r["sid"]}.npz')
        obs = observed_window(a['measured'], r['end'], 16, tau=a['baseline_tau'][r['end']-1])
        retry = DepthFit(spatial=True, moment=True, wrench_nuisance='bias_tangent', max_nfev=100).fit(obs, r['baseline_kraw'])
        budget.append(dict(family=fam, end=r['end'], original_k=r['k'], retry_k=retry['k'],
                           original_residual=r['residual_rms'], retry_residual=retry['residual_rms'],
                           retry_success=retry['solver_success'], retry_rank=retry['data_rank']))
    report['solver_budget_repair'] = budget
    # Geometric and sensor calibration sensitivity on a fixed, causal window.
    # This is a mechanism test, not a trajectory robustness claim.
    a = np.load(HERE/'results/finite_patch_input_S2_110000.npz'); end = 480
    original = observed_window(a['measured'], end, 16, tau=a['baseline_tau'][end-1])
    sensitivity = []
    for radius, gain in ((.022, 1.), (.028, 1.), (.034, 1.), (.028, .8), (.028, 1.2)):
        o = {k: v.copy() if isinstance(v, np.ndarray) else v for k,v in original.items()}
        o['my'] *= gain
        fitter = DepthFit(spatial=True, moment=True)
        fitter.model = FinitePatch(Geometry(radius=radius), spatial=True)
        result = fitter.fit(o, float(a['baseline_kraw'][end-1]))
        sensitivity.append(dict(radius=radius, moment_calibration_gain=gain, k=result['k'],
                                data_rank=result['data_rank'], residual=result['residual_rms']))
    report['local_geometry_calibration_sensitivity'] = sensitivity
    # Scientific reproduction of the original negative v0 instance, using the
    # preserved implementation options. Ignore elapsed times and source hashes.
    old = json.loads((HERE/'results/finite_patch_pilot_inverse_v0.json').read_text(encoding='utf-8'))
    old_rows = [r for r in old['rows'] if r.get('family')=='S2' and r.get('model')=='patch_wrench_arrival']
    a = np.load(HERE/'results/finite_patch_input_S2_110000.npz')
    legacy = DepthFit(spatial=True, moment=True, arrival=True, objective='inverse_unweighted',
                      legacy_arrival=True, wrench_nuisance='bias_tangent')
    differences = []
    for r in old_rows:
        o = observed_window(a['measured'], r['end'], 64, tau=a['baseline_tau'][r['end']-1])
        new = legacy.fit(o, r['baseline_kraw'])
        differences.append(abs(new['k']-r['k']))
    assert max(differences) < 1e-6
    report['legacy_v0_reproduction'] = dict(family='S2', model='patch_wrench_arrival', windows=len(old_rows),
                                           max_absolute_k_difference=max(differences))
    target = HERE/'results/finite_patch_audit.json'
    target.write_text(json.dumps(report, indent=2), encoding='utf-8')
    write_summary(report)
    print(json.dumps(report, indent=2))


def write_summary(report):
    names = {'uniform_p1':'同窗口线性模型', 'uniform_arrival':'线性模型＋到达项',
             'patch_normal_arrival':'空间模型＋到达项（法向）',
             'patch_wrench_arrival':'空间模型＋到达项（完整力矩干扰）'}
    lines = ['# 接触模型机制验证汇总', '',
             '数据：v4 只有 S2/T3m 各一条诊断轨迹；窗口不是独立场景。全部为离线因果回放，未改变控制器。', '',
             '增益误差为 |ln(k估计/k评分参考)|。条件力 RMSE 在新编码器输入已知时预测窗口外一拍的滤波力，',
             '不是未来未知运动预测或闭环跟踪误差。没有新确认集或真实机器人证据。', '',
             '|族|模型|窗口数|增益误差中位/P90|条件力RMSE N|局部数值秩不足|未收敛|',
             '|---|---|---:|---:|---:|---:|---:|']
    for fam, methods in report.get('task_data_validation', {}).items():
        for name, r in methods.items():
            lines.append('|%s|%s|%d|%.3f/%.3f|%.4f|%d|%d|' % (fam,names[name],r['n'],
                r['log_error_median'],r['log_error_P90'],r['conditional_force_RMSE'],r['deficient_rank'],r['solver_failures']))
    lines += ['', '局部数值秩按数据雅可比相对奇异值阈值 1e-4 定义，不是全局可辨识证明；',
              '后验包含先验与近似噪声假设，不能直接作为控制安全保证。', '',
              '完整物理/干扰/求解预算检查见 finite_patch_audit.json；解释与后续见 research/contact_model_study.md。',
              '主结果 finite_patch_ablation_v4.json 和五个 finite_patch_input_*.npz 可复现。']
    (HERE/'results/finite_patch_summary.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')


def score_results(path):
    """Held-out FORCE consistency given the new encoder sample, not closed-loop RMSE.

    Fit and surface-motion polynomial see only prefix [:end]. New x/y are the
    conditional input at end; new force is hidden until scoring. This does not
    claim a forecast of unknown future encoder positions or surface motion.
    """
    data = json.loads(path.read_text(encoding='utf-8'))
    rows = data['rows']; caches = {}
    grouped = {}
    decay = math.exp(-2*math.pi*30*.008)
    for r in rows:
        if 'k' not in r:
            continue
        key = r['family'], r['sid']
        if key not in caches:
            caches[key] = np.load(HERE/'results'/f'finite_patch_input_{key[0]}_{key[1]}.npz')
        a = caches[key]; z = a['measured']; end = r['end']
        if end >= len(z):
            continue
        tau = a['baseline_tau'][end-1]
        o = observed_window(z, end, data['protocol']['fit_window'], tau=tau)
        m = FinitePatch(Geometry(exponent=1.5 if r['model']=='uniform_p1.5' else 1.),
                        spatial=r['model'].startswith('patch'))
        theta = np.array(r['params'])
        def force_prediction(model, th):
            dd = model.inverse(o['x'], o['f'], th)
            kk = model.response(o['x'], dd, th)[2]
            beta = np.linalg.lstsq(o['B']*kk[:, None], (o['y']-dd)*kk, rcond=None)[0]
            # Same within-interval timestamp convention as observed_window.
            mm = np.arange(8); af=1-math.exp(-2*math.pi*30*.001)
            ww = af*(1-af)**(7-mm)
            frac = -np.expm1(-mm*.001/tau)/(-math.expm1(-.008/tau))
            q = float(ww@frac/ww.sum())
            xn = z[end-1, 1]+q*(z[end, 1]-z[end-1, 1])
            yn = z[end-1, 2]+q*(z[end, 2]-z[end-1, 2])
            tn = (z[end, 0]-z[end-1, 0])/(data['protocol']['fit_window']*.008)
            sn = np.array([tn**i for i in range(4)])@beta
            fn = float(model.response(xn, max(yn-sn, 0.), th, stiffness=False)[0])
            return decay*z[end-1, 3]+(1-decay)*fn
        pred = force_prediction(m, theta)
        base = force_prediction(FinitePatch(), np.array([math.log(max(r['baseline_kgood'], 1.))]))
        entry = grouped.setdefault((r['family'],r['model']), [])
        entry.append(dict(logerr=abs(math.log(r['k']/r['scoring_k'])) if r['scoring_k'] > 0 else None,
                          force_error=pred-z[end, 3], baseline_force_error=base-z[end, 3],
                          persistence_error=z[end-1, 3]-z[end, 3], rank=r['data_rank'],
                          solver_success=r['solver_success'], sd=r['sd_logk'], transition=r['transition']))
    out = {}
    for (fam, name), items in grouped.items():
        err = [r['logerr'] for r in items if r['logerr'] is not None]
        out.setdefault(fam, {})[name] = dict(n=len(items), log_error_median=float(np.median(err)),
            log_error_P90=float(np.percentile(err,90)),
            conditional_force_RMSE=float(np.sqrt(np.mean([r['force_error']**2 for r in items]))),
            accepted_gain_conditional_RMSE=float(np.sqrt(np.mean([r['baseline_force_error']**2 for r in items]))),
            persistence_force_RMSE=float(np.sqrt(np.mean([r['persistence_error']**2 for r in items]))),
            deficient_rank=sum(r['rank'] < (4 if name.startswith('patch') else 1) for r in items),
            solver_failures=sum(not r['solver_success'] for r in items),
            approximate_1sd_coverage=float(np.mean([r['logerr'] <= r['sd'] for r in items if r['logerr'] is not None])))
    return out


if __name__ == '__main__':
    main()
