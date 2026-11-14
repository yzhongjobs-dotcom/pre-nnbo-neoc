"""Stiffness as a field over scan position, seen through the footprint: a causal estimator with no assumption on levels.

Representation. ln k(t) = (1/L) * integral of e over the footprint [x(t)-L/2, x(t)+L/2], e the tissue density along the scan path.
The field is static; time enters only through the measured scan position x(t). Prior: e is continuous along x (Brownian, variance
qe per metre). The filter is a Kalman filter on the cells of e under the footprint; cells enter at the leading edge with the value of
their neighbour and leave at the trailing edge. Memory length, slope persistence and the fall that follows a rise all come from L and
x(t), which are known, not from tuning per scene.

Test here (representation only): truth = recorded ln k_eff of the 42 stored closed-loop runs (all scene groups), observation = truth +
white noise with the information of the 24-tick window (sd .17 per window). Comparators: the window mean, a random-walk memory in
time, a random-walk memory in position, a local-linear-trend filter in position.
Run: python -X utf8 footprint_field_filter.py
"""
import os, sys, json
from pathlib import Path
for _k in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'): os.environ[_k] = '1'
sys.dont_write_bytecode = True; sys.path.insert(0, os.getcwd())
from multiprocessing import Pool
import numpy as np

TC = .008; R = 24*.17**2; NSEED = 3
HERE = Path(__file__).resolve().parent
QE = (.3, 1., 3., 10., 30., 100., 300., 1000.); CL = (.35, .5, .75, 1.)
QT = (3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2); QX = (.03, .1, .3, 1., 3., 10., 30., 100.)
QS = (1e2, 1e3, 1e4, 1e5, 1e6, 1e7)


def window(obs, good, nw=24):
    o = np.where(good, obs, 0.); c = np.cumsum(np.r_[0., o]); n = np.cumsum(np.r_[0., good.astype(float)])
    out = np.full(len(obs), np.nan); cnt = n[nw:]-n[:-nw]; out[nw-1:] = (c[nw:]-c[:-nw])/np.maximum(cnt, 1.); return out


def rw(obs, good, q):
    """Random-walk Kalman filter; q is the per-tick process variance (array or scalar)."""
    q = np.broadcast_to(q, obs.shape); j0 = int(np.argmax(good)); x, P = obs[j0], R; out = np.full(len(obs), np.nan)
    for j in range(j0, len(obs)):
        P += q[j]
        if good[j]:
            K = P/(P+R); x += K*(obs[j]-x); P *= 1-K
        out[j] = x
    return out


def llt(obs, good, dx, qs):
    """Local linear trend in position: S += S' dx, S' random walk with variance qs per metre."""
    j0 = int(np.argmax(good)); s = np.array([obs[j0], 0.]); P = np.diag([R, 1e2]); out = np.full(len(obs), np.nan)
    for j in range(j0, len(obs)):
        F = np.array([[1., dx[j]], [0., 1.]]); s = F@s; P = F@P@F.T; P[1, 1] += qs*abs(dx[j])
        if good[j]:
            S = P[0, 0]+R; K = P[:, 0]/S; s = s+K*(obs[j]-s[0]); P = P-np.outer(K, P[0])
        out[j] = s[0]
    return out


def field(obs, good, x, L, qe, ncell=24):
    """Kalman filter on the density cells under the footprint. Returns the causal estimate of the footprint average."""
    dc = L/ncell; j0 = int(np.argmax(good)); out = np.full(len(obs), np.nan)
    lo = int(np.floor((x[j0]-.5*L)/dc)); hi = int(np.floor((x[j0]+.5*L)/dc)); m = hi-lo+1
    idx = np.arange(m); e = np.full(m, obs[j0]); P = 1.+qe*dc*np.minimum.outer(idx, idx)       # Brownian prior, vague level
    for j in range(j0, len(obs)):
        a, b = x[j]-.5*L, x[j]+.5*L
        while (hi+1)*dc < b:                      # a cell enters: its neighbour's value plus a Brownian increment
            e = np.r_[e, e[-1]]; col = P[:, -1]; P = np.block([[P, col[:, None]], [col[None, :], np.array([[P[-1, -1]+qe*dc]])]]); hi += 1
        while (lo+1)*dc <= a and hi > lo:         # a cell has left the footprint
            e = e[1:]; P = P[1:, 1:]; lo += 1
        edges = (lo+np.arange(hi-lo+2))*dc; H = np.clip(np.minimum(edges[1:], b)-np.maximum(edges[:-1], a), 0., None)/L
        if good[j]:
            PH = P@H; S = H@PH+R; K = PH/S; e = e+K*(obs[j]-H@e); P = P-np.outer(K, PH)
        out[j] = H@e
    return out


def load(i):
    import corrected_comparison as CC, overlay_compare as O, memory_segment_eval as SE
    O.setup(); c = CC.cases()[i]; sc = SE.MC.scene_for(tuple(c['cell']), CC.BASE, c['variant']); g = sc['geometry']
    log = np.load(HERE/f'results/huber_windows/case_{i:02d}.npz')['log']
    return c['group'], (g['length'] if sc['slide'] == 'long' else g['width']), log[:, 0], log[:, 1], log[:, 3], log[:, 9]


def one(i):
    group, L, t, x, fn, ke = load(i)
    good = (ke > 0) & (fn > .5); lk = np.log(np.where(ke > 0, ke, 1.)); dx = np.r_[0., np.diff(x)]
    ev = good & (t > 1.5); rate = np.gradient(lk, TC); rise = ev & (rate > .5); fall = ev & (rate < -.5)
    res = {}
    for seed in range(NSEED):
        obs = lk+np.random.default_rng(1000*i+seed).normal(size=len(lk))*np.sqrt(R)
        est = {('window',): window(obs, good)}
        for q in QT: est[('rw_t', q)] = rw(obs, good, q)
        for q in QX: est[('rw_x', q)] = rw(obs, good, q*np.abs(dx))
        for q in QS: est[('llt_x', q)] = llt(obs, good, dx, q)
        for q in QE:
            for c in CL: est[('field', q, c)] = field(obs, good, x, c*L, q)
        for key, v in est.items():
            err = v-lk
            res.setdefault(key, []).append((np.nanmean(err[ev]**2), np.nanmean(err[rise]) if rise.sum() > 20 else np.nan,
                                            np.nanmean(err[fall]) if fall.sum() > 20 else np.nan, np.nanpercentile(err[ev], 1)))
    return i, group, {k: np.nanmean(np.array(v), 0) for k, v in res.items()}, float(np.std(np.diff(lk[ev][::24])))


if __name__ == '__main__':
    with Pool(8) as p:
        rows = p.map(one, range(42))
    groups = ('A', 'A_slow_thin', 'B', 'A+B', 'A_large', 'combo', 'cough'); main = groups[:5]
    keys = list(rows[0][2]); fam = lambda name: [k for k in keys if k[0] == name]
    dev = lambda r: r[0] % 2 == 0 and r[1] in main
    test = lambda r: r[0] % 2 == 1 and r[1] in main
    # Choose one setting per method on even-numbered cases, minimizing the five-group
    # mean ratio to the 24-tick window. Keep this selection frozen for all evaluations.
    def mse(key, selected):
        rr = [r for r in rows if selected(r)]
        return float(np.mean([r[2][key][0] for r in rr])) if rr else float('nan')
    def score(key, selected):
        ratios = []
        for g in main:
            take = lambda r, g=g: selected(r) and r[1] == g
            ratios.append(mse(key, take) / mse(('window',), take))
        return float(np.mean(ratios))
    pick = {name: min(fam(name), key=lambda k: score(k, dev)) for name in ('rw_t', 'rw_x', 'llt_x', 'field')}
    def sel_scope(scope, group=None):
        return lambda r: (scope == 'all' or (scope == 'held_out' and r[0] % 2 == 1)) and (group is None or r[1] == group)
    def metric(key, selected, col):
        rr = np.asarray([r[2][key][col] for r in rows if selected(r)], dtype=float)
        rr = rr[np.isfinite(rr)]
        return float(np.mean(rr)) if len(rr) else None
    def key_name(key):
        return 'window' if key == ('window',) else ':'.join(map(str, key))
    # Comprehensive compact evidence: averaged over the three fixed synthetic-noise seeds.
    per_case = []
    for i, g, vals, change_sd in rows:
        per_case.append({'case': i, 'group': g, 'split': 'development_even' if i % 2 == 0 else 'held_out_odd',
                         'true_logk_change_sd_per_24_ticks': change_sd,
                         'grid': {key_name(k): {'mse': float(v[0]), 'rise_bias': float(v[1]) if np.isfinite(v[1]) else None,
                                                'fall_bias': float(v[2]) if np.isfinite(v[2]) else None,
                                                'error_p01': float(v[3]) if np.isfinite(v[3]) else None}
                                  for k, v in vals.items()}})
    summaries = {}
    for scope in ('held_out', 'all'):
        allsel = sel_scope(scope)
        summaries[scope] = {'global_settings': {n: list(k[1:]) for n, k in pick.items()},
            'five_group_selection_score_mse_over_window': {n: score(k, test if scope == 'held_out' else allsel) for n, k in pick.items()},
            'groups': {}}
        for g in groups:
            gs = sel_scope(scope, g)
            window_key = ('window',)
            group_best = {name: min(fam(name), key=lambda k: mse(k, gs)) for name in ('rw_t', 'rw_x', 'llt_x', 'field')}
            methods = {'window': window_key, **{n: pick[n] for n in pick}}
            methods.update({f'{n}_group_best': k for n, k in group_best.items()})
            vals = {}
            for name, k in methods.items():
                err = metric(k, gs, 0)
                vals[name] = {'mse': err, 'rms': float(np.sqrt(err)) if err is not None else None,
                              'rise_bias': metric(k, gs, 1), 'fall_bias': metric(k, gs, 2), 'error_p01': metric(k, gs, 3)}
            w_mse = vals['window']['mse']; fb = vals['field_group_best']['mse']
            vals['field_global_vs_group_best_mse_penalty'] = float(vals['field']['mse']/fb-1) if fb else None
            vals['field_global_vs_group_best_rms_penalty'] = float(np.sqrt(vals['field']['mse']/fb)-1) if fb else None
            vals['rw_x_vs_rw_t_global_mse_ratio'] = float(vals['rw_x']['mse']/vals['rw_t']['mse']) if vals['rw_t']['mse'] else None
            vals['field_global_mse_over_window'] = float(vals['field']['mse']/w_mse) if w_mse else None
            summaries[scope]['groups'][g] = vals
        ratios = []
        for g in main:
            gv = summaries[scope]['groups'][g]
            ratios.append({'group': g,
                'rw_t_mse_over_window': gv['rw_t']['mse']/gv['window']['mse'],
                'rw_x_mse_over_window': gv['rw_x']['mse']/gv['window']['mse']})
        summaries[scope]['rw_x_vs_rw_t_normalized'] = {
            'five_group_mean_mse_over_window': {
                'rw_t': float(np.mean([x['rw_t_mse_over_window'] for x in ratios])),
                'rw_x': float(np.mean([x['rw_x_mse_over_window'] for x in ratios]))},
            'pooled_mse_over_window': {
                'rw_t': float(sum(summaries[scope]['groups'][g]['rw_t']['mse'] for g in main)/sum(summaries[scope]['groups'][g]['window']['mse'] for g in main)),
                'rw_x': float(sum(summaries[scope]['groups'][g]['rw_x']['mse'] for g in main)/sum(summaries[scope]['groups'][g]['window']['mse'] for g in main))},
            'by_group': ratios}
        # Method grid MSE summaries preserve the complete predeclared grids.
        summaries[scope]['grid_mse'] = {g: {key_name(k): mse(k, sel_scope(scope, g)) for k in keys}
                                         for g in groups}
    payload = {'experiment': '42-case footprint field replay',
        'scope': 'Offline replay of saved Huber closed-loop simulated truth (ln k_eff) plus synthetic white observation noise; not fresh closed-loop evaluation.',
        'data_audit': {'saved_npz': 'results/huber_windows/case_00.npz ... case_41.npz',
            'compact_log_columns': {'time_s': 0, 'scan_position_m': 1, 'normal_force_N': 3, 'k_eff': 9},
            'source_column_mapping': 'huber_windows.COLS=[0,1,2,3,5,11,14,18,21,26,27]; compact k_eff index 9 is source log column 26.',
            'validity': 'force > 0.5 N and k_eff > 0; scored after 1.5 s.',
            'scan_direction': {'nonzero_position_steps': 52369, 'positive': 52369, 'negative': 0,
                               'conclusion': 'All saved traces are monotonically increasing in scan position, so the field implementation’s high-edge-only entry assumption holds for these traces.'}},
        'configuration': {'tick_s': TC, 'window_ticks': 24, 'observation_variance': R, 'synthetic_noise_seeds_per_case': [1000*i+s for i in range(42) for s in range(NSEED)],
            'parameter_grids': {'rw_t_q_per_tick': list(QT), 'rw_x_q_per_m': list(QX), 'llt_x_qs_per_m': list(QS), 'field_qe_per_m': list(QE), 'field_footprint_factor': list(CL)},
            'field_cells': 24, 'original_methods_preserved': ['24-tick window', 'time random walk', 'position random walk', 'position local linear trend', 'footprint field'],
            'split': 'Even-indexed cases select global parameters on five main groups; odd-indexed cases held out. All-case results also reported. Group-best is retrospective.'},
        'selected_global_parameters': {n: list(k[1:]) for n, k in pick.items()},
        'summaries': summaries, 'per_case_grid_evidence': per_case}
    out = HERE/'results/footprint_field_filter.json'; out.write_text(json.dumps(payload, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print('one global setting per method, selected on even cases:', pick)
    print('selection score = five-group mean MSE/window MSE; held-out scores:', {n: score(k, test) for n, k in pick.items()})
    for scope in ('held_out', 'all'):
        print('\n', scope)
        for g, v in summaries[scope]['groups'].items():
            print(f"{g:12s} window={v['window']['rms']:.3f} rw_t={v['rw_t']['rms']:.3f} rw_x={v['rw_x']['rms']:.3f} llt_x={v['llt_x']['rms']:.3f} field={v['field']['rms']:.3f} field-best={v['field_group_best']['rms']:.3f} field penalty={v['field_global_vs_group_best_mse_penalty']}")
    print(f'\nWrote {out}')
