"""Paired diagnostics for the saved Huber replay, with frozen estimator settings."""
import json
from time import perf_counter

import footprint_field_filter as F
import numpy as np
from scipy.linalg import cho_factor, cho_solve, solve_triangular
from scipy.optimize import brentq

import corrected_comparison as CC
import memory_segment_eval as SE
import overlay_compare as O

P = SE.P3
CASES = (7, 13, 23, 31, 39)
OUT = F.HERE / 'results' / 'footprint_field_diagnosis.json'
ORIGINAL = F.HERE / 'results' / 'footprint_field_filter.json'
DT = .0005


def rms(values):
    return float(np.sqrt(np.mean(np.square(values))))


def logged_states(sc, seed, compact):
    """Recover omitted diagnostic states without changing the plant or commands."""
    platform = dict(P.PLATFORM, tau_a=sc.get('tau_a_true', P.PLATFORM['tau_a']), log_dt=DT)
    full = P.replay(sc, compact[:, 5], seed=seed, platform=platform, sensor=True)
    tick = int(round(F.TC / DT))
    log = full[::tick]
    saved = full[::tick][:, [0, 1, 3, 26]]
    replay_error = np.max(np.abs(saved - compact[:, [0, 1, 3, 9]]), axis=0)
    np.testing.assert_allclose(saved, compact[:, [0, 1, 3, 9]], rtol=1e-11, atol=1e-12)
    state, _ = P._init(sc, platform)
    depth = full[:, 6]
    with_initial = np.r_[state[22], state[13], depth]
    vd = (3 * depth - 4 * with_initial[1:-1] + with_initial[:-2]) / (2 * DT)
    xc = state[7] + DT * np.cumsum(full[:, 16])
    activation = np.interp(log[:, 0], np.arange(len(sc['act'])) * .001, sc['act'])
    return log, xc[::tick], vd[::tick], np.clip(activation, 0, 1), replay_error


def physical_paths(sc, log, xc, vd, activation, evaluate):
    geometry = P._geometry(sc)
    c = P._constants(sc, P.PLATFORM)
    alpha, tau, _ = P.memory_arrays(sc)
    w = 2 * np.pi * 20.8
    memory = 1 + np.sum(alpha * (1j * w * tau) / (1 + 1j * w * tau))
    zb = complex(c[2] - c[0] * w**2, w * c[1])
    dynamic = c[36] >= 0 and c[33] > 0
    arm_c = 2 * c[36] * np.sqrt(c[3] * c[33]) if dynamic else 0.
    zarm = complex(c[3] - c[33] * w**2, w * arm_c)
    carm = 1 / zarm if dynamic else 1 / c[3]
    lift = complex(c[3], w * arm_c) / zarm if dynamic else 1 + c[33] * w**2 / c[3]
    tissue = (geometry, sc['q_field'], sc['E_field'], sc['ribs'])

    def force(pos, d, theta, v, a):
        return P._contact(pos, d, theta, v, *tissue, a)[0]

    def gain(pos, d, theta, v, a=0., rate=1., hereditary=None):
        fe, k, _, length, _ = P._contact(pos, d, theta, v, *tissue, a)
        if geometry[5]:
            k = (force(pos, d + 1e-6, theta, v, a) - force(pos, d - 1e-6, theta, v, a)) / 2e-6
        if hereditary is None:
            hereditary = fe
        kc = rate * k * memory + 1j * w * c[4] * hereditary
        friction = P._friction(1., v, length, c, geometry) - c[17] * v
        ct, st = np.cos(theta), np.sin(theta)
        axial, lateral = ct + friction * st, friction * ct - st
        series = ct * (axial * (carm + 1 / zb) + c[26] * lateral)
        return abs(kc * axial * ct * lift / (1 + kc * series)), k, length

    depth, v, theta, fn = log[:, 6], log[:, 16], log[:, 17], log[:, 3]
    d0, v0, theta0 = (float(np.median(a[evaluate])) for a in (depth, v, theta))
    keep = (fn > .5) & (log[:, 26] > 0)
    load0 = float(np.median([force(xc[j], depth[j], theta[j], v[j], activation[j])
                             for j in np.flatnonzero(evaluate)]))
    names = ('fixed_depth', 'fixed_load', 'fixed_load_encoder', 'actual_depth',
             'depth_activation', 'depth_activation_pose_speed', 'all', 'elastic_only_no_calibration',
             'fixed_load_contact_length')
    paths = {name: np.full(len(log), np.nan) for name in names}
    for j in np.flatnonzero(keep):
        pos = xc[j]
        rate = max(1 + c[4] * vd[j], 0.)
        if rate <= 0:
            raise ValueError('Positive normal force requires a positive rate multiplier')

        def at_load(x):
            # Frozen pose/speed and activation=0, across the saved spatial profile.
            d = brentq(lambda d: force(x, d, theta0, v0, 0.) - load0, -.02, .10, xtol=1e-13)
            return gain(x, d, theta0, v0)

        paths['fixed_depth'][j] = gain(pos, d0, theta0, v0)[0]
        paths['fixed_load'][j], _, paths['fixed_load_contact_length'][j] = at_load(pos)
        paths['fixed_load_encoder'][j] = at_load(log[j, 1])[0]
        paths['actual_depth'][j] = gain(pos, depth[j], theta0, v0)[0]
        paths['depth_activation'][j] = gain(pos, depth[j], theta0, v0, activation[j])[0]
        (paths['depth_activation_pose_speed'][j], paths['elastic_only_no_calibration'][j],
         _) = gain(pos, depth[j], theta[j], v[j], activation[j])
        paths['all'][j] = gain(pos, depth[j], theta[j], v[j], activation[j], rate, fn[j] / rate)[0]
    np.testing.assert_allclose(paths['all'][keep], log[keep, 26], rtol=2e-6, atol=1e-7)
    meta = dict(depth_reference_m=d0, elastic_force_reference_N=load0,
                theta_reference_rad=theta0, speed_reference_m_per_s=v0,
                hc_lambda=float(c[4]), memory_complex=[memory.real, memory.imag],
                body_support_complex_N_per_m=[zb.real, zb.imag], arm_compliance=[carm.real, carm.imag],
                max_relative_gain_reconstruction=float(np.max(abs(paths['all'][keep] / log[keep, 26] - 1))))
    return paths, meta


def weights(x, length, cell):
    lo = int(np.floor((x.min() - length / 2) / cell))
    hi = int(np.floor((x.max() + length / 2) / cell))
    edges = np.arange(lo, hi + 2) * cell
    left, right = x[:, None] - length / 2, x[:, None] + length / 2
    return np.maximum(np.minimum(edges[1:], right) - np.maximum(edges[:-1], left), 0.) / length


def information(x, length, j, qe):
    """Data-only and Brownian-conditional information near an actual fast rise."""
    first = np.searchsorted(x, x[j] - length)
    cell = length / 24
    A = weights(x[first:j+1], length, cell)
    past, h = A[:-1], A[-1]
    _, s, V = np.linalg.svd(past, full_matrices=False)
    rank = int(np.sum(s > 1e-10 * s[0]))
    observable = V[:rank] @ h
    residual = np.linalg.norm(h - V[:rank].T @ observable) / np.linalg.norm(h)
    _, s_now, V_now = np.linalg.svd(A, full_matrices=False)
    now_rank = int(np.sum(s_now > 1e-10 * s_now[0]))
    naive_sd = np.sqrt(F.R) * np.linalg.norm((V_now[:now_rank] @ h) / s_now[:now_rank])

    index = np.arange(A.shape[1])
    C = 1 + qe * cell * np.minimum.outer(index, index)
    D = past @ C
    cov = D @ past.T + F.R * np.eye(len(past))
    Cminus = C - D.T @ cho_solve(cho_factor(cov), D)
    endpoint_var = h @ Cminus @ h
    Cplus = Cminus - np.outer(Cminus @ h, Cminus @ h) / (endpoint_var + F.R)
    # An adjacent-cell contrast tests a mode largely rejected by the footprint mean.
    center = np.flatnonzero(h)[len(np.flatnonzero(h)) // 2]
    local = np.eye(len(h))[center] - np.eye(len(h))[center - 1]
    return dict(past_ticks=len(past), cells=A.shape[1], past_rank=rank,
                past_nonzero_singular_values=s[:rank].tolist(),
                endpoint_relative_null_projection=float(residual),
                endpoint_data_only_sd_with_current=float(naive_sd),
                endpoint_model_sd_prediction=float(np.sqrt(endpoint_var)),
                endpoint_model_sd_filtered=float(np.sqrt(h @ Cplus @ h)),
                endpoint_noise_equivalent_samples=float(F.R / (h @ Cplus @ h)),
                adjacent_contrast_prior_sd=float(np.sqrt(local @ C @ local)),
                adjacent_contrast_posterior_sd=float(np.sqrt(local @ Cplus @ local)))


def metrics(err, masks):
    return {name: dict(n=int(mask.sum()), rms=rms(err[mask]), bias=float(np.mean(err[mask])),
                       sd=float(np.std(err[mask])))
            for name, mask in masks.items() if mask.any()}


def masks_for(y, evaluate):
    rate = np.gradient(y, F.TC)
    cutoff = np.quantile(rate[evaluate], .90)
    # Offline labels only. Exclude neighborhoods that cross an evaluation gap.
    kernel = np.ones(25) / 25
    smooth_rate = np.gradient(np.convolve(y, kernel, mode='same'), F.TC)
    complete = np.convolve(evaluate.astype(int), np.ones(27), mode='same') == 27
    smooth_cutoff = np.quantile(smooth_rate[complete], .90)
    return dict(all=evaluate, rise=evaluate & (rate > .5), top_positive_decile=evaluate & (rate > max(0., cutoff)),
                top_smoothed_rise=complete & (smooth_rate > max(0., smooth_cutoff)),
                quiet=evaluate & (np.abs(rate) <= .5))


def comparison(ell, good, time, settings, case):
    evaluate = good & (time > 1.5)
    masks = masks_for(ell, evaluate)
    noise = [np.sqrt(F.R) * np.random.default_rng(1000 * case + j).normal(size=len(ell))
             for j in range(F.NSEED)]
    out = {}
    for name, estimate in settings.items():
        err = estimate(ell) - ell
        errors = np.array([estimate(ell + e) - ell for e in noise])
        # Linear in observations: difference isolates the noise contribution.
        noise_only = errors - err
        out[name] = dict(noiseless=metrics(err, masks),
                         noisy={label: dict(n=int(mask.sum()), rms=rms(errors[:, mask]),
                                            bias=float(np.mean(errors[:, mask])),
                                            noise_component_rms=rms(noise_only[:, mask]),
                                            cross_mse=float(2 * np.mean(err[mask] * noise_only[:, mask])))
                                for label, mask in masks.items() if mask.any()})
    return out


def exact_paths(x, length, level, qe):
    H = weights(x, length, length / 24)
    n = H.shape[1]
    rng = np.random.default_rng(20261114)
    cells = np.r_[level, level + np.cumsum(rng.normal(0, np.sqrt(qe * length / 24), n - 1))]
    span = np.arange(n) * (length / 24)
    ribs = level + .60 * ((span % (3 * length)) < .6 * length)
    return {'matched_brownian': H @ cells, 'boxcar_ribs': H @ ribs}


def noise_match(x, length, qe, ell, good, evaluate):
    """Full-prefix batch check and exact output noise, including noisy initialization."""
    A = weights(x[good], length, length / 24)
    i = np.arange(A.shape[1])
    C = 1. + qe * length / 24 * np.minimum.outer(i, i)
    covariance = (A @ C) @ A.T + F.R * np.eye(len(A))
    chol = np.linalg.cholesky(covariance)
    inverse = solve_triangular(chol, np.eye(len(A)), lower=True)
    # For each prefix ending at j, only row j of its inverse covariance is needed.
    coef = np.eye(len(A)) - F.R * inverse / np.diag(chol)[:, None]
    coef[:, 0] += 1. - coef.sum(axis=1)
    err = np.max(np.abs(coef @ ell[good] - F.field(ell, good, x, length, qe)[good]))
    assert err < 1e-9, err
    target = float(F.R * np.mean(np.sum(coef[evaluate[good]]**2, axis=1)))
    gain = 2. * target / (F.R + target)
    q = F.R * gain**2 / (1. - gain)
    return dict(field_noise_sd_exact=float(np.sqrt(target)), time_q_steady_state=float(q),
                max_batch_vs_recursive_error=float(err))


def slope_correlation(x, ell, good, length):
    # Symmetric 24-tick mean, offline only, with all ticks of each neighborhood valid.
    average = np.convolve(np.where(good, ell, 0.), np.ones(24) / 24, mode='same')
    smoothed = np.gradient(average, x)
    valid = (np.convolve(good.astype(int), np.ones(24), mode='same') == 24)
    valid &= np.r_[False, valid[:-1]] & np.r_[valid[1:], False]
    src = np.flatnonzero(valid & (x + length <= x[-1]))
    dest = np.searchsorted(x, x[src] + length)
    frac = (x[src] + length - x[dest - 1]) / (x[dest] - x[dest - 1])
    both = valid[dest] & valid[dest - 1]
    a = smoothed[src[both]]
    b = ((1 - frac) * smoothed[dest - 1] + frac * smoothed[dest])[both]
    return dict(n=len(a), correlation=float(np.corrcoef(a, b)[0, 1])) if len(a) > 2 and a.std() and b.std() else None


def grid_audit(data):
    """Reuse all saved grid values. No fit or additional parameter evaluations."""
    global_key = {name: name + ':' + ':'.join(map(str, values))
                  for name, values in data['selected_global_parameters'].items()}
    rows = data['per_case_grid_evidence']
    groups = {}
    for group in data['summaries']['all']['groups']:
        selected = [c for c in rows if c['group'] == group]
        scopes = {}
        for label, take in [('all', selected), ('held_out', [c for c in selected if c['case'] % 2])]:
            def score(key):
                return float(np.sqrt(np.mean([c['grid'][key]['mse'] for c in take])))

            settings = data['selected_global_parameters']
            q0, cf = settings['field']
            grid = {k: score(k) for k in take[0]['grid'] if k.startswith('field:')}
            delta = {str(c['case']): float(np.sqrt(c['grid'][global_key['field']]['mse'])
                                           - np.sqrt(c['grid'][global_key['rw_t']]['mse'])) for c in take}
            scopes[label] = dict(best_field=min(grid, key=grid.get), best_rms=min(grid.values()),
                                 best_rw_t_rms=min(score(k) for k in take[0]['grid'] if k.startswith('rw_t:')),
                                 field_wins=sum(v < 0 for v in delta.values()), field_minus_rw_t=delta,
                                 rms_by_q_fixed_length={str(q): grid[f'field:{q}:{cf}'] for q in F.QE},
                                 rms_by_length_fixed_q={str(f): grid[f'field:{q0}:{f}'] for f in F.CL})
        groups[group] = scopes
    return groups


def diagnose_case(case, settings):
    i = case['index']
    sc = SE.MC.scene_for(tuple(case['cell']), CC.BASE, case['variant'])
    compact = np.load(F.HERE / f'results/huber_windows/case_{i:02d}.npz')['log']
    log, xc, vd, act, reconstruction = logged_states(sc, case['seed'], compact)
    t, x, fn, k = log[:, 0], log[:, 1], log[:, 3], log[:, 26]
    good = (fn > .5) & (k > 0)
    evaluate = good & (t > 1.5)
    ell = np.log(np.where(k > 0, k, 1.))
    footprint = sc['geometry']['length' if sc['slide'] == 'long' else 'width']
    qe, factor = settings['field']
    length = footprint * factor
    q_time, = settings['rw_t']
    paths, phys = physical_paths(sc, log, xc, vd, act, evaluate)
    pairs = [('fixed_depth', 'all'), ('fixed_load', 'all'), ('actual_depth', 'all'),
             ('depth_activation', 'all'), ('depth_activation_pose_speed', 'all'),
             ('fixed_load', 'fixed_load_encoder'), ('actual_depth', 'depth_activation'),
             ('depth_activation', 'depth_activation_pose_speed'),
             ('depth_activation_pose_speed', 'elastic_only_no_calibration')]
    decomposition = {f'{a}_vs_{b}': metrics(np.log(paths[a] / paths[b]), masks_for(ell, evaluate))
                     for a, b in pairs}
    filters = {'window': lambda obs: F.window(obs, good),
               'time_global': lambda obs: F.rw(obs, good, q_time),
               **{f'field_q{q:g}': lambda obs, q=q: F.field(obs, good, x, length, q)
                  for q in (1., qe, 30.)}}
    if np.all(np.diff(xc) >= 0):
        filters['field_true_center'] = lambda obs: F.field(obs, good, xc, length, qe)
    matched = None
    if i == 39:
        matched = noise_match(x, length, qe, ell, good, evaluate)
        filters['time_noise_matched'] = lambda obs: F.rw(obs, good, matched['time_q_steady_state'])
    real_comparison = comparison(ell, good, t, filters, i)
    out = dict(case=i, group=case['group'], law=sc['contact_law'], nominal_footprint_m=footprint,
               filter_footprint_m=length, reconstruction_max_abs=reconstruction.tolist(),
               speed_m_per_s=float(np.median(np.diff(x)) / F.TC),
               center_offset_rms_mm=1000 * rms((xc - x)[evaluate]), nonforward_center_ticks=int(sum(np.diff(xc) < 0)),
               actual_contact_scan_length_mm=np.quantile(
                   (log[:, 23] if sc['slide'] == 'long' else np.full(len(log), footprint))[evaluate] * 1000,
                   [.05, .50, .95]).tolist(),
               n_evaluated=int(evaluate.sum()), good_timeline=phys,
               physical_decomposition=decomposition, filter_real=real_comparison,
               lag_correlations={name: {label: slope_correlation(x, e, good & (t > 1.5), L)
                                               for label, L in [('effective', length), ('nominal', footprint)]}
                                 for name, e in [('real', ell), ('fixed_load', np.log(paths['fixed_load']))]})
    if i in (31, 39):
        alternative = exact_paths(x, length, float(np.median(ell[evaluate])), qe) if i == 39 else {}
        alternative['fixed_load'] = np.where(good, np.log(paths['fixed_load']), ell)
        out['filter_counterfactuals'] = {name: comparison(e, good, t, filters, i)
                                        for name, e in alternative.items()}
    if i == 39:
        out['noise_match'] = matched
        rate = np.gradient(ell, F.TC)
        j = int(np.argmax(np.where(evaluate, rate, -np.inf)))
        out['information'] = dict(t=float(t[j]), rise_logk_per_s=float(rate[j]),
                                  **information(x, length, j, qe))
        start = j - 23
        # A step in the previously unseen field. The signal is a ramp in the
        # footprint output until the step has moved through the whole footprint.
        step = np.clip((x[start:j+1] - x[start]) / length, 0., 1.)
        delta = .30 * step / step[-1]
        out['information']['new_edge_ramp'] = dict(
            ticks=len(step), endpoint_difference_logk=float(delta[-1]), traversed_fraction=float(step[-1]),
            required_local_step_log=float(.30 / step[-1]), kl=float(delta @ delta / (2 * F.R)))
        # Direct feedthrough in the R -> 0 limit, separately from the noiseless
        # input arm whose estimator still believes the original R.
        original_variance = F.R
        try:
            F.R *= 1e-6
            out['near_exact_observations_rms'] = rms((F.field(ell, good, x, length, qe) - ell)[evaluate])
        finally:
            F.R = original_variance
    return out


def main():
    stamp = perf_counter()
    data = json.loads(ORIGINAL.read_text(encoding='utf-8'))
    O.setup()
    cases = CC.cases()
    output = dict(settings=data['selected_global_parameters'], fixed_noise_variance=F.R,
                  existing_grid=grid_audit(data), diagnostic_cases=[])
    scene_options = [SE.MC.scene_for(tuple(case['cell']), CC.BASE, case['variant']) for case in cases]
    output['scene_option_audit'] = dict(
        n=len(scene_options), hc_enabled=sum(sc['phi_hc'] != 0. for sc in scene_options),
        dynamic_arm_enabled=sum(sc.get('arm_zeta', -1.) >= 0. and sc.get('tool_mass', 0.) > 0.
                                for sc in scene_options),
        active_muscle_cases=[i for i, sc in enumerate(scene_options) if np.any(sc['act'] > 0)])
    for i in CASES:
        row = diagnose_case(cases[i], data['selected_global_parameters'])
        original_grid = data['per_case_grid_evidence'][i]['grid']
        for arm, name in [('window', 'window'), ('time_global', 'rw_t'), ('field_q10', 'field')]:
            key = (name if name == 'window' else
                   name + ':' + ':'.join(map(str, data['selected_global_parameters'][name])))
            np.testing.assert_allclose(row['filter_real'][arm]['noisy']['all']['rms']**2,
                                       original_grid[key]['mse'], rtol=1e-10, atol=1e-12)
        output['diagnostic_cases'].append(row)
        print(f'{i}: {cases[i]["group"]} checked', flush=True)
    output['wall_seconds_including_plant_replay'] = perf_counter() - stamp
    OUT.write_text(json.dumps(output, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print(OUT)


if __name__ == '__main__':
    main()
