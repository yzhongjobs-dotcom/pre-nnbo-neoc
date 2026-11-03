"""在 V1+Proposed 闭环轨迹上回放三种估计，并用零噪声孪生分解误差。"""
from __future__ import annotations

import argparse
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / 'motion_aware_scan'))
sys.path.insert(0, str(HERE.parent / 'contact_information'))
from cohort_eval import FINAL, cohort  # noqa: E402
import common  # noqa: E402
from joint_policy import JointPolicy  # noqa: E402
from local_contact import fit as lean_fit  # noqa: E402
from conditional_filter import ConditionalLinearContact, MeasuredOnly  # noqa: E402
from v1_servo_id import ServoID  # noqa: E402


def masked(log):
    a = np.full_like(log, np.nan)
    a[:, [0, 1, 2, 5, 11, 12, 13]] = log[:, [0, 1, 2, 5, 11, 12, 13]]
    return a


def replay_new(log, params):
    p = ConditionalLinearContact(single_model=False, **params)
    p.reset(dict(common.P.PLATFORMS['PC125']))
    a = masked(log); rows = []
    for n in range(1, len(log) + 1):
        p(MeasuredOnly(p).inner and common.P.Obs(a, n))
        k = float(p.k_hat) if np.isfinite(p.k_hat) else np.nan
        sig = float(math.sqrt(max(p.cl_filter.P[0, 0], 0.))) if p.cl_filter.initialized else np.nan
        rows.append((k, sig))
    return np.asarray(rows), dict(accepted=p.accepted_count, trials=len(p.accept_history))


def replay_baselines(log):
    a = masked(log); b = ServoID(single_model=False, **dict(FINAL)); b.reset(dict(common.P.PLATFORMS['PC125']))
    rows16, rows32 = [], []
    lean = JointPolicy(wave=None, new_fit=True, rate=False, amplitude_scale=1., **dict(FINAL))
    lean.reset(dict(common.P.PLATFORMS['PC125']))
    for n in range(1, len(log) + 1):
        obs = common.P.Obs(a, n)
        b(obs)
        f = common.window_fit(a, n, FINAL.get('nw', 16), FINAL.get('degree', 3), FINAL.get('noise', .04),
                              b.ns, b.eb, b.afc, b.ndel, FINAL['enc_y'])
        if np.isfinite(f.get('k', np.nan)) and f.get('info', 0.) > 0:
            rows16.append((float(np.clip(f['k'], 200., 3e5)), float(1. / math.sqrt(f['info'])),
                           bool(f['chi2'] <= FINAL['chi2_max'] and 1. / math.sqrt(f['info']) / max(f['k'], 200.) <= FINAL['rel_max'])))
        else:
            rows16.append((np.nan, np.nan, False))
        m = a[:n, [0, 1, 2, 5, 11, 12, 13]]
        r = lean_fit(m, 32, b.tau, b.fc, FINAL.get('noise', .04), rate=False, constant=False)
        if r is None:
            rows32.append((np.nan, np.nan, False))
        else:
            rows32.append((float(np.clip(r['k'], 200., 3e5)), float(r['sigma']),
                           bool(r['risk'] <= math.log(1.2) and r['chi2'] <= FINAL['chi2_max'])))
        lean(obs)
    return np.asarray(rows16, dtype=float), np.asarray(rows32, dtype=float)


def twin(p, baseline_log, T):
    """同一控制指令、零噪声，直接调用 plant.advance。"""
    pl = dict(common.P.PLATFORMS['PC125']); dt = 5e-4; NC = int(round(pl['Tc'] / dt)); nt = int(round(T / pl['Tc']))
    state = common.P.init_state(np.asarray(p, float), dt)
    out = np.zeros_like(baseline_log)
    noise = np.zeros((int(T * 1000) + 10, 2)); extra = np.zeros((int(T * 1000) + 10, 4))
    for n in range(nt):
        cmd = float(baseline_log[n, 11])
        g = np.array([1., 100., 0., 0., 0., cmd])
        got = common.P.advance(state, np.asarray(p, float), g, 1, NC, int(pl['ndel']),
                               float(pl['tau_a']), float(pl['fc']), float(pl['VMAX']), dt,
                               noise, out, n, True, extra)
        if got < 0:
            return None
    return out[:nt]


def one(job):
    scene, params = job
    fam, sid, p, T = scene
    base = ServoID(single_model=False, **dict(FINAL))
    noisy, _ = common.run(MeasuredOnly(base), p, T, sid)
    twin_log = twin(p, noisy, T)
    if twin_log is None:
        return dict(family=fam, sid=sid, error='twin_diverged')
    est_n16, est_n32 = replay_baselines(noisy)
    est_t16, est_t32 = replay_baselines(twin_log)
    pp = dict(FINAL); pp.update(params)
    if 'noise_scale' in pp:
        pp['noise'] = .04 * pp.pop('noise_scale')
    est_nn, qn = replay_new(noisy, pp)
    est_nt, qt = replay_new(twin_log, pp)
    true_n, true_t = noisy[:, 8], twin_log[:, 8]
    out = dict(family=fam, sid=sid, groups=('calm' if fam == 'S2' else 'fast_motion' if fam == 'T2f-stiff' else 'stiffness_change'))
    for name, en, et, truthn, trutht, acc in (
            ('projection16', est_n16, est_t16, true_n, true_t, (est_n16[:, 2] > .5).mean()),
            ('lean32', est_n32, est_t32, true_n, true_t, (est_n32[:, 2] > .5).mean()),
            ('conditional_v0', est_nn, est_nt, true_n, true_t, qn['accepted'] / max(qn['trials'], 1))):
        mn = np.isfinite(en[:, 0]) & np.isfinite(truthn)
        mt = np.isfinite(et[:, 0]) & np.isfinite(trutht)
        mn[:125] = False; mt[:125] = False
        e_struct = et[mt, 0] - trutht[mt]
        both = mn & mt & (np.abs(en[:, 0] - et[:, 0]) < np.inf)
        e_noise = en[both, 0] - et[both, 0]
        coverage = np.mean(np.abs(et[mt, 0] - trutht[mt]) <= 1.96 * np.maximum(et[mt, 1], 1e-12)) if np.any(mt) else np.nan
        out[name] = dict(structure_bias=float(np.mean(e_struct)) if len(e_struct) else None,
                         structure_rmse=float(np.sqrt(np.mean(e_struct ** 2))) if len(e_struct) else None,
                         noise_rmse=float(np.sqrt(np.mean(e_noise ** 2))) if len(e_noise) else None,
                         twin_coverage95=float(coverage), acceptance_rate=float(acc),
                         n_structure=int(len(e_struct)), n_noise=int(len(e_noise)))
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--params', default=None); ap.add_argument('--workers', type=int, default=14); ap.add_argument('--n', type=int, default=50)
    args = ap.parse_args()
    params = json.loads(Path(args.params).read_text(encoding='utf8')) if args.params else {}
    scenes = cohort('diag', args.n)
    rows = []
    with ProcessPoolExecutor(max_workers=min(args.workers, 14)) as pool:
        fs = [pool.submit(one, (s, params)) for s in scenes]
        for i, f in enumerate(as_completed(fs), 1):
            rows.append(f.result())
            if i % 10 == 0: print(i, '/', len(fs), flush=True)
    path = HERE / 'results' / 'open_loop.json'
    path.write_text(json.dumps(dict(protocol=dict(trajectory='V1+Proposed diagnostic closed loop',
                                                   seed_split='110000+', n=args.n,
                                                   confirmation_accessed=False, noise_free_twin='plant.advance',
                                                   params=params), rows=rows), indent=2, default=float), encoding='utf8')
    print(path)


if __name__ == '__main__':
    main()
