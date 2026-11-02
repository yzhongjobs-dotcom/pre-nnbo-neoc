"""Contact-only mechanism pilot: causal inverse-depth consistency, not a new controller.

python -X utf8 experiments/contact_model/run_finite_patch_pilot.py --verify-only
python -X utf8 experiments/contact_model/run_finite_patch_pilot.py --scenes 1

At most five diagnosis families by --scenes; no final cohort, parameter search or
plant modification. JSON retains per-window outcomes including failures. npz is
the reproducible measured input plus separately named scoring truth. Model code
cannot access truth, scene parameters, pressure kernels or true indentation.
"""
import argparse
import hashlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from finite_patch_model import FinitePatch, Geometry, orthogonal_residual


def observed_window(measured, end, window, tau=.019, fc=30.):
    """measured columns: t,x,y,filtered Fn,sent command,filtered Ft,filtered My."""
    j = np.arange(max(1, end-window), end)
    decay = math.exp(-2*math.pi*fc*.008)
    F = (measured[j, 3]-decay*measured[j-1, 3])/(1-decay)
    Ft = (measured[j, 5]-decay*measured[j-1, 5])/(1-decay)
    My = (measured[j, 6]-decay*measured[j-1, 6])/(1-decay)
    m = np.arange(8)
    af = 1-math.exp(-2*math.pi*fc*.001)
    w = af*(1-af)**(7-m)
    frac = -np.expm1(-m*.001/tau)/(-math.expm1(-.008/tau))
    q = float(w@frac/w.sum())
    x = measured[j-1, 1]+q*(measured[j, 1]-measured[j-1, 1])
    y = measured[j-1, 2]+q*(measured[j, 2]-measured[j-1, 2])
    # Keep decisions depend only on the available force channel.
    good = (measured[j, 3] > .5) & (measured[j-1, 3] > .5) & (F > .1)
    if good.sum() < max(12, .8*window):
        return None
    tt = (measured[j, 0]-measured[end-1, 0])/(window*.008)
    B = np.column_stack([tt**i for i in range(4)])[good]
    return dict(x=x[good], y=y[good], f=F[good], ft=Ft[good], my=My[good], B=B,
                xend=float(x[-1]), fend=float(F[-1]), tend=float(measured[end-1, 0]))


class DepthFit:
    def __init__(self, exponent=1., spatial=False, moment=False, arrival=False, objective='force_linearized',
                 legacy_arrival=False, wrench_nuisance='full', max_nfev=25, contact='foundation'):
        if contact == 'line':
            from line_patch_model import LinePatch
            self.model = LinePatch(spatial=spatial)
        else:
            self.model = FinitePatch(Geometry(exponent=exponent), spatial)
        self.spatial, self.moment, self.arrival = spatial, moment, arrival
        self.previous = None
        self.objective = objective
        self.legacy_arrival = legacy_arrival
        self.wrench_nuisance, self.max_nfev = wrench_nuisance, max_nfev

    def fit(self, obs, k_init):
        tic = time.perf_counter()
        x, y, F, B = (obs[k] for k in ('x', 'y', 'f', 'B'))
        xe, fe = obs['xend'], max(obs['fend'], .1)
        k_init = float(np.clip(k_init, 200., 3e5))
        p = self.model.g.exponent
        d0 = p*fe/k_init
        lc = math.log(fe/d0**p)
        if hasattr(self.model, 'initial_log_scale'):
            lc = self.model.initial_log_scale(fe, k_init)
        lo = np.array([math.log(20.), -math.log(60), -6., .1]) if self.spatial else np.array([math.log(20.)])
        hi = np.array([math.log(1e11), math.log(60), 6., 6.]) if self.spatial else np.array([math.log(1e11)])
        starts = [np.array([lc, math.log(3), -.7, 1.4]),
                  np.array([lc, -math.log(3), -1., 2.]),
                  np.array([lc, math.log(5), .5, 1.4])] if self.spatial else [np.array([lc])]
        # Offset and lever are nuisance parameters; SVD treats near-collinearity explicitly.
        N = np.column_stack([np.ones(len(F)), obs['ft']] + ([F] if self.wrench_nuisance == 'full' else []))
        prior, precision = None, None
        if self.arrival and self.previous is not None:
            prior, cov, old_x, old_t = self.previous
            prior = prior.copy()
            if self.spatial:
                prior[2] += (old_x-xe)/.01
            # Physical location is transported with encoder motion; stochastic
            # drift scales with ACTUAL elapsed time, including skipped windows.
            dt_scale = 1. if self.legacy_arrival else max((obs['tend']-old_t)/.256, 1.)
            cov = cov+np.diag([.15**2, .3**2, .05**2, .05**2] if self.spatial else [.15**2])*dt_scale
            if self.spatial:
                old_theta = np.array([prior[0], prior[1], xe+.01*prior[2], .01*prior[3]])
                dlast = self.model.inverse(xe, fe, old_theta)
                a = math.sqrt(2*self.model.g.radius*float(dlast))
                if not self.legacy_arrival and xe > old_theta[2]+old_theta[3]+1.2*a:
                    # Local patch has passed out of the footprint. Reopen a new
                    # prospective patch, without reading a scene boundary.
                    prior[1:] = (0., .5, 1.4)
                    cov = np.diag([max(cov[0, 0], .3**2), 3.**2, 4.**2, 4.**2])
            vals, vec = np.linalg.eigh(cov)
            precision = (vec/np.sqrt(np.maximum(vals, 1e-6)))@vec.T
            starts = [prior, *starts]

        def physical(z):
            return np.array([z[0], z[1], xe+.01*z[2], .01*z[3]]) if self.spatial else np.array([z[0]])

        def data_residual(z):
            theta = physical(z)
            d = self.model.inverse(x, F, theta)
            if self.objective == 'inverse_unweighted':
                r = orthogonal_residual(B, y-d)*k_init/.04
            else:
                # Inverting noisy force makes displacement noise sigma_F/K(theta),
                # not a fixed sigma_F/k_init. Normalise in the same force domain.
                # This recovers FORWARD regression exactly for uniform p=1.
                _, _, kk = self.model.response(x, d, theta)
                r = orthogonal_residual(B*kk[:, None], (y-d)*kk)/.04
            if self.moment:
                _, m1, _ = self.model.response(x, d, theta, stiffness=False)
                rm = orthogonal_residual(N, obs['my']-m1)/.001
                r = np.r_[r, rm]
            return r

        def residual(z):
            r = data_residual(z)
            if prior is not None:
                # Overlap correction is approximate; reported as such, not a calibrated posterior.
                # Windows in this sparse pilot are disjoint (stride >= length).
                r = np.r_[r*math.sqrt(.5 if self.legacy_arrival else 1.), precision@(z-prior)]
            elif self.spatial:
                r = np.r_[r, z[1]/3., (z[3]-1.5)/5.]
            return r

        best = None
        for start in starts:
            opt = least_squares(residual, np.clip(start, lo+1e-6, hi-1e-6), bounds=(lo, hi),
                                max_nfev=self.max_nfev, ftol=2e-4, xtol=2e-4, gtol=2e-4,
                                x_scale='jac')
            cost = float(opt.fun@opt.fun)
            if best is None or cost < best[0]:
                best = cost, opt
        _, opt = best; theta = physical(opt.x)
        de = float(self.model.inverse(xe, fe, theta))
        _, _, ke = self.model.response(xe, de, theta)
        rdata = data_residual(opt.x)
        # Numerical data Jacobian is separate from the regularising arrival information.
        J = np.empty((len(rdata), len(opt.x)))
        gradk = np.empty(len(opt.x))
        for i in range(len(opt.x)):
            h = 1e-4; zp, zm = opt.x.copy(), opt.x.copy(); zp[i] += h; zm[i] -= h
            J[:, i] = (data_residual(zp)-data_residual(zm))/(2*h)
            vals = []
            for z in (zp, zm):
                th = physical(z); dd = self.model.inverse(xe, fe, th)
                vals.append(float(self.model.response(xe, dd, th)[2]))
            gradk[i] = (math.log(max(vals[0], 1e-9))-math.log(max(vals[1], 1e-9)))/(2*h)
        sv = np.linalg.svd(J, compute_uv=False)
        rank = int(np.sum(sv > max(sv[0]*1e-4, 1e-8)))
        H = opt.jac.T@opt.jac
        cov = np.linalg.pinv(H, rcond=1e-10)
        sd = float(np.sqrt(max(gradk@cov@gradk, 0.)))
        if self.arrival:
            self.previous = opt.x.copy(), cov, xe, obs['tend']
        return dict(k=float(ke), sd_logk=sd, depth=de, params=theta.tolist(),
                    residual_rms=float(np.sqrt(np.mean(rdata**2))), data_singular_values=sv.tolist(),
                    data_rank=rank, n_parameters=len(opt.x), solver_success=bool(opt.success),
                    evaluations=int(opt.nfev), elapsed=time.perf_counter()-tic)


def verify():
    checks = {}
    for p in (1., 1.5):
        m = FinitePatch(Geometry(exponent=p), spatial=True)
        theta = np.array([math.log(8e4), math.log(3.), .001, .013])
        xs = np.linspace(-.015, .025, 7); ds = np.linspace(.0003, .003, 7)
        f, mom, k = m.response(xs, ds, theta)
        hd = 1e-7
        fd = (m.response(xs, ds+hd, theta)[0]-m.response(xs, ds-hd, theta)[0])/(2*hd)
        derivative_error = float(np.max(np.abs(fd-k)/k))
        roundtrip = float(np.max(np.abs(m.inverse(xs, f, theta)-ds)/ds))
        checks[str(p)] = dict(derivative_relative_error=derivative_error, inverse_relative_error=roundtrip,
                              positive_stiffness=bool(np.all(k > 0)))
        assert derivative_error < .03, checks
        assert roundtrip < 1e-5, checks
        th = np.array([math.log(8e4), 0., .001, .013])
        f0, m0, k0 = m.response(xs, ds, th)
        np.testing.assert_allclose(f0, 8e4*ds**p, rtol=1e-10)
        np.testing.assert_allclose(k0, p*8e4*ds**(p-1), rtol=1e-10)
        assert np.max(np.abs(m0)) < 1e-12
        assert m.response(0., 0., th)[0] == 0
    B = np.c_[np.ones(20), np.ones(20)*2.]
    np.testing.assert_allclose(orthogonal_residual(B, np.arange(20.)), np.arange(20.)-9.5, atol=1e-12)
    # AuSoScan printed static limit: E* and L are arbitrary positive test values.
    K = math.pi*.04*6300/4
    checks['ausoscan_static_limit'] = dict(force_at_1mm=K*.001, derivative_from_13_15=K,
                                         derivative_printed_17_at_v_zero=0.)
    return checks


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--verify-only', action='store_true')
    ap.add_argument('--scenes', type=int, default=1); ap.add_argument('--max-windows', type=int, default=10)
    ap.add_argument('--window', type=int, choices=(16, 32, 64), default=16)
    ap.add_argument('--families', nargs='+', default=None)
    ap.add_argument('--models', nargs='+', default=None)
    ap.add_argument('--stride', type=int, default=32)
    ap.add_argument('--legacy-arrival', action='store_true', help='reproduce the v0/v1 sparse-arrival implementation')
    ap.add_argument('--label', default='finite_patch_pilot')
    ap.add_argument('--wrench-nuisance', choices=('full', 'bias_tangent'), default='full')
    ap.add_argument('--max-nfev', type=int, default=25)
    ap.add_argument('--reuse', action='store_true')
    ap.add_argument('--objective', choices=('inverse_unweighted', 'force_linearized'), default='force_linearized')
    args = ap.parse_args()
    if args.stride < args.window and not args.legacy_arrival:
        ap.error('this pilot requires nonoverlapping windows, stride >= window')
    if not args.label.replace('_', '').isalnum():
        ap.error('label must contain only letters, digits and underscores')
    checks = verify()
    if args.verify_only:
        print(json.dumps(checks, indent=2)); return
    sys.path.insert(0, str(HERE.parent/'motion_aware_scan'))
    sys.path.insert(0, str(HERE.parent/'method_redesign'))
    import common, cohort_eval as CE
    from v1_servo_id import ServoID
    outdir = HERE/'results'; outdir.mkdir(exist_ok=True)
    class Record(ServoID):
        def reset(self, platform):
            super().reset(platform); self.k_history = []; self.raw_history = []; self.tau_history = []
        def step(self, o):
            cmd = super().step(o)
            self.k_history.append(self.k_good); self.raw_history.append(self.k_hat)
            self.tau_history.append(self.tau)
            return cmd
    records = []; scene_metrics = []
    for fam, sid, scene, T in CE.cohort('diag', args.scenes):
        if args.families and fam not in args.families:
            continue
        cache = outdir/f'finite_patch_input_{fam}_{sid}.npz'
        if args.reuse and cache.exists():
            a = np.load(cache); measured, truth = a['measured'], a['scoring_truth_k']
            kgood, kraw, taus = a['baseline_kgood'], a['baseline_kraw'], a['baseline_tau']
        else:
            policy = Record(single_model=False, **CE.FINAL)
            log, _ = common.run(policy, scene, T, sid)
            measured = log[:, [0, 1, 2, 5, 11, 12, 13]].copy()
            truth = log[:, 8].copy(); kgood = np.array(policy.k_history); kraw = np.array(policy.raw_history)
            taus = np.array(policy.tau_history)
            np.savez_compressed(cache, measured=measured, scoring_truth_k=truth,
                                baseline_kgood=kgood, baseline_kraw=kraw, baseline_tau=taus)
        shared = dict(objective=args.objective, legacy_arrival=args.legacy_arrival,
                      wrench_nuisance=args.wrench_nuisance, max_nfev=args.max_nfev)
        fits = {'uniform_p1': DepthFit(**shared), 'uniform_p1.5': DepthFit(1.5, **shared),
                'uniform_arrival': DepthFit(arrival=True, **shared),
                'patch_normal': DepthFit(spatial=True, **shared),
                'patch_normal_arrival': DepthFit(spatial=True, arrival=True, **shared),
                'patch_wrench': DepthFit(spatial=True, moment=True, **shared),
                'patch_wrench_arrival': DepthFit(spatial=True, moment=True, arrival=True, **shared),
                'line_uniform_arrival': DepthFit(contact='line', arrival=True, **shared),
                'line_patch_normal_arrival': DepthFit(contact='line', spatial=True, arrival=True, **shared),
                'line_patch_wrench_arrival': DepthFit(contact='line', spatial=True, moment=True, arrival=True, **shared)}
        if args.models:
            fits = {name: fits[name] for name in args.models}
        ends = np.arange(128, len(measured)+1, args.stride)
        ends = ends[measured[ends-1, 0] >= 1.]
        # Evenly spaced observation times, fixed before inspecting any scoring truth.
        if args.max_windows > 0 and len(ends) > args.max_windows:
            ends = ends[np.linspace(0, len(ends)-1, args.max_windows).astype(int)]
        for end in ends:
            obs = observed_window(measured, end, args.window, tau=taus[end-1])
            if obs is None:
                records.append(dict(family=fam, sid=sid, end=int(end), kind='unavailable_contact_window')); continue
            initial = kraw[end-1]
            if not np.isfinite(initial) or initial <= 0:
                initial = kgood[end-1]
            if not np.isfinite(initial) or initial <= 0:
                initial = 3e4
            base = dict(family=fam, sid=sid, end=int(end), time=obs['tend'], scoring_k=float(truth[end-1]),
                        baseline_kgood=float(kgood[end-1]), baseline_kraw=float(kraw[end-1]))
            # Truth is introduced only here for labels and scores, never inside fit/observed_window.
            kt = truth[max(end-17, 0):end]
            base['transition'] = bool(np.all(kt > 0) and np.max(np.abs(np.log(kt/kt[-1]))) > .1)
            for name, fitter in fits.items():
                try:
                    r = fitter.fit(obs, initial); records.append(dict(base, model=name, **r))
                except Exception as exc:
                    records.append(dict(base, model=name, error=repr(exc)))
            print(f'{fam} sid={sid} t={obs["tend"]:.3f} finished', flush=True)
        scene_metrics.append(dict(family=fam, sid=sid, input=str(cache.name), samples=len(measured)))
    output = dict(protocol=dict(cohort='diag', scenes_per_family=args.scenes, max_windows=args.max_windows,
                                fit_window=args.window, moment_sigma=.001, force_sigma=.04, objective=args.objective,
                                legacy_arrival=args.legacy_arrival, stride=args.stride, families=args.families,
                                wrench_nuisance=args.wrench_nuisance, max_nfev=args.max_nfev,
                                control_changed=False, final_used=False), verification=checks,
                  scenes=scene_metrics, rows=records,
                  source_hashes={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in
                                 (Path(__file__), HERE/'finite_patch_model.py', HERE/'line_patch_model.py')})
    target = outdir/(args.label+'.json')
    target.write_text(json.dumps(output, indent=2, allow_nan=True), encoding='utf-8')
    for fam in CE.FAMS:
        print(fam)
        rr = [r for r in records if r.get('family') == fam and 'scoring_k' in r and r['scoring_k'] > 0]
        for name in ['baseline_kgood', 'baseline_kraw', *fits]:
            reference_model = 'uniform_p1' if 'uniform_p1' in fits else next(iter(fits))
            rows = [r for r in rr if r.get('model') == (reference_model if name.startswith('baseline') else name)]
            vals = [(r.get(name) if name.startswith('baseline') else r.get('k')) for r in rows]
            err = [abs(math.log(v/r['scoring_k'])) for r, v in zip(rows, vals) if v is not None and np.isfinite(v) and v > 0]
            if err:
                print('  %s %d/%d median/P90 %.3f/%.3f' % (name, len(err), len(rows), np.median(err), np.percentile(err, 90)))
    print(target)


if __name__ == '__main__':
    main()
