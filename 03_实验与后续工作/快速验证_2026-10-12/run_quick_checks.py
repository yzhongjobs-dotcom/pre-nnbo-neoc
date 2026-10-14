"""Bounded R15 feasibility probes; not the full robot or PublishOK certificate.

Sources: current appendix D.1 and M.3. No manuscript or prior result is modified.
Run: python -X utf8 run_quick_checks.py
"""
from pathlib import Path
import csv
import hashlib
import json
import platform
import time
import numpy as np
import scipy
from scipy.linalg import lu_factor, lu_solve
from scipy.optimize import least_squares

OUT = Path(__file__).resolve().parent
SEED = 20261012
U = 0.7

def exact(eta):
    w = (-1 + np.sqrt(1 + 4 * eta**2)) / 2
    if w < U:
        return np.array([w, -w, 0.0])
    w = (eta**2 + U**2) / (1 + 2 * U)
    return np.array([w, -U, 2 * (w - U)])

def F(z, eta):
    w, u, lam = z
    return np.array([-w + eta**2 + u*u + 2*w*u,
                     2*u + 2*w - lam, lam * (u + U)])

def J(z):
    w, u, lam = z
    return np.array([[-1+2*u, 2*u+2*w, 0.],
                     [2., 2., -1.], [0., lam, u+U]])

ANCHOR_ETA = 0.8
ANCHOR = exact(ANCHOR_ETA)
FACTOR = lu_factor(J(ANCHOR))
SENS = lu_solve(FACTOR, np.array([-2*ANCHOR_ETA, 0., 0.]))

def candidate(method, eta):
    if method == 'ANEOC':
        return ANCHOR + SENS * (eta - ANCHOR_ETA)
    if method == 'cached_chord':
        return ANCHOR - lu_solve(FACTOR, F(ANCHOR, eta))
    if method == 'one_Newton':
        return ANCHOR - np.linalg.solve(J(ANCHOR), F(ANCHOR, eta))
    if method == 'structured_ANEOC':
        w = ANCHOR[0] + SENS[0] * (eta-ANCHOR_ETA)
        return np.array([w,-w,0.])
    if method == 'structured_chord':
        w = ANCHOR[0] + (eta*eta-ANCHOR_ETA**2)/(1+2*ANCHOR[0])
        return np.array([w,-w,0.])
    raise ValueError(method)

# In this inactive example F_Z at the old root has no direct eta dependence.
# Thus one Newton at the old root/new parameter and cached chord are identical.

def screen(z, eta, cost_budget):
    """Common analytic toy certificate relative to the finite-root policy.

    Requires inactive branch, 0 <= W <= U and exact inactive stationarity.
    B(W)=eta^2-W-W^2, so |W-W*|=|B|/(1+W+W*) <= |B|/(1+W).
    Recovered feedback u=-(1-t)W z is feasible on z in [0,1].
    z(t)=exp[-W(t-t^2/2)] and |dJ/dW| <=
    2U/3 + 2eta^2/3 + 2U^2/15, by integrating polynomial upper bounds.
    This checks cost error to the finite-root policy, NOT its HJB approximation
    error to the continuous optimum, and is NOT the R15 PublishOK routine.
    """
    w, u, lam = z
    if eta*eta >= U+U*U or not (0 <= w <= U):
        return False, float('inf'), float('inf')
    if abs(u+w) > 1e-12 or abs(lam) > 1e-12:
        return False, float('inf'), float('inf')
    ew = abs(eta*eta-w-w*w)/(1+w)
    lj = 2*U/3 + 2*eta*eta/3 + 2*U*U/15
    return lj*ew <= cost_budget, np.sqrt(2)*ew, lj*ew

NODES, WEIGHTS = np.polynomial.legendre.leggauss(96)
TS, WS = (NODES+1)/2, WEIGHTS/2

def policy_cost(w, eta):
    z = np.exp(-w*(TS-TS*TS/2))
    return float(np.sum(WS * (eta*eta + (1-TS)**2*w*w)*z*z))

def write_csv(name, rows):
    with (OUT/name).open('w', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

def root_probe():
    methods = ['ANEOC', 'cached_chord', 'one_Newton',
               'structured_ANEOC', 'structured_chord']
    for eta in [.72,.8,.82,.88]:
        assert np.linalg.norm(candidate('structured_ANEOC',eta)-candidate('ANEOC',eta))<1e-12
        assert np.linalg.norm(candidate('structured_chord',eta)-candidate('cached_chord',eta))<1e-12
    rows = []
    for h in [.005, .01, .02, .08, .2, .4]:
        eta = ANCHOR_ETA+h
        for method in methods:
            z = candidate(method, eta)
            ok, root_bound, cost_bound = screen(z, eta, 1e-3)
            err = float(np.linalg.norm(z-exact(eta)))
            if np.isfinite(root_bound):
                gap = abs(policy_cost(z[0], eta)-policy_cost(exact(eta)[0], eta))
                assert err <= root_bound+1e-12
                assert gap <= cost_bound+1e-12
                # eta*tanh(eta) is the exact continuous scalar LQR value at t=0,z=1.
                # Its free optimum is feasible in the inactive parameter range used here.
                optimum = eta*np.tanh(eta)
                # The exact unconstrained optimum can violate the action limit
                # before the collocation root changes branch. Do not call it
                # the constrained optimum in that region.
                conf = (policy_cost(exact(eta)[0], eta)-optimum
                        if optimum <= U else None)
                if conf is not None:
                    assert conf >= -1e-12
            else:
                gap, conf = None, None
            rows.append(dict(step=h, eta=eta, method=method, root_error=err,
                             root_bound=root_bound, cost_difference_to_root=gap,
                             cost_bound=cost_bound, accepted_budget_1e_3=bool(ok),
                             finite_root_cost_gap_to_continuous_optimum=conf))
    write_csv('root_steps.csv', rows)
    rng = np.random.default_rng(SEED)
    summary = []
    for radius in [.02, .08, .2]:
        query = ANCHOR_ETA + rng.uniform(-radius, radius, 2000)
        for budget in [1e-4, 1e-3]:
            for method in methods:
                count = sum(screen(candidate(method, eta), eta, budget)[0] for eta in query)
                summary.append(dict(radius=radius, cost_budget=budget, method=method,
                                    query_count=len(query), accepted_count=int(count),
                                    acceptance_rate=count/len(query)))
    write_csv('acceptance.csv', summary)
    # Interleaved repeated batches. Candidate generation plus same analytic check.
    query = ANCHOR_ETA + np.random.default_rng(SEED+1).uniform(-.02, .02, 3000)
    for method in methods:
        for eta in query[:100]:
            screen(candidate(method, eta), eta, 1e-3)
    timings = {m: [] for m in methods}
    for batch in range(7):
        for method in methods[batch%len(methods):]+methods[:batch%len(methods)]:
            t0 = time.perf_counter_ns()
            count = sum(screen(candidate(method, eta), eta, 1e-3)[0] for eta in query)
            elapsed = (time.perf_counter_ns()-t0)/len(query)/1000
            assert count == len(query)
            timings[method].append(elapsed)
    # Incremental sensitivity preparation with existing factorization, actual measurement.
    prep = []
    for _ in range(7):
        t0 = time.perf_counter_ns()
        for j in range(3000):
            lu_solve(FACTOR, np.array([-2*ANCHOR_ETA, 0., 0.]))
        prep.append((time.perf_counter_ns()-t0)/3000/1000)
    timing = {m: dict(median_us=float(np.median(v)), batches_us=v) for m,v in timings.items()}
    timing['ANEOC_incremental_sensitivity_prepare_median_us'] = float(np.median(prep))
    diff = timing['cached_chord']['median_us']-timing['ANEOC']['median_us']
    timing['empirical_prepare_break_even_queries_vs_chord'] = (
        float(np.median(prep))/diff if diff > 0 else None)
    # External correctness check for sensitivity and quadrature, separate from timing.
    fd = (exact(.800001)-exact(.799999))/.000002
    assert np.linalg.norm(fd-SENS) < 1e-8
    from scipy.integrate import quad
    integral = quad(lambda t: (.82**2+(1-t)**2*ANCHOR[0]**2)*
                    np.exp(-2*ANCHOR[0]*(t-t*t/2)), 0, 1, epsabs=1e-12)[0]
    assert abs(integral-policy_cost(ANCHOR[0],.82)) < 1e-11
    return dict(step_results=rows, acceptance=summary, timing=timing,
                chord_newton_identical=True, anchor=ANCHOR.tolist())

def observation(positions, ef, boundary, halfwidth):
    r = np.clip(boundary-positions, -halfwidth, halfwidth)
    root = np.sqrt(np.maximum(halfwidth*halfwidth-r*r,0))
    area = .5*(halfwidth**2*np.arccos(r/halfwidth)-r*root)
    first = root**3/3
    c, e0 = .02/(2*.03), 30000.
    force = c*(e0*np.pi*halfwidth**2/2+(ef-e0)*area)
    moment = c*(ef-e0)*first
    return force, moment

def sensing_probe():
    rng = np.random.default_rng(SEED+2)
    positions = np.linspace(-.002,-.0002,12)
    true_ef, true_b, a = 45000., .0005, .003
    truth = observation(positions,true_ef,true_b,a)
    rows = []
    scenarios = [('nominal_noise',.003,2e-6,1.),
                 ('higher_noise',.01,8e-6,1.),
                 ('width_5percent_mismatch',.003,2e-6,1.05)]
    for name,sf,sm,widthscale in scenarios:
        for trial in range(100):
            yf = truth[0]+rng.normal(0,sf,len(positions))
            ym = truth[1]+rng.normal(0,sm,len(positions))
            for use_moment in [False,True]:
                def residual(q):
                    f,m = observation(positions,q[0]*30000.,q[1]*.003,a*widthscale)
                    rf=(f-yf)/sf
                    return np.r_[rf,(m-ym)/sm] if use_moment else rf
                fit=least_squares(residual,[1.4,.2],bounds=([.8,-1.],[3.,1.]),
                                  xtol=1e-10,ftol=1e-10,gtol=1e-10,max_nfev=150)
                rows.append(dict(scenario=name,trial=trial,method='F+M_history' if use_moment else 'F_history',
                                 boundary_error_mm=(fit.x[1]*.003-true_b)*1000,
                                 modulus_error_percent=(fit.x[0]*30000-true_ef)/true_ef*100,
                                 solver_success=bool(fit.success),nfev=int(fit.nfev)))
    write_csv('sensing_trials.csv',rows)
    summaries=[]
    for name,*_ in scenarios:
        for method in ['F_history','F+M_history']:
            subset=[r for r in rows if r['scenario']==name and r['method']==method]
            b=np.array([r['boundary_error_mm'] for r in subset])
            e=np.array([r['modulus_error_percent'] for r in subset])
            summaries.append(dict(scenario=name,method=method,trials=len(subset),
                                  boundary_rmse_mm=float(np.sqrt(np.mean(b*b))),
                                  boundary_abs_p95_mm=float(np.percentile(abs(b),95)),
                                  modulus_rmse_percent=float(np.sqrt(np.mean(e*e))),
                                  boundary_bias_mm=float(np.mean(b)),
                                  successful_fits=sum(r['solver_success'] for r in subset)))
    write_csv('sensing_summary.csv',summaries)
    return dict(summary=summaries,settings=dict(seed=SEED+2,samples=12,trials=100,
                true_E0=30000,true_Ef=true_ef,true_boundary_m=true_b,halfwidth_m=a,
                positions_m=positions.tolist(),normal_moment_is_idealized=True))

def main():
    started=time.perf_counter()
    print('Running scalar root/feedback checks...',flush=True)
    root=root_probe()
    print('Running paired contact-history estimation probes...',flush=True)
    sensing=sensing_probe()
    result=dict(date='2026-10-12',seed=SEED,python=platform.python_version(),
                numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform(),
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                root=root,sensing=sensing,duration_seconds=time.perf_counter()-started)
    (OUT/'results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(dict(timing=root['timing'],sensing=sensing['summary'],
                         duration_seconds=result['duration_seconds']),ensure_ascii=True,indent=2))

if __name__=='__main__':
    main()
