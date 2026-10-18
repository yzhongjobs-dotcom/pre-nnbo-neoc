"""Finite-footprint force/moment scan: information vs feedback, not an US simulator.

Run with project numpy/scipy/numba requirements. All parameters are synthetic.
The Winkler foundation is an independent carrier, not B0's Hertz model.
"""
from pathlib import Path
import hashlib
import json
import platform
import time

import numpy as np
from numba import njit
from scipy.optimize import root

HERE = Path(__file__).resolve().parent
A, L, FD, SPEED, SCAN = .015, .010, 3., .008, 10.
DT, DELAY = .01, .02
NAMES = ["fixed_PI", "passive_fit", "preprobe_fit", "embedded_probe_fit", "local_J_oracle"]


@njit
def field(s, t, u, p):
    # p: kappa0, contrast, boundary, width, surface amplitude, wavelength,
    # phase, breathing amplitude, frequency, offset, curvature amplitude
    xx = s + u
    k = p[0] * np.exp(p[1] * np.tanh((xx - p[2]) / p[3]))
    h = (p[9] + p[4] * np.sin(2*np.pi*xx/p[5] + p[6])
         + p[7]*np.sin(2*np.pi*p[8]*t + p[6])
         + p[10]*(xx/.08)**2)
    return k, h


@njit
def wrench(x, s, t, p, nq):
    out = np.zeros(2)
    jac = np.zeros((2, 2))
    for j in range(nq):
        u = -A + (j+.5)*2*A/nq
        r = u/A
        k, h = field(s, t, u, p)
        # x=[normal displacement, A*roll], positive normal = indentation.
        penetration = x[0] + x[1]*r - u*u/(2*.050) - h
        if penetration > 0:
            weight = k*L*2*A/nq
            out[0] += weight*penetration
            out[1] += weight*penetration*r  # moment/A, same units as force
            jac[0, 0] += weight
            jac[0, 1] += weight*r
            jac[1, 0] += weight*r
            jac[1, 1] += weight*r*r
    return out, jac


def environments(seed, n, family):
    rng = np.random.default_rng(seed)
    ps = []
    for _ in range(n):
        ps.append([np.exp(rng.uniform(np.log(1e6), np.log(12e6))),
                   rng.uniform(-.65, .65) if family else 0.,
                   rng.uniform(.025, .055), rng.uniform(.004, .012),
                   rng.uniform(.0004, .0018), rng.uniform(.04, .10),
                   rng.uniform(-np.pi, np.pi),
                   rng.uniform(.0002, .001) if family == 2 else 0.,
                   rng.uniform(.15, .4), rng.uniform(-.002, .002),
                   rng.uniform(-.0004, .0004) if family == 2 else 0.])
    return np.array(ps)


def initial_poses(ps, nq=61):
    xs = []
    for p in ps:
        sol = root(lambda x: wrench(x, 0., 0., p, nq)[0]-[FD, 0], [.020+p[9], 0.],
                   jac=lambda x: wrench(x, 0., 0., p, nq)[1])
        if np.linalg.norm(sol.fun) > 1e-8:
            raise RuntimeError("Initial equilibrium failed")
        xs.append(sol.x)
    return np.array(xs)


@njit
def spd_project(j):
    a, b, c = j[0, 0], .5*(j[0, 1]+j[1, 0]), j[1, 1]
    angle = .5*np.arctan2(2*b, a-c)
    co, si = np.cos(angle), np.sin(angle)
    d = np.sqrt((a-c)**2+4*b*b)
    hi, lo = min(12000., max(80., (a+c+d)/2)), min(12000., max(80., (a+c-d)/2))
    return np.array([[hi*co*co+lo*si*si, (hi-lo)*co*si],
                     [(hi-lo)*co*si, hi*si*si+lo*co*co]])


@njit
def estimate(xh, yh, i, old):
    # Window affine regression with independent linear time drift nuisance.
    # Priors regularize; a full-rank Hessian is NOT data identifiability.
    n = min(i+1, 60)
    if n < 30:
        return old
    xx, yy = xh[i-n+1:i+1].copy(), yh[i-n+1:i+1].copy()
    tm = (n-1)/2
    var_t = n*(n*n-1)/12
    for k in range(2):
        mx, my, sx, sy = 0., 0., 0., 0.
        for q in range(n):
            mx += xx[q, k]/n
            my += yy[q, k]/n
            sx += (q-tm)*xx[q, k]/var_t
            sy += (q-tm)*yy[q, k]/var_t
        for q in range(n):
            xx[q, k] -= mx + sx*(q-tm)
            yy[q, k] -= my + sy*(q-tm)
    gram, rhs = xx.T@xx, yy.T@xx
    # Minimal singular information threshold, explicitly a synthetic choice.
    mineig = .5*(gram[0, 0]+gram[1, 1]-np.sqrt(
        (gram[0, 0]-gram[1, 1])**2+4*gram[0, 1]**2))
    if mineig < 1e-8:
        return old
    ridge = 1e-7
    new = (rhs+ridge*old)@np.linalg.inv(gram+ridge*np.eye(2))
    return .8*old+.2*spd_project(new)


@njit
def episode(p, x0, noise, method, bandwidth, integral, dt=.01, nq=61, substeps=5):
    prep = 1. if method == 2 else 0.
    count = int(round((SCAN+prep)/dt))
    lag = int(round(DELAY/dt))
    x, v, ei = x0.copy(), np.zeros(2), np.zeros(2)
    # Nominal prior shared by every deployable controller.
    jhat = np.array([[1500., 0.], [0., 500.]])
    xh, yh = np.zeros((count+1, 2)), np.zeros((count+1, 2))
    jh = np.zeros((count+1, 2, 2))
    trace = np.zeros((count, 9))
    result = np.zeros(9)
    gain_error = 0.
    for i in range(count):
        t = i*dt
        s = SPEED*max(0., t-prep)
        f, jtrue = wrench(x, s, max(0., t-prep), p, nq)
        xh[i] = x
        yh[i] = f+noise[i]
        jh[i] = jtrue
        obs_idx = max(0, i-lag)
        if method == 4:
            jhat = spd_project(jh[obs_idx])
        elif method > 0 and i % 10 == 0:
            jhat = estimate(xh, yh, obs_idx, jhat)
        err = np.array([FD, 0.])-yh[obs_idx]
        ei = np.minimum(np.maximum(ei+dt*err, -2.), 2.)
        cmd = np.linalg.solve(jhat, bandwidth*err+integral*bandwidth**2*ei)
        probing = (method == 2 and t < prep) or (method == 3 and t < 1.)
        if probing:
            cmd += np.array([.0002*2*np.pi*3*np.cos(2*np.pi*3*t),
                             .00015*2*np.pi*5*np.cos(2*np.pi*5*t)])
        cmd = np.minimum(np.maximum(cmd, np.array([-.020, -.0075])), np.array([.020, .0075]))
        # Exact integral of first-order velocity actuators for held cmd.
        # Only changing contact force/geometry is sampled at substeps.
        h = dt/substeps
        force2, moment2, peak, amin, band = 0., 0., 0., 1e9, 0.
        for sub in range(substeps):
            decay = np.exp(-h/np.array([.020, .030]))
            x += cmd*h+(v-cmd)*np.array([.020, .030])*(1-decay)
            v = cmd+(v-cmd)*decay
            ts = t+(sub+1)*h
            fs, _ = wrench(x, SPEED*max(0., ts-prep), max(0., ts-prep), p, nq)
            force2 += h*(fs[0]-FD)**2
            moment2 += h*(A*fs[1])**2
            peak, amin = max(peak, fs[0]), min(amin, fs[0])
            band += h*(abs(fs[0]-FD) <= .3 and abs(A*fs[1]) <= .01)
        result[0] += force2
        result[1] += moment2
        result[2] = max(result[2], peak)
        result[3] += (amin <= .01)*dt
        if t >= prep:
            result[4] += force2
            result[5] += moment2
            result[6] += band/SCAN
        gain_error += np.linalg.norm(jhat-jtrue)**2*dt
        trace[i] = np.array([t, s, x[0], x[1]/A, fs[0], A*fs[1], jhat[0, 0], jtrue[0, 0], probing])
    result[7] = gain_error/(SCAN+prep)
    result[8] = SCAN+prep
    return result, trace


def noise_for(seed, n, dt=DT):
    # Sensor force/moment noises .03 N and .0003 Nm, shared paired realizations.
    rng = np.random.default_rng(seed)
    return rng.normal(size=(n, int(round(11/dt))+1, 2))*[.03, .0003/A]


def score(r):
    # Normalized wrench squared error, includes preparation; no gain preference.
    return r[..., 0]+r[..., 1]/A**2


def batch(ps, xs, noise, method, bw, integ, **kw):
    return np.array([episode(p, x, nn, method, bw, integ, **kw)[0]
                     for p, x, nn in zip(ps, xs, noise)])


def main():
    start = time.perf_counter()
    # Freeze all tuning before drawing separate held-out outcomes.
    cal = np.concatenate([environments(310+f, 12, f) for f in range(3)])
    xc, nc = initial_poses(cal), noise_for(411, len(cal))
    selected, calibration = {}, {}
    for m, name in enumerate(NAMES):
        rows = []
        for bw in (2., 4., 8., 12., 16., 24.):
            for integ in (0., .1, .3):
                r = batch(cal, xc, nc, m, bw, integ)
                rows.append({"bandwidth": bw, "integral": integ, "mean_total_loss": float(score(r).mean())})
        best = min(rows, key=lambda x: x['mean_total_loss'])
        selected[name] = best
        calibration[name] = rows
        print("selected", name, best, flush=True)
    output, arrays = {}, {}
    for f, label in enumerate(("uniform_static", "heterogeneous_static", "heterogeneous_moving")):
        ps = environments(7000+f, 128, f)
        xs, noise = initial_poses(ps), noise_for(8000+f, len(ps))
        arrays[label+"_parameters"] = ps
        arrays[label+"_initial_poses"] = xs
        rr = []
        for m, name in enumerate(NAMES):
            cfg = selected[name]
            rr.append(batch(ps, xs, noise, m, cfg['bandwidth'], cfg['integral']))
        rr = np.array(rr)
        arrays[label+"_metrics"] = rr
        rows = {}
        for m, name in enumerate(NAMES):
            # Paired difference inference over independent environment episodes.
            difference = score(rr[m])-score(rr[0])
            se = float(difference.std(ddof=1)/np.sqrt(len(ps)))
            rows[name] = dict(total_loss=float(score(rr[m]).mean()),
                scan_force_rms=float(np.sqrt(rr[m, :, 4].mean()/SCAN)),
                scan_moment_rms=float(np.sqrt(rr[m, :, 5].mean()/SCAN)),
                mean_peak_force=float(rr[m, :, 2].mean()),
                worst_peak_force=float(rr[m, :, 2].max()),
                mean_contact_loss_time=float(rr[m, :, 3].mean()),
                illustrative_band_fraction=float(rr[m, :, 6].mean()),
                mean_J_frobenius_error2=float(rr[m, :, 7].mean()),
                total_time=float(rr[m, 0, 8]),
                paired_total_loss_difference=float(difference.mean()),
                approximate_95ci_difference=[float(difference.mean()-1.96*se), float(difference.mean()+1.96*se)])
        output[label] = rows
        # The first held-out trajectory, selected by index before seeing outcomes.
        for m, name in enumerate(NAMES):
            cfg = selected[name]
            arrays[label+"_trace_"+name] = episode(ps[0], xs[0], noise[0], m, cfg['bandwidth'], cfg['integral'])[1]
        print(label, rows, flush=True)
    config = dict(model="independent unilateral Winkler finite footprint; elastic; no ultrasound",
        a=A, length=L, radius=.050, force_target=FD, moment_target=0,
        scan_speed=SPEED, scan_time=SCAN, dt=DT, delay=DELAY,
        actuator_time_constants=[.020, .030], velocity_caps=[.020, .5],
        force_noise_sd=.03, moment_noise_sd=.0003, exact_position_observation=True,
        quadrature_midpoints=61, substeps=5, calibration_episodes=36,
        test_episodes_per_family=128, calibration_seeds=[310,311,312],
        test_seeds=[7000,7001,7002], noise_calibration_seed=411,
        noise_test_seeds=[8000,8001,8002], methods=NAMES,
        tuning_objective="mean total integral ((F-3)^2+(M/a)^2), includes preprobe",
        band_for_illustration_only=[.3,.01], known_environment_parameters_to_controllers=False,
        oracle="exact local wrench Jacobian at delayed observed state; upper reference, no full map",
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (HERE/'results.json').write_text(json.dumps(dict(config=config, selected=selected,
        calibration=calibration, test=output, runtime_s=time.perf_counter()-start,
        python=platform.python_version(), numpy=np.__version__), indent=2)+'\n')
    np.savez_compressed(HERE/'paired.npz', **arrays)


if __name__ == '__main__':
    main()
