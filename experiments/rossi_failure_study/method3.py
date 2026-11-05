"""The Method on the scan-scene plant: force setpoint per scene, probe geometry known to the friction forms, metrics at the scene's force."""
import math
import time

import numpy as np

import contact_candidates as C
import physics_plant3 as P3
import sweeps as S
import common

TOL_F, FMIN, F_OVER = .5, .5, 3.
GATE_ON_PRIOR = True
FF_JUMP = 1.5        # the Method's rule that opens the feed-forward when the gain is plainly stale; without it a thin-cover rib crossing runs away


class ForceSensing:
    """What a real wrist sensor stage does before the Method sees the data (Piedra 2026 section 3): aligns force and encoder, removes the force accelerating the probe.

    Uses only calibrated values (probe mass, servo lag, force time offset), the encoder, the commands and the known force filter.
    """
    def __init__(self, pol, mass, delay_ms):
        self.pol = pol; self.mass = mass; self.delay = max(delay_ms, 0.); self.step = pol.step; self.reset_inner = pol.reset
        if self.delay > 0. and hasattr(pol, '_identify'):
            # Servo identification pairs the encoder with unshifted sent commands.
            identify = pol._identify
            pol._identify = lambda log, n: identify(self.raw_log, n)

    def reset(self, platform):
        self.reset_inner(platform)
        # the calibrated actuator (lag and fast pole), not the Method's online servo estimate, which has a different job
        self.t2 = P3.ACTUATOR_T2; self.t1 = platform['tau_a']-self.t2
        self.y = None; self.vel = 0.; self.x = 0.; self.a = np.zeros(4); self.done = 0; self.comp = np.zeros(4096)

    def __call__(self, o):
        p = self.pol; n = o.n; log = o.log
        if self.delay > 0.:
            self.raw_log = log.copy()
        lag = min(int(self.delay), 2); frac = self.delay-lag
        dt = .0005; nsub = int(round(p.Tc/dt)); t1, t2 = self.t1, self.t2
        for j in range(self.done, n):
            if self.y is None:
                self.y = float(log[j, 2])
            i = j-1-p.ndel
            if j > 0:
                # the same model the axis follows, driven by the sent commands from the first encoder reading
                applied = log[i, 11] if i >= 0 else self.y
                for k in range(nsub):
                    acc = (applied-self.y-(t1+t2)*self.vel)/(t1*t2)
                    if k % 2 == 0:
                        self.a[1:] = self.a[:-1]; self.a[0] = acc
                        self.x += p.afc*(self.mass*((1.-frac)*self.a[lag]+frac*self.a[lag+1])-self.x)
                    self.vel += dt*acc; self.y += dt*self.vel
            self.comp[j] = self.x
        self.done = n
        log[:n, 5] -= self.comp[:n]
        if self.delay > 0. and n > 1:
            # the force is older than the encoder sample, so compare it with the position at the same moment
            shift = self.delay*1.e-3/p.Tc
            y = log[:n, 2].copy()
            log[1:n, 2] = y[1:]-shift*(y[1:]-y[:-1])
        return self.step(o)


TRACKER = dict(r_scale=1., q_maneuver=3.335)     # chosen with equal budgets in select4.py and d_compare.py


def build(name, sc, overrides=None):
    pol = C.build(name, dict(TRACKER, **(overrides or {})))
    if sc.get('tool_mass', 0.) > 0. or sc.get('sensor_delay_ms', 0.) > 0.:
        sensing = ForceSensing(pol, sc.get('tool_mass_hat', 0.), sc.get('sensor_delay_hat', 0.))
        pol.step = sensing; pol.reset = sensing.reset
    pol.force_set = float(sc['FD'])
    pol.gate_on_prior = GATE_ON_PRIOR
    pol.params['ff_jump'] = FF_JUMP
    pol.R_probe = sc['geometry']['radius'] if math.isfinite(sc['geometry']['radius']) else 1.e3
    return pol


def run(name, sc, seed=0, sensor=True):
    pol = build(name, sc)
    t0 = time.perf_counter()
    log = P3.episode(sc, pol, seed, sensor=sensor)
    return log, pol, time.perf_counter()-t0


def metric(log, pol, sc, f_exc=20.8):
    if log is None or len(log) < 10:
        return dict(J=1000., rmse=np.nan, normal_rmse=np.nan, mean_err=np.nan, mean_normal_err=np.nan, peak=np.nan, osc=np.nan, ripple=np.nan, long_loss=True, failed=True,
                    event_rmse=np.nan, call_p50=np.nan, call_p99=np.nan, accept=0., est_re=np.nan, est_abs=np.nan, phi=np.nan, phi_true=np.nan,
                    k_true=np.nan, depth=np.nan)
    fd = sc['FD']; t = log[:, 0]; warm = t >= 1.
    # The Method controls and measures axial force, so J scores it. Normal force differs only through pose and friction, Fn = (Fa - Ft sin th)/cos th,
    # which no contact model can correct at the setpoint, and is reported beside it.
    f = log[warm, 18]; e = f-fd
    normal = log[warm, 3]-fd
    loss = f < FMIN
    longest = cur = 0
    for b in loss:
        cur = cur+1 if b else 0; longest = max(longest, cur)
    J = float(np.mean((e/TOL_F)**2)+20.*np.mean(loss)+np.mean(np.maximum(0., f-fd-F_OVER)**2))
    tt = t[warm]; w = 12
    hp = f-np.convolve(f, np.ones(w)/w, mode='same'); hp = hp[w:-w]; tw = tt[w:-w]
    A = np.c_[np.cos(2*np.pi*f_exc*tw), np.sin(2*np.pi*f_exc*tw)]
    coef = np.linalg.lstsq(A, hp, rcond=None)[0]; rest = hp-A@coef
    ev = np.zeros(len(t), bool)
    for _, a0, a1 in sc.get('events', []):
        ev |= (t >= a0) & (t <= a1)
    em = warm & ev
    F = log[:, 18]
    calls = np.asarray(pol.t_call, float)*1e6
    trace = getattr(pol, 'accept_trace', [])
    good = []; good_abs = []
    for tk, kh, _, _, ok in trace:
        if ok:
            j = int(np.clip(np.searchsorted(t, tk), 0, len(t)-1))
            if log[j, 20] > 0.:
                good.append(math.log(kh/log[j, 20])); good_abs.append(math.log(kh/log[j, 21]))
    return dict(J=min(J, 1000.), rmse=float(np.sqrt(np.mean(e**2))), normal_rmse=float(np.sqrt(np.mean(normal**2))),
                mean_err=float(np.mean(e)), mean_normal_err=float(np.mean(normal)),
                peak=float(f.max()), peak_normal=float(log[warm, 3].max()), over_band=float(np.mean(f > fd+F_OVER)*(t[-1]-1.)), over_15=float(np.mean(f > 15.)*(t[-1]-1.)),
                over_36=bool(f.max() > 36.5), loss_time=float(np.mean(loss)*(t[-1]-1.)), osc=float(np.sqrt(np.mean(rest**2))), ripple=float(np.hypot(*coef)/math.sqrt(2)),
                long_loss=bool(longest*.008 > .5), failed=False, event_rmse=float(np.sqrt(np.mean((F[em]-fd)**2))) if em.sum() else np.nan,
                call_p50=float(np.percentile(calls, 50)), call_p99=float(np.percentile(calls, 99)), accept=float(np.mean([x[4] for x in trace])) if trace else 0.,
                est_re=float(np.median(good)) if good else np.nan, est_abs=float(np.median(good_abs)) if good_abs else np.nan,
                phi=float(getattr(pol, 'phi', 0.)), phi_true=float(np.median(log[warm, 22])), k_true=float(np.median(log[warm, 20])), depth=float(np.median(log[warm, 6])))
