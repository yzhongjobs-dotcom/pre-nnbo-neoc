"""V4: C1 GN decision with FINAL motion measurements/IMM and V1 servo ID.
The copied numerical kernels retain the original objective, probes and solver.
Only the nominal rollout drift is replaced by a tabulated force increment.
No true scene/contact fields are read. Run v4_evaluate.py for reproduction.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import sys, math, time
from pathlib import Path
import numpy as np
from numba import njit
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'motion_aware_scan'))
sys.path.insert(0, str(HERE.parent/'cobot_scan_benchmark'))
import common, variants, cohort_eval as CE
import plant as P
from contact_update import ContactUpdate, profile_contact
from candidate1_model import _sync, XLO, XHI, LDMIN, LDMAX, NX, FDH, FD
from v1_servo_id import ServoID

@njit(cache=True)
def _roll(x, H, per, st, q, dith, pr, kmult, dh, wsc, out, o0, mode, vd):
    """closed-loop predictor of the exact robot side (admittance ZOH, KPF path, u_ff, ndel delay,
    first-order servo lag at 1 kHz, 1 kHz first-order F/T filter) through the estimated contact
    model with the surface motion extrapolated. Starts at the last logged tick t_{n-1}; interval 0
    is already committed. Writes 3 residuals per predicted tick t_{n+1}..t_{n+H} and one terminal term.
    mode 0: from the current (replicated) state with the estimated drift; dh = surface step [m].
    mode 1: robustness probe from equilibrium (F = FD, at rest, no feedforward/dither) against an
            unpredicted surface ramp of velocity vd on a contact of stiffness kmult * k_hat."""
    M = pr[0]; K = pr[1]; ki = pr[2]; Tc = pr[3]; VMAX = pr[4]; eb = pr[5]; afc = pr[6]; ns = int(pr[7])
    ga = pr[8]; tv = pr[9]; Fc = pr[10]; wc = math.sqrt(pr[11]); Fl = pr[12]; wl = math.sqrt(pr[13])
    D = math.exp(x[0]); KP = x[1] * 1e-4; u0 = x[2] * 1e-3; u1 = x[3] * 1e-2
    yc = st[0]; v = st[1]; I = st[2]; yref = st[3]; y = st[4]; Ff = st[5]; tgt = st[6]
    g0 = st[7]; k = st[8] * kmult; vs = st[9]; acc = st[10]; ye = st[11]
    nd = q.shape[0]
    buf = np.empty(max(nd, 1))
    for i in range(nd): buf[i] = q[i]
    if mode == 1:
        yc = ye; v = 0.; I = 0.; yref = ye; y = ye; Ff = FD; tgt = ye; g0 = FD; acc = 0.; vs = vd; u0 = 0.; u1 = 0.
        tv = 1e9
        for i in range(nd): buf[i] = ye
    e = math.exp(-D * Tc / M)
    sH = wsc / math.sqrt(H)
    m = o0; f = 0.
    for i in range(H + 1):
        for mm in range(1, ns + 1):
            y = tgt + (y - tgt) * eb
            tau = (i * ns + mm) * 1e-3
            tt = min(tau, tv)
            hh = vs * tt + 0.5 * ga * acc * tt * tt + dh
            if mode == 0 and st.shape[0] > 16:
                # Table is in N, not metres: the measurement gain, not the
                # feedback/model stiffness, converts equivalent motion to force.
                f = g0 + k * (y - ye) + kmult * st[16+i*ns+mm-1]
            else:
                f = g0 + k * (y - ye + hh)
            flin = f
            if f < 0.: f = 0.
            Ff += afc * (f - Ff)
        if i >= 1:
            out[m] = sH * (f - FD) / 0.5
            out[m + 1] = sH * wc * (f - Fc if f > Fc else 0.)
            fl_ = flin if pr[16] > 0 else f      # lossfix: barrier on the unclipped model force (gradient in loss)
            out[m + 2] = sH * wl * ((Fl - fl_) / 0.5 if fl_ < Fl else 0.)
            m += 3
        if i == H: break
        ef = FD - Ff
        I = min(max(I + ef * Tc, -20.), 20.)
        uu = ef + ki * I - K * (yc - yref)
        vn = v * e + (1 - e) * uu / D
        vn = min(max(vn, -VMAX), VMAX)
        yc += 0.5 * (v + vn) * Tc; v = vn
        blk = (i // per) * per
        uf = u0 + u1 * min(blk * Tc, pr[15]) + (dith[i] if mode == 0 else 0.)
        if pr[16] > 0: uf = min(max(uf, -0.02), 0.02)     # lossfix: plant clip of u_ff known to the predictor
        ycmd = yc + KP * ef + uf
        if nd > 0:
            tgt = buf[0]
            for kk in range(nd - 1): buf[kk] = buf[kk + 1]
            buf[nd - 1] = ycmd
        else:
            tgt = ycmd
    # terminal velocity-consistency term: the admittance (not the feedforward offset) should carry the
    # extrapolated surface velocity at the horizon end (prevents integrator / feedforward cancellation)
    tend = (H + 1) * ns * 1e-3
    vse = (vs + ga * acc * min(tend, tv)) if tend < tv else 0.0
    out[m] = wsc * math.sqrt(pr[14]) * (v + vse) / 0.01
    return m + 1


@njit(cache=True)
def _resid(x, H, per, st, q, dith, pr, scen, prev, prior, rw, out):
    """stacked residual vector: per scenario (kmult, surface step dh, weight) 3H force/barrier terms,
    then 7 regularisation terms (gain/feedforward rate to the previously applied decision, gain prior,
    feedforward prior)."""
    m = 0
    for s in range(scen.shape[0]):
        m = _roll(x, H, per, st, q, dith, pr, scen[s, 0], scen[s, 1], math.sqrt(scen[s, 2]), out, m,
                  int(scen[s, 3]), scen[s, 4])
    out[m] = rw[0] * (x[0] - prev[0]); out[m + 1] = rw[1] * (x[1] - prev[1]); out[m + 2] = rw[2] * (x[2] - prev[2])
    out[m + 3] = rw[3] * (x[0] - prior[0]); out[m + 4] = rw[4] * (x[1] - prior[1])
    out[m + 5] = rw[5] * x[2]; out[m + 6] = rw[6] * x[3]
    return m + 7


@njit(cache=True)
def _jac(x, r, H, per, st, q, dith, pr, scen, prev, prior, rw, lo, hi, J, r2):
    nr = r.shape[0]
    for j in range(NX):
        h = FDH[j]
        if x[j] + h > hi[j]: h = -h
        xj = x[j]; x[j] = xj + h
        _resid(x, H, per, st, q, dith, pr, scen, prev, prior, rw, r2)
        x[j] = xj
        for i in range(nr): J[i, j] = (r2[i] - r[i]) / h


@njit(cache=True)
def _boxqp(Hm, g, lo, hi):
    """min 0.5 d'Hd + g'd on lo <= d <= hi (convex, small n): projected-Newton active set, then
    projected Gauss-Seidel polish."""
    n = g.shape[0]
    d = np.zeros(n)
    for i in range(n): d[i] = min(max(0., lo[i]), hi[i])
    free = np.ones(n, dtype=np.bool_)
    for it in range(8):
        gr = Hm @ d + g
        nf = 0
        for i in range(n):
            free[i] = not ((d[i] <= lo[i] + 1e-12 and gr[i] > 0) or (d[i] >= hi[i] - 1e-12 and gr[i] < 0))
            if free[i]: nf += 1
        if nf == 0: break
        idx = np.empty(nf, dtype=np.int64); t = 0
        for i in range(n):
            if free[i]: idx[t] = i; t += 1
        Hf = np.empty((nf, nf)); rf = np.empty(nf)
        for a in range(nf):
            rf[a] = -gr[idx[a]]
            for c in range(nf): Hf[a, c] = Hm[idx[a], idx[c]]
        step = np.linalg.solve(Hf, rf)
        moved = 0.
        for a in range(nf):
            i = idx[a]; dn = min(max(d[i] + step[a], lo[i]), hi[i]); moved += abs(dn - d[i]); d[i] = dn
        if moved < 1e-12: break
    for sweep in range(20):
        for i in range(n):
            gi = g[i]
            for c in range(n): gi += Hm[i, c] * d[c]
            if Hm[i, i] > 0:
                d[i] = min(max(d[i] - gi / Hm[i, i], lo[i]), hi[i])
    return d


@njit(cache=True)
def _cost(r):
    return 0.5 * np.sum(r * r)


@njit(cache=True)
def _solve_ord(x, iters, H, per, st, q, dith, pr, scen, prev, prior, rw, lo, hi, Jc, rc, info):
    """ordinary re-solve: warm-started projected Gauss-Newton / Levenberg-Marquardt with cost
    acceptance. Returns final cost; info[0] = rollouts used. Leaves J and r at the final point in Jc, rc
    (the last Jacobian evaluated, i.e. at the start point of the last iteration)."""
    nr = rc.shape[0]; nsc = scen.shape[0]
    r = np.empty(nr); r2 = np.empty(nr); rn = np.empty(nr); J = np.empty((nr, NX))
    _resid(x, H, per, st, q, dith, pr, scen, prev, prior, rw, r); nroll = nsc
    c = _cost(r); lamb = 1e-3
    for it in range(iters):
        _jac(x, r, H, per, st, q, dith, pr, scen, prev, prior, rw, lo, hi, J, r2); nroll += NX * nsc
        Hm = J.T @ J; g = J.T @ r
        for i in range(NX): Hm[i, i] += lamb * (Hm[i, i] + 1e-6)
        d = _boxqp(Hm, g, lo - x, hi - x)
        xn = x + d
        _resid(xn, H, per, st, q, dith, pr, scen, prev, prior, rw, rn); nroll += nsc
        cn = _cost(rn)
        Jc[:, :] = J
        if cn < c:
            x[:] = xn; r[:] = rn; c = cn; lamb = max(lamb * 0.3, 1e-6)
        else:
            lamb *= 10.
    rc[:] = r
    info[0] = nroll
    return c


class V4(ContactUpdate):
    W = ServoID.W

    def __init__(self, motion=True, motion_weight=1., ru=0., **kw):
        super().__init__(ru=ru, **kw)
        self.motion = bool(motion)
        self.motion_weight = float(motion_weight)

    def reset(self, platform):
        super().reset(platform)
        self.tau_hat, self.ndel_hat = self.tau, self.ndel
        self.observer = variants.build('Proposed-IMM', dict(CE.FINAL))
        self.observer.reset(platform)
        self.motion_ready = False
        self.counts = dict(ticks=0, motion_ticks=0, recovery_ticks=0,
                           gain_ticks=0, uff_bound_ticks=0)
        self.force_drift_max = 0.

    def _apply(self, tau, ndel):
        self.tau, self.ndel = tau, int(ndel)
        self.eb = math.exp(-1e-3/tau) if tau > 0 else 0.
        self.pr[5] = self.eb
        ServoID._apply(self.observer, tau, ndel)

    def _observe(self, o):
        ob = self.observer
        F = o.log[o.n-1, 5]
        self.motion_ready = False
        if self.mode or F < .3:
            ob._clear_motion()
            return
        ob.kbar, _ = ob.stiffness(o, o.n)
        ob.drift_pred(o, o.n)
        self.motion_ready = (np.isfinite(ob.k_good) and ob.k_good > 0
                             and ob.tracker.initialized)

    def _motion_state(self, st):
        # The motion measurement is an interval-weighted quantity. As in
        # FINAL's (ndel+0.5)*Tc prediction, half a tick aligns its timestamp.
        # No additional tau*velocity: the rollout already simulates the servo.
        ob = self.observer
        drift = np.zeros((self.H+1)*self.ns)
        st[9:11] = 0.
        if self.motion_ready:
            h = (np.arange(len(drift))+1)*.001
            v, a = ob.tracker.x[1:3]
            h0 = .5*self.Tc
            drift = self.motion_weight*ob.k_good*(v*h+a*(h0*h+.5*h*h))
            self.force_drift_max = max(self.force_drift_max, float(np.max(abs(drift))))
            self.th[2] = self.motion_weight*ob.k_good/self.kuse*(v+a*h0)
            self.th[3] = self.motion_weight*ob.k_good/self.kuse*a
        return np.concatenate((st, drift))

    def _record(self, o, started, cost, rollouts):
        super()._record(o, started, cost, rollouts)
        self.counts['ticks'] += 1
        self.counts['motion_ticks'] += int(self.motion_ready and not self.mode)
        self.counts['gain_ticks'] += int(np.isfinite(self.observer.k_good))
        self.counts['recovery_ticks'] += int(bool(self.mode))
        self.counts['uff_bound_ticks'] += int(abs(self.g[5]) >= .0199)

    def __call__(self, o):
        started = time.perf_counter()
        n, log = o.n, o.log
        if n > 60 and n % 25 == 0:
            ServoID._identify(self, log, n)
        if self.motion and n:
            self._observe(o)
        if n == 0:
            D, KP = self.rule(self.kuse)
            self.x = np.array([math.log(D), KP/1e-4, 0., 0.])
            self.g = (self.M, D, 0., self.ki, KP, 0.)
            return self.g
        if self.rs is None:
            self.rs = np.array([log[0, 2], 0., 0., log[0, 2]])
        _sync(log, self.m, n, self.rs, *self.g[:4], self.Tc, self.VMAX)
        self.m = n
        F, y = log[n-1, 5], log[n-1, 2]
        if self.mode == 0 and F < .3:
            self.mode = 1
            self.motion_before_loss = np.clip(self.th[2], -.04, .04)
            self.Ioff = self.rs[2]
            self.inband = 0
        if self.mode:
            # Keep the contact-stiffness rule at re-contact. Fixed gains tuned
            # on a different tissue are not used. During loss close relative
            # to the last estimated drift, rather than removing all damping.
            D = np.clip(self.kuse/self.recovery_damping, 100., P.GMAX[1])
            KP = min(self.recovery_cp/self.kuse, P.GMAX[4])
            if F < .5:
                v = np.clip(-self.motion_before_loss+self.approach, -.04, .06)
                target = log[n-1, 11]+v*self.Tc
                self.inband = 0
            else:
                # A force-level correction advances the commanded reference,
                # so the first touch cannot introduce a proportional jump.
                v = np.clip((P.FD-F)/(self.kuse*.04), -.03, .03)
                target = log[n-1, 11]+v*self.Tc
                self.inband = self.inband+1 if 1. <= F <= 6. else 0
            target = np.clip(target, y-.003, y+.003)
            yn = self._next_admittance(F, D, KP, 0., 0.)
            uff = float(np.clip(target-yn-KP*(P.FD-F), -.02, .02))
            self.g = (self.M, D, 0., 0., KP, uff)
            self.th[:] = [math.log(self.kuse), F, 0., 0.]
            self.x = np.array([math.log(D), KP/1e-4, uff/1e-3, 0.])
            if self.inband >= 12:
                self.mode = 0
                self.rc = n-1
                self.Ioff = self.rs[2]
            self._record(o, started, 0., 0.)
            return self.g

        nw = self.Nw if self.rc == 0 else min(self.Nw, n-1-self.rc)
        legacy = None
        if not self.projection or self.projection == 'adaptive':
            from candidate1_model import _est
            legacy = self.th.copy()
            mu = legacy.copy()
            mu[2] += legacy[3]*self.Tc
            pw = np.array([1/self.sk**2,0.,1/self.sv**2,1/self.sa**2])
            info = np.zeros(2)
            valid = _est(log,n,nw,self.lam,self.ndel,self.eb,self.afc,self.ns,
                         legacy,mu,pw,self.nfree,self.est_lo,self.est_hi,
                         self.iters_est,self.sz,info)
            if valid == 0:
                legacy[:] = [mu[0],F,0.,0.]
            if abs(legacy[2]) > self.motion_threshold:
                self.project_until = o.t[-1]+self.gate_release
        project = bool(self.projection) and (self.projection != 'adaptive' or o.t[-1]<self.project_until)
        if project:
            fit = profile_contact(log, n, nw, self.lam, self.ndel,
                                  self.eb, self.afc, self.ns, self.kuse,
                                  self.relative_prior, self.noise, self.degree)
            k = max(fit[0], self.kuse*math.exp(-self.Tc/self.trel)) if self.trel else fit[0]
            self.kuse = k
            self.th[:] = [math.log(k), np.clip(fit[1], -5., 30.),
                          np.clip(fit[2], -.1, .1), np.clip(fit[3], -5., 5.)]
            self.information = fit[4]
        else:
            self.th[:] = legacy
            self.kuse = max(math.exp(self.th[0]), self.kuse*math.exp(-self.Tc/self.trel))

        # Conditional integral: remove the integral accumulated before/during
        # recovery without resetting the actual plant's hidden state.
        In = np.clip(self.rs[2]+(P.FD-F)*self.Tc, -20., 20.)
        ki = self.ki if self.Ioff == 0. else float(np.clip(self.ki*(In-self.Ioff)/In, 0., self.ki)) if abs(In)>1e-9 else 0.
        I0 = ki/self.ki*In-(P.FD-F)*self.Tc if self.ki else 0.
        nd = self.ndel
        tgt = log[n-1-nd,11] if n-1-nd >= 0 else log[0,2]
        q = np.array([log[j,11] if j>=0 else log[0,2] for j in range(n-nd,n)])
        st = np.array([self.rs[0],self.rs[1],I0,self.rs[3],y,F,tgt,
                       self.th[1],self.kuse,self.th[2],self.th[3],y,0.,0.,0.,0.])
        if self.motion:
            st = self._motion_state(st)
        self.dith_a = np.clip(self.Fdith/self.kuse,self.amin,self.amax) if self.Fdith else 0.
        dith = np.array([self._dither(n+i) for i in range(self.H)])
        D, KP = self.rule(self.kuse)
        prior = np.array([math.log(D),KP/1e-4,0.,0.])
        prev = self.x.copy()
        prev[2] = np.clip(prev[2]+prev[3]*1e-2*min(self.Tc,self.tu)/1e-3,XLO[2],XHI[2])
        lo, hi = XLO.copy(), XHI.copy()
        lo[0], hi[0] = max(LDMIN,prev[0]-self.rate_slow),min(LDMAX,prev[0]+self.rate_fast)
        lo[1], hi[1] = max(0.,prev[1]-self.rateK_fast),min(10.,prev[1]+self.rateK_slow)
        if not self.reference:
            lo[2:] = 0.; hi[2:] = 0.; prev[2:] = 0.
        x = np.clip(prev, lo, hi)
        info = np.zeros(4)
        cost = _solve_ord(x,self.iters,self.H,self.per,st,q,dith,self.pr,
                          self.scen,prev,prior,self.rw,lo,hi,self.Jc,self.rprev,info)
        self.x = x
        self.g = (self.M,math.exp(x[0]),0.,ki,x[1]*1e-4,x[2]*1e-3+dith[0])
        self._record(o, started, cost, info[0])
        return self.g
