"""V4b: direct FINAL motion increments plus C1 in residual coordinates.

Only columns 0,2,5,11 are read. Future motion increments are held at zero.
No contact truth or future observation is used. See v4b_evaluate.py.
"""
import math
import time
import numpy as np
from v4_architecture import ContactUpdate, ServoID, CE, P, common, variants
from candidate1_model import _sync, XLO, XHI, LDMIN, LDMAX
from contact_update import profile_contact
from v4b_kernels import _solve_ord


class ResidualC1(ContactUpdate):
    W = ServoID.W

    def reset(self, platform):
        super().reset(platform)
        self.tau_hat, self.ndel_hat = self.tau, self.ndel
        self.known_ff_increment = 0.

    def _apply(self, tau, ndel):
        self.tau, self.ndel = tau, int(ndel)
        self.eb = math.exp(-.001/tau)
        self.pr[5] = self.eb

    def __call__(self, o):
        started = time.perf_counter()
        n, log = o.n, o.log
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
        st = np.r_[st, self.known_ff_increment]
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

class V4b:
    direct = True

    def __init__(self, motion=True, startup_prior=30000., gate_settling=False):
        self.motion = motion
        self.startup_prior = float(startup_prior)
        self.gate_settling = bool(gate_settling)

    def reset(self, platform):
        self.platform = dict(platform)
        self.Tc, self.VMAX = platform['Tc'], platform['VMAX']
        self.c1 = ResidualC1(k0=self.startup_prior)
        self.c1.reset(platform)
        self.observer = ServoID(single_model=False, **CE.FINAL)
        self.observer.reset(platform)
        self.ff_offset = self.ff_servo = 0.
        self.offsets, self.responses = [], []
        self.virtual = None
        self.exec_rs = None
        self.t_call = []
        self.counts = dict(ticks=0, motion_ticks=0, gain_ticks=0,
                           recovery_ticks=0, recontacts=0, rate_clip_ticks=0,
                           contact_without_gain_ticks=0, nonzero_ff_ticks=0,
                           ready_in_recovery_ticks=0)
        self.first_motion_time = None
        self.max_ff_increment = 0.
        self.max_command_identity_error = 0.
        self.sum_rebase_abs = 0.
        self.sum_ff_abs = self.sum_opt_abs = 0.
        self.max_residual_velocity = 0.

    def _reset_motion(self):
        # Preserve accepted stiffness as in FINAL loss handling; reset IMM level
        # and increment anchor. Clearing k_good would require a fresh excitation
        # window after every touch, introducing another startup gate.
        self.observer._clear_motion()

    def __call__(self, o):
        started = time.perf_counter()
        k = o.n-1                      # completed ticks; row k is provisional
        log = o.log
        ob, c = self.observer, self.c1
        F = float(log[k, 5])
        if self.virtual is None:
            self.virtual = np.zeros_like(log)
        if k:
            # The command sent at k-1-ndel acted over the previous interval.
            j = k-1-c.ndel
            target = self.offsets[j] if j >= 0 else 0.
            rho = math.exp(-self.Tc/c.tau)
            self.ff_servo = target+(self.ff_servo-target)*rho
        self.responses.append(self.ff_servo)
        if k > 60 and k % 25 == 0:
            ServoID._identify(ob, log, k+1)
            c.tau_hat, c.ndel_hat = ob.tau_hat, ob.ndel_hat
            c._apply(ob.tau_hat, ob.ndel_hat)
        # Rebuild only the newly completed row. Both physical and residual logs
        # retain the same causal force samples; no copied truth columns.
        if k:
            j = k-1
            self.virtual[j, 0] = log[j, 0]
            self.virtual[j, 2] = log[j, 2]-self.responses[j]
            self.virtual[j, 5] = log[j, 5]
            self.virtual[j, 11] = log[j, 11]-self.offsets[j]
        # Motion estimation is allowed whenever force is in contact, including
        # C1's 12-tick settling phase. Output resumes at F>=0.5, independently of C1 settling.
        # Thus observer readiness is independent of the C1 recovery-mode flag;
        # sufficient motion-measurement information is still required.
        ready = False
        can_ff = (not c.mode) if self.gate_settling else F >= .5
        increment = 0.
        if self.motion and F >= .3:
            ob.kbar, _ = ob.stiffness(o, o.n)
            q = ob.drift_pred(o, o.n)
            ready = np.isfinite(ob.k_good) and ob.k_good > 0 and ob.tracker.initialized
            if ready and q is not None and ob.q_prev is not None and can_ff:
                # Exact FINAL tuple-path conversion: dq=k_good*Delta eta_p,
                # -dq/k_good is added directly to the sent position command.
                dq = ob.k_good*(q[1]-ob.q_prev)
                increment = -dq/ob.k_good
            if q is not None:
                ob.q_prev = q[1] if ready else None
        else:
            self._reset_motion()
        if not can_ff:
            # Disable increments, keep the position offset continuous at loss.
            increment = 0.
            ob.q_prev = None
        self.ff_offset += increment
        c.known_ff_increment = increment
        oldmode = c.mode
        g = c(P.Obs(self.virtual, k))
        if oldmode == 1 and c.mode == 0:
            self._reset_motion()
            self.counts['recontacts'] += 1
        # C1 gains are chosen from the completed log, as in its gain interface.
        # The admittance itself consumes the fresh force, as plant.advance does.
        if self.exec_rs is None:
            y0 = float(log[k, 2])
            self.exec_rs = np.array([y0, 0., 0., y0])
        _sync(log, k, k+1, self.exec_rs, *g[:4], self.Tc, self.VMAX)
        base = self.exec_rs[0]+g[4]*(P.FD-F)+g[5]
        requested = base+self.ff_offset
        sent = float(log[k, 11])
        command = float(np.clip(requested, sent-self.VMAX*self.Tc, sent+self.VMAX*self.Tc))
        # FINAL integrates from the actually sent command. Mirror that anchor:
        # the unsent part of an FF offset must not survive a physical slew cap.
        if self.motion:
            previous_offset = self.ff_offset
            self.ff_offset = command-base
            self.sum_rebase_abs += abs(self.ff_offset-previous_offset)
            self.max_command_identity_error = max(self.max_command_identity_error,
                                                  abs(base+self.ff_offset-command))
        self.offsets.append(self.ff_offset)
        self.counts['ticks'] += 1
        self.counts['motion_ticks'] += int(ready and can_ff)
        self.counts['gain_ticks'] += int(np.isfinite(ob.k_good))
        self.counts['nonzero_ff_ticks'] += int(abs(increment)>1e-12)
        self.counts['ready_in_recovery_ticks'] += int(ready and F>=.5 and bool(c.mode))
        self.counts['recovery_ticks'] += int(bool(c.mode))
        self.counts['rate_clip_ticks'] += int(abs(command-requested)>1e-12)
        self.counts['contact_without_gain_ticks'] += int(F >= .5 and not ready)
        if ready and can_ff and self.first_motion_time is None:
            self.first_motion_time = float(log[k, 0])
        self.max_ff_increment = max(self.max_ff_increment, abs(increment))
        self.sum_ff_abs += abs(increment)
        self.sum_opt_abs += abs(g[5]-c._dither(k))
        self.max_residual_velocity = max(self.max_residual_velocity, abs(c.th[2]))
        self.t_call.append(time.perf_counter()-started)
        return command


class PriorC1(ResidualC1):
    """C1 estimator/recovery copied unchanged; normal decisions use its prior.

    No GN solver or prediction rollouts are executed. The copy below differs
    only at the decision block; this avoids the legacy iters=0 frozen-gain bug.
    """
    def __call__(self, o):
        started = time.perf_counter()
        n, log = o.n, o.log
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
        # No objective/Jacobian/rollout evaluation in the noOpt path.
        x = prior.copy()
        info = np.zeros(4)
        cost = 0.
        self.x = x
        self.g = (self.M,D,0.,ki,KP,dith[0])
        self._record(o, started, cost, info[0])
        return self.g

class V4bNoOpt(V4b):
    def reset(self, platform):
        super().reset(platform)
        self.c1 = PriorC1(k0=self.startup_prior)
        self.c1.reset(platform)

