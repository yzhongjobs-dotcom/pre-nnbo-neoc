"""Step 3: force-domain MPC; estimation/supervision inherited verbatim from FINAL.

z_i = k_good * delta_u_i [N]; cost sum (F_filtered-FD)^2 +
R*(kbar/k_good*z_i)^2, last force weight multiplied by terminal.
The force-domain Hessian depends only on k_good/kbar. A command selected now
acts after one complete period. Force samples precede servo updates, matching
the reference linearisation, with exact exponential servo propagation.
"""
from pathlib import Path
import sys
import time
import os
from functools import lru_cache
import numpy as np
from scipy import sparse
from scipy.linalg import cho_factor, cho_solve
import osqp

sys.dont_write_bytecode = True
os.environ.setdefault('NUMBA_CACHE_DIR', str(Path(__file__).resolve().parent/'.numba_cache'))

SCAN = Path(__file__).resolve().parents[1] / 'motion_aware_scan'
sys.path.insert(0, str(SCAN))
import common
import cohort_eval as CE
import variants
from numba import njit

TC, TAU, FC, FMAX = .008, .019, 30., 15.
NS = 8
AF = 1-np.exp(-2*np.pi*FC*.001)
DEC = (1-AF)**NS
WEIGHTS = AF*(1-AF)**np.arange(NS-1, -1, -1)
AS = np.exp(-TC/TAU)
AM = np.exp(-np.arange(NS)*.001/TAU)
CY = float(WEIGHTS@AM)
CU = float(WEIGHTS@(1-AM))


@njit(cache=True)
def oracle_surface_corrections(p, x0, t, H, kh, h0, ht):
    """True prescribed scanner/surface trajectory, including the plant x-servo.
    This oracle may read scenario parameters, unlike either deployable MPC.
    """
    rex, fex = np.zeros(H), np.zeros(H)
    x, ff = x0, 0.
    for it in range(H*16):
        dt = it*.0005
        xc = common.P.X0+p[5]*max(t+dt-TC, 0.)
        x += .0005/TAU*(xc-x)
        hh = common.P.surface(x, t+dt, p)[0]
        if it % 2 == 0:
            ff = (1-AF)*ff+AF*kh*(hh-h0-ht*dt)
        if it % 16 == 15:
            i = it//16
            fex[i] = ff
            dt_end = (i+1)*TC
            rex[i] = kh*(common.P.surface(x, t+dt_end, p)[0]-h0-ht*dt_end)
    return rex, fex


@lru_cache(None)
def prediction_matrices(H):
    """Baseline map from [khat*(sent-y), model_F, model_raw_F,
    khat*eta_velocity, khat*eta_acceleration]; response to future z.
    eta=-surface, so positive eta increases force at fixed encoder position.
    """
    def rollout(v, z):
        U, filt, raw, vel, acc = v
        Y = 0.
        ff, rr = [], []
        for i in range(H):
            tm = i*TC+np.arange(NS)*.001
            filt = DEC*filt+CY*Y+CU*U+float(WEIGHTS@(raw+vel*tm+.5*acc*tm**2))
            Y = AS*Y+(1-AS)*U
            U += z[i]
            te = (i+1)*TC
            ff.append(filt)
            rr.append(raw+Y+vel*te+.5*acc*te**2)
        return np.array(ff), np.array(rr)
    BF, BR = [], []
    for j in range(5):
        v = np.eye(5)[j]
        f, r = rollout(v, np.zeros(H)); BF.append(f); BR.append(r)
    GF, GR = [], []
    for j in range(H):
        f, r = rollout(np.zeros(5), np.eye(H)[j]); GF.append(f); GR.append(r)
    return np.array(BF).T, np.array(BR).T, np.array(GF).T, np.array(GR).T


class ForceMPC:
    def __init__(self, H=15, qr=.1, terminal=1., constrained=True):
        self.H, self.qr, self.terminal = int(H), float(qr), float(terminal)
        self.constrained = constrained
        self.BF, self.BR, self.GF, self.GR = prediction_matrices(self.H)
        self.weights = np.ones(self.H); self.weights[-1] = self.terminal
        self.GTG = self.GF.T@(self.weights[:, None]*self.GF)
        self.qp = None
        self.times = []
        self.statuses = {}
        self.infeasible = 0
        self.max_violation = 0.
        self._iu = np.triu_indices(self.H)
        # OSQP uses column-major upper triangular CSC storage.
        self._order = np.lexsort((self._iu[0], self._iu[1]))
        self._warm = np.zeros(self.H)
        self._A = sparse.csc_matrix(np.vstack([np.eye(self.H), self.GR]))

    def solve(self, state, bias, khat, kbar, vmax=.25, raw_extra=None, filt_extra=None):
        t0 = time.perf_counter()
        f = self.BF@state+bias
        raw = self.BR@state
        if raw_extra is not None: raw = raw+raw_extra
        if filt_extra is not None: f = f+filt_extra
        ratio = kbar/khat
        mat = self.GTG+np.eye(self.H)*(ratio*ratio/self.qr)
        q = self.GF.T@(self.weights*(f-common.FD))
        if not self.constrained:
            sol = -cho_solve(cho_factor(mat, check_finite=False), q, check_finite=False)
            status = 'analytic'
        else:
            P = 2*mat
            px = P[self._iu][self._order]
            A = self._A
            limit = vmax*TC*khat
            lo = np.r_[np.full(self.H, -limit), np.full(self.H, -np.inf)]
            hi = np.r_[np.full(self.H, limit), (FMAX-bias)/ratio-raw]
            if self.qp is None:
                pp = sparse.csc_matrix((px, (self._iu[0][self._order], self._iu[1][self._order])), shape=P.shape)
                self.qp = osqp.OSQP()
                self.qp.setup(P=pp, q=2*q, A=A, l=lo, u=hi, verbose=False,
                              eps_abs=1e-5, eps_rel=1e-5, max_iter=2000,
                              polishing=False, adaptive_rho=True, warm_starting=True)
            else:
                self.qp.update(Px=px, q=2*q, l=lo, u=hi)
            self.qp.warm_start(x=self._warm)
            ans = self.qp.solve(raise_error=False)
            status = ans.info.status
            if ans.info.status_val in (1, 2) and ans.x is not None and np.all(np.isfinite(ans.x)):
                sol = ans.x
                self.max_violation = max(self.max_violation, float(np.max(A@sol-hi)))
                self._warm = np.r_[sol[1:], 0.]
            else:
                # A hard cap may already be impossible in the delayed first interval.
                # Keep the hard problem, report failure, and retract at the speed limit.
                self.infeasible += 1
                sol = np.full(self.H, -limit)
                self._warm[:] = 0.
        elapsed = time.perf_counter()-t0
        self.times.append(elapsed)
        self.statuses[status] = self.statuses.get(status, 0)+1
        return float(sol[0]/khat)


class EstimationAwareMPC(variants.MotionAwareIMM):
    """super().step owns estimation, IMM, dither, contact detection and recovery.

    Only admittance.increment and the returned feed-forward are suppressed.
    The same drift_pred still updates the unmodified IMM once per contact tick.
    Its state is relative eta=-s, with an unobservable constant; the force anchor
    fixes that gauge. A constant measured-minus-model filtered-force residual
    is added to all horizon predictions (offset-free output disturbance).
    """
    def __init__(self, H=15, qr=.1, terminal=1., constrained=True, oracle_scene=None):
        super().__init__(**dict(CE.FINAL))
        self.mpc_params = dict(H=H, qr=qr, terminal=terminal, constrained=constrained)
        self.oracle_scene = oracle_scene
        self.name = 'MPC-oracle' if oracle_scene is not None else ('MPC-C' if constrained else 'MPC-U')

    def reset(self, platform):
        super().reset(platform)
        if self.Tc != TC or self.ndel != 1 or self.tau != TAU or self.fc != FC:
            raise ValueError('prototype prediction is defined for nominal PC125 only')
        self.adm.increment = lambda e: 0.
        self.mpc = ForceMPC(**self.mpc_params)
        self._anchor = None
        self._pred_f = None
        self._last_k = None
        self._controlled = False

    def stiffness(self, o, n):
        self._controlled = True
        return super().stiffness(o, n)

    def drift_pred(self, o, n):
        super().drift_pred(o, n)
        return None  # tracking preserved; old feed-forward replaced by MPC

    def step(self, o):
        self._controlled = False
        inherited = super().step(o)
        if not self._controlled:  # identical inherited approach and recontact supervision
            self._anchor = None; self._pred_f = None
            return inherited
        j = o.n-1
        t, y, F, sent = o.log[j, [0, 2, 5, 11]]
        kh = self.k_good if np.isfinite(self.k_good) else self.kbar
        kb = self.kbar
        eta = self.tracker.x.copy() if self.tracker.initialized else np.zeros(3)
        # The interval-mean observation precedes the encoder endpoint. Shift the
        # CA state by the same nominal interval timestamp, without servo inversion.
        lag = (1-self._encoder_fraction)*TC
        eta = np.array([eta[0]+lag*eta[1]+.5*lag*lag*eta[2], eta[1]+lag*eta[2], eta[2]])
        if self.oracle_scene is not None:
            kh = kb = max(float(o.log[j, 8]), common.K_MIN)
            h, ht, _ = common.P.surface(float(o.log[j, 1]), t, self.oracle_scene)
            eta = np.array([h, ht, 0.])  # plant indentation=y+h, eta=-s=h
        # Reanchor on gain changes: physical force-origin uncertainty must not
        # produce a jump when stiffness or the equivalent-motion gauge changes.
        if self._anchor is None or kh != self._last_k:
            fin = F if j == 0 else (F-DEC*o.log[j-1, 5])/(1-DEC)
            yp = y if j == 0 else o.log[j-1, 2]
            ybar = yp+self._encoder_fraction*(y-yp)
            f0 = float(fin+kh*((y-ybar)+lag*eta[1]))
            self._anchor = (y, eta[0], f0)
            self._pred_f = float(F)
        ya, ea, fa = self._anchor
        f0 = fa+kh*(y-ya+eta[0]-ea)
        bias = F-self._pred_f
        state = np.array([kh*(sent-y), self._pred_f, f0, kh*eta[1], kh*eta[2]])
        rex = fex = None
        if self.oracle_scene is not None:
            # Exact future surface trajectory, constant current true tangent k.
            # Future stiffness is not read from another controller's trajectory.
            rex, fex = oracle_surface_corrections(np.asarray(self.oracle_scene),
                                                  float(o.log[j, 1]), float(t), self.mpc.H, kh, h, ht)
        du = self.mpc.solve(state, bias, kh, kb, self.VMAX, rex, fex)
        # inherited contains ONLY unchanged external force-amplitude dither.
        result = inherited+du
        # Predict the next filtered measurement under the command pipeline.
        self._pred_f = float(self.mpc.BF[0]@state+(0. if fex is None else fex[0]))
        self._last_k = kh
        return result
