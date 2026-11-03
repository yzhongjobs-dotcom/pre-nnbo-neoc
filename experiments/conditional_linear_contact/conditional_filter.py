"""条件线性接触响应滤波器与 V0 控制接口。

本文件只在新实验目录中工作；被继承的控制器、场景和被控对象均保持只读。
状态为 [k, D, r, a]。每个 8 ms 区间使用一个因果测量：

    z_j = F_{j+1} - a_f^8 F_j
        = [Y_j, W0, W1, W2] x_j + eps_j .

Y_j 由 common._build 的 enc_y 分支计算，W0/W1/W2 使用区间起点相对时间
tau_m = 0,...,7 ms。IMM 的四个模式是运动 {q_s,q_m} 与刚度
{q_lnk0,q_lnk1} 的笛卡尔积。
"""
from __future__ import annotations

import math
import os
import sys
from collections import deque

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "motion_aware_scan"))
sys.path.insert(0, os.path.join(HERE, "..", "method_redesign"))

import common  # noqa: E402
from common import K_MAX, K_MIN  # noqa: E402
from v1_servo_id import ServoID  # noqa: E402


def _psd(P: np.ndarray, floor: float = 1e-18) -> np.ndarray:
    """对称化并将极小的数值负特征值截到正半定。"""
    P = (np.asarray(P, dtype=float) + np.asarray(P, dtype=float).T) * 0.5
    if not np.all(np.isfinite(P)):
        raise FloatingPointError('non-finite covariance')
    # Joseph form is PSD up to round-off.  Keep the hot path algebraic; use an
    # eigendecomposition only when a grossly negative diagonal is observed.
    d = np.diag(P).copy()
    d[d < floor] = floor
    np.fill_diagonal(P, d)
    if np.min(d) >= floor:
        return P
    else:
        vals, vecs = np.linalg.eigh(P)
        vals = np.maximum(vals, floor)
        return (vecs * vals) @ vecs.T


class FourModeIMM:
    """四模式、四状态的标准 IMM。

    update_interval() 在测量 z_j 对应的 x_j 上做 Joseph 更新，然后用
    F_j 将每个模式推进到 x_{j+1}。这样 Delta-y 明确进入 F_j，k-D
    交叉协方差不会被丢掉。
    """

    def __init__(self, dt=0.008, q_smooth=1e-4, q_maneuver=1.,
                 p_stay_motion=.999, q_lnk0=1e-4, q_lnk1=1e-2,
                 p_stay_stiff=.98, R=.04 ** 2, r_scale=1., prob_floor=1e-12):
        self.dt = float(dt)
        self.q_motion = (float(q_smooth), float(q_maneuver))
        self.q_lnk = (float(q_lnk0), float(q_lnk1))
        self.R = float(R) * float(r_scale)
        self.prob_floor = float(prob_floor)
        self.modes = [(m, s) for m in range(2) for s in range(2)]
        tm = np.array([[p_stay_motion, 1. - p_stay_motion],
                       [1. - p_stay_motion, p_stay_motion]], dtype=float)
        ts = np.array([[p_stay_stiff, 1. - p_stay_stiff],
                       [1. - p_stay_stiff, p_stay_stiff]], dtype=float)
        self.transition = np.array([[tm[mi, mj] * ts[si, sj]
                                     for mj, sj in self.modes]
                                    for mi, si in self.modes], dtype=float)
        self.A_ca = np.array([[1., self.dt, self.dt ** 2 / 2.],
                              [0., 1., self.dt], [0., 0., 1.]])
        self.Q_ca = np.array([[self.dt ** 5 / 20., self.dt ** 4 / 8., self.dt ** 3 / 6.],
                              [self.dt ** 4 / 8., self.dt ** 3 / 3., self.dt ** 2 / 2.],
                              [self.dt ** 3 / 6., self.dt ** 2 / 2., self.dt]])
        self.reset()

    def reset(self, k0=3e4, D0=0.):
        self.initialized = False
        self.mu = np.full(4, .25, dtype=float)
        self.x_models = np.tile([float(k0), float(D0), 0., 0.], (4, 1))
        self.P_models = np.tile(np.diag([max(float(k0), 1.) ** 2, 1., 1., 1.]), (4, 1, 1))
        self.x = self.x_models[0].copy()
        self.P = self.P_models[0].copy()
        self.last_nis = np.nan
        self.last_innovation = np.nan
        self.last_S = np.nan
        self.last_mu = self.mu.copy()
        self.last_mode_states = self.x_models.copy()

    def freeze_contact(self, inflate=4.):
        """失联时保持 k，重置 D/r/a 并放大 k 的协方差。"""
        if not self.initialized:
            return
        for i in range(4):
            k = self.x_models[i, 0]
            P = self.P_models[i].copy()
            P[0, 0] *= float(inflate)
            self.x_models[i] = [k, 0., 0., 0.]
            self.P_models[i] = _psd(P)
        self.x = self.x_models @ self.mu
        d = self.x_models - self.x
        self.P = _psd(np.einsum('i,ijk->jk', self.mu, self.P_models)
                      + np.einsum('i,ij,ik->jk', self.mu, d, d))

    def recontact(self):
        if not self.initialized:
            return
        for i in range(4):
            self.x_models[i, 1:] = 0.
            self.P_models[i, 1:, 1:] = np.diag([1., 1., 1.])
        self.x = self.x_models @ self.mu

    def _Q(self, mode, k_ref):
        motion, stiff = self.modes[mode]
        k_ref = float(np.clip(abs(k_ref), K_MIN, K_MAX))
        qF = k_ref * k_ref * self.q_motion[motion]
        qk = (k_ref * self.q_lnk[stiff]) ** 2 * self.dt
        Q = np.zeros((4, 4), dtype=float)
        Q[0, 0] = qk
        Q[1:, 1:] = qF * self.Q_ca
        return Q

    def _combine(self):
        self.x = self.mu @ self.x_models
        d = self.x_models - self.x
        self.P = _psd(np.einsum('i,ijk->jk', self.mu, self.P_models)
                      + np.einsum('i,ij,ik->jk', self.mu, d, d))

    def initialize_if_needed(self, k0, D0):
        if self.initialized:
            return
        self.reset(k0, D0)
        self.initialized = True

    def update_interval(self, z, H, F, k_ref):
        """混合、并行 Joseph 更新、似然合并和时间传播。"""
        if (not np.isfinite(z) or not np.all(np.isfinite(H)) or not np.all(np.isfinite(F))
                or np.max(np.abs(F)) > 1e8 or np.max(np.abs(H)) > 1e8):
            return False
        if not self.initialized:
            self.initialize_if_needed(k_ref, 0.)

        joint = self.mu[:, None] * self.transition
        prior_mu = joint.sum(axis=0)
        mixing = joint / np.maximum(prior_mu[None, :], np.finfo(float).tiny)
        mixed_x = mixing.T @ self.x_models
        mixed_P = np.empty_like(self.P_models)
        for j in range(4):
            d = self.x_models - mixed_x[j]
            mixed_P[j] = _psd(np.einsum('i,ijk->jk', mixing[:, j], self.P_models)
                              + np.einsum('i,ij,ik->jk', mixing[:, j], d, d))

        post_x = np.empty_like(self.x_models)
        post_P = np.empty_like(self.P_models)
        log_lik = np.empty(4)
        nis = np.empty(4)
        for j in range(4):
            xp = mixed_x[j]
            try:
                Pp = _psd(mixed_P[j])
                with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
                    S0 = float(H @ Pp @ H + self.R)
                    innov = float(z - H @ xp)
                    if not np.isfinite(S0) or S0 <= 0. or S0 > 1e30 or not np.isfinite(innov):
                        return False
                    K = (Pp @ H) / max(S0, np.finfo(float).tiny)
                    xu = xp + K * innov
                    I_KH = np.eye(4) - np.outer(K, H)
                    Praw = I_KH @ Pp @ I_KH.T + self.R * np.outer(K, K)
                    # The measurement is at x_j; only then propagate to x_{j+1}.
                    xnext = F @ xu
                    Q = self._Q(j, max(xu[0], k_ref))
                    Pnext_raw = F @ Praw @ F.T + Q
                if (not np.all(np.isfinite(xu)) or not np.all(np.isfinite(Praw))
                        or not np.all(np.isfinite(xnext)) or not np.all(np.isfinite(Pnext_raw))):
                    return False
                Pu = _psd(Praw)
                Pnext = _psd(Pnext_raw)
            except (FloatingPointError, ValueError, np.linalg.LinAlgError):
                return False
            post_x[j], post_P[j] = xnext, Pnext
            nis[j] = min(float(innov * innov / max(S0, np.finfo(float).tiny)), 1e30)
            log_lik[j] = -0.5 * (math.log(2. * math.pi) + math.log(max(S0, np.finfo(float).tiny)) + nis[j])

        lw = np.log(np.maximum(prior_mu, np.finfo(float).tiny)) + log_lik
        lw -= np.max(lw)
        w = np.exp(lw)
        w /= max(float(w.sum()), np.finfo(float).tiny)
        self.mu = self.prob_floor + (1. - 4. * self.prob_floor) * w
        self.x_models, self.P_models = post_x, post_P
        self._combine()
        self.last_nis = float(self.mu @ nis)
        self.last_innovation = float(np.sum(self.mu * np.array([
            z - H @ (mixing[:, j] @ self.x_models) for j in range(4)])))
        self.last_S = float(self.R + max(self.P[0, 0], 0.) * H[0] ** 2)
        self.last_mu = self.mu.copy()
        self.last_mode_states = self.x_models.copy()
        return True


def interval_measurement(log, j, ns, eb, afc, tau):
    """严格构造一拍的 Y、W 和 z；Y 与 common._build(enc=True) 同源。"""
    if j < 0 or j + 1 >= len(log):
        return None
    if min(float(log[j, 5]), float(log[j + 1, 5])) <= .5:
        return None
    # _build 的 enc_y 分支只用端点和 eb；用单行调用避免重复实现插值细节。
    X, Y, Z, keep = common._build(log, j, 1, 0, int(ns), float(eb), float(afc),
                                  1, float(log[j, 0]), float(log[j, 2]),
                                  float(ns * 1e-3), True)
    if not bool(keep[0]):
        return None
    a = 1. - float(afc)
    mm = np.arange(1, int(ns) + 1, dtype=float)
    w = float(afc) * a ** (int(ns) - mm)
    tau_m = (mm - 1.) * 1e-3
    W0 = float(np.sum(w))
    W1 = float(w @ tau_m)
    W2 = float(w @ (tau_m ** 2 / 2.))
    H = np.array([float(Y[0]), W0, W1, W2], dtype=float)
    z = float(Z[0])
    return H, z


class ConditionalLinearContact(ServoID):
    """V0：以条件线性 IMM 后验替换刚度估计源，其余控制接口继承 V1。"""

    name = "conditional-linear-v0"

    def __init__(self, q_lnk0=1e-4, q_lnk1=1e-2, p_stay_stiff=.98,
                 gain_rel_max=math.log(1.2), consistency_max=2.68,
                 consistency_window=8, switch_time=1.0, k_cov_inflate=4.,
                 **kw):
        super().__init__(**kw)
        self.conditional_params = dict(q_lnk0=float(q_lnk0), q_lnk1=float(q_lnk1),
                                       p_stay_stiff=float(p_stay_stiff),
                                       gain_rel_max=float(gain_rel_max),
                                       consistency_max=float(consistency_max),
                                       consistency_window=int(consistency_window),
                                       switch_time=float(switch_time),
                                       k_cov_inflate=float(k_cov_inflate))
        self.params.update(self.conditional_params)

    def reset(self, platform):
        super().reset(platform)
        p = self.params
        self.cl_filter = FourModeIMM(
            dt=self.Tc, q_smooth=p['q_smooth'], q_maneuver=p['q_maneuver'],
            p_stay_motion=p['p_stay'], q_lnk0=p['q_lnk0'], q_lnk1=p['q_lnk1'],
            p_stay_stiff=p['p_stay_stiff'], R=p['noise'] ** 2,
            r_scale=p.get('r_scale', 1.))
        self._next_interval = 1
        self._contact_lost = False
        self._filter_ready = False
        self._handoff_pending = False
        self._nis_history = deque(maxlen=max(1, int(p['consistency_window'])))
        self.mode_history = []
        self.filter_trace = []
        self.accept_history = []
        self.accepted_count = 0
        self.valid_count = 0
        self.filter_initialized_count = 0

    def _base_fit(self, o, n):
        p = self.params
        fit = common.window_fit(o.log, n, p['nw'], p['degree'], p['noise'],
                                self.ns, self.eb, self.afc, self.ndel, p['enc_y'])
        self.fit = fit  # V0 保留原有 drift_pred 的多项式前馈接口。
        return fit

    def _update_filter(self, o, n):
        Fnow = float(o.log[n - 1, 5])
        if Fnow <= .5:
            if not self._contact_lost:
                self.cl_filter.freeze_contact(self.params['k_cov_inflate'])
            self._contact_lost = True
            self._next_interval = max(self._next_interval, n - 1)
            return
        if self._contact_lost:
            self.cl_filter.recontact()
            self._contact_lost = False
            self._next_interval = max(self._next_interval, n - 1)

        # ServoID._identify has already run in this call, so self.tau is current V1 tau.
        last = n - 2
        if last < self._next_interval:
            return
        for j in range(self._next_interval, last + 1):
            item = interval_measurement(o.log, j, self.ns, self.eb, self.afc, self.tau)
            if item is None:
                continue
            H, z = item
            if not self.cl_filter.initialized:
                fit0 = common.window_fit(o.log, min(n, j + 2), self.params['nw'],
                                         self.params['degree'], self.params['noise'],
                                         self.ns, self.eb, self.afc, self.ndel,
                                         self.params['enc_y'])
                k0 = fit0['k'] if np.isfinite(fit0.get('k', np.nan)) else self.k_seek
                k0 = float(np.clip(k0, K_MIN, K_MAX))
                D0 = float(o.log[j, 5])
                self.cl_filter.initialize_if_needed(k0, D0)
                self.filter_initialized_count += 1
            dy = float(o.log[j + 1, 2] - o.log[j, 2])
            Fmat = np.array([[1., 0., 0., 0.],
                             [dy, 1., self.Tc, self.Tc ** 2 / 2.],
                             [0., 0., 1., self.Tc],
                             [0., 0., 0., 1.]], dtype=float)
            self.cl_filter.update_interval(z, H, Fmat,
                                           max(float(self.cl_filter.x[0]), K_MIN))
            if np.isfinite(self.cl_filter.last_nis):
                self._nis_history.append(float(self.cl_filter.last_nis))
            score = float(np.mean(self._nis_history)) if self._nis_history else np.inf
            self.valid_count += 1
            self.mode_history.append(dict(time=float(o.log[j + 1, 0]),
                                          mu=self.cl_filter.mu.tolist(),
                                          nis=float(self.cl_filter.last_nis),
                                          consistency=score))
            self.filter_trace.append(dict(time=float(o.log[j + 1, 0]),
                                          k=float(self.cl_filter.x[0]),
                                          sigma=float(math.sqrt(max(self.cl_filter.P[0, 0], 0.))),
                                          consistency=score,
                                          mu=self.cl_filter.mu.tolist()))
        self._next_interval = last + 1

    def _stiffness_est(self, o, n):
        p = self.params
        # Keep the original fit for the unchanged V0 drift feed-forward path.
        base_fit = self._base_fit(o, n)
        self._update_filter(o, n)
        t = float(o.log[n - 1, 0])
        if t < p['switch_time'] or not self.cl_filter.initialized:
            # Exactly the V1+Proposed path during startup; the new filter is parallel only.
            return super()._stiffness_est(o, n)
        if not self._filter_ready:
            self._filter_ready = True
            self._handoff_pending = True
        if self._handoff_pending:
            # Preserve the already-running V1 feedback stiffness at the hand-off instant.
            self._handoff_pending = False
            return self.kbar, False

        k_hat = float(np.clip(self.cl_filter.x[0], K_MIN, K_MAX))
        sigma = float(math.sqrt(max(self.cl_filter.P[0, 0], 0.)))
        rel = sigma / max(k_hat, K_MIN)
        consistency = float(np.mean(self._nis_history)) if self._nis_history else np.inf
        precise = rel <= p['rel_max']
        consistent = consistency <= p['consistency_max']
        self.k_hat, self.rel = k_hat, rel
        self.k_check = k_hat
        self.valid_count += 0
        if rel <= p['gain_rel_max'] and consistent:
            self.k_good = k_hat
            self.f_set = max(float(o.log[n - 1, 5]), .3)
            self.kg_rel = rel
            self.accepted_count += 1
            accepted_gain = True
        else:
            accepted_gain = False
        self.accept_history.append(dict(time=t, accepted=accepted_gain, precise=precise,
                                        consistent=consistent, rel=rel,
                                        consistency=consistency, k=k_hat,
                                        sigma=sigma))
        # Feedback rule is the MotionAware rule with the IMM posterior and the
        # innovation-consistency average replacing the window chi2.
        if not np.isfinite(k_hat) or not np.isfinite(sigma):
            return self.kbar, False
        if self.first:
            if precise and consistent:
                self.first = False
                cand = np.clip(k_hat + p['c'] * sigma, K_MIN, K_MAX)
                return cand, True
            return self.kbar, False
        if precise and consistent:
            cand = np.clip(k_hat + p['c'] * sigma, K_MIN, K_MAX)
            return max(cand, (1. - p['rho_d']) * self.kbar), True
        if precise:
            # Inconsistent windows can only raise the conservative feedback stiffness.
            return max(k_hat, self.kbar), False
        return self.kbar, False


class MeasuredOnly:
    """只向控制器暴露 run_joint.py 约定的测量列。"""
    direct = True

    def __init__(self, inner):
        self.inner = inner

    def reset(self, platform):
        self.inner.reset(dict(common.P.PLATFORMS['PC125']))

    def __call__(self, obs):
        a = np.full((obs.n, 16), np.nan)
        a[:, [0, 1, 2, 5, 11, 12, 13]] = obs.log[:obs.n, [0, 1, 2, 5, 11, 12, 13]]
        return self.inner(common.P.Obs(a, obs.n))
