"""Candidate 2 (learning route, YOPO-style division of information): privileged teacher ->
causal student by DAgger, then task-cost fine-tuning by evolution strategies.

DIVISION OF INFORMATION
  Training only (privileged): the true footprint stiffness k_true (log column 8, no delay), used
  by the teacher (truth-stiffness rule, oracle B3b at 0 s latency) to label the student's own
  states; and the plant's task cost J (scenarios.metrics on the TRUE contact force), the
  stage-2 objective.
  Deployment (causal): the student reads only Obs-readable signals: filtered measured force
  (log col 5), encoder y (col 2) and scan x (col 1). It never reads k_true, log cols 3/4/6/7/8 or
  the scene parameters p. (y_cmd, col 11, and the applied gains, cols 9/10, are computed but
  masked, see OFF_FEATS: they copy the executed action and leak the teacher's k in DAgger.)

POLICY (class Cand2Policy, interface identical to the baselines: reset(pl), __call__(obs) -> g)
  Sensing front end (fixed): the EtR-A lock-in stiffness estimator and its 20.8 Hz force-
  targeted command dither, re-implemented in numba with the EtR-A T3-tuned settings (ncyc 2,
  f_dither 20 Hz -> 6 ticks/cycle, ndiff 3, rho 0.3, peak hold trel 0.1 s, F_dith 0.1 N,
  a in [2, 30] um).
  Observation (35 slots, 32 used, every control tick, causal window 0.4 s):
    0 log k_hat (raw lock-in)  1 log k_hat (peak hold)  2-13 ef = FD - F_filt at lags 0..11 ticks
    14-25 encoder increments y[i]-y[i-1] at lags 0..11   27 y(n-1)-y(n-13)  28 y(n-1)-y(n-51)
    29 leaky integral of ef (tau 0.3 s)  30 FD - mean F over 50 ticks  31 force std over 25 ticks
    32 scan speed from x;  26 (y_cmd - y), 33/34 (applied log D, log KPF): masked to 0.
  Features standardised with mean/sd from rollouts of the prior policy on train scenes, clip +-6.
  Network: MLP 35-64-64-5, tanh hidden, numpy/numba (torch is not installed); first-order
  filter (tau 0.03 s) on the network output.
  Output = bounded residual around the structured prior "stiffness rule on the causal lock-in
  estimate" (rule = EtR-A's tuned rule = the teacher's rule, wI 80, cP 2, ki 30, Dmin 100,
  M 0.3, K 0):
    k_inf = k_hold * exp(2.0 tanh z1)              (inferred contact stiffness)
    D = max(100, k_inf/80);  KPF = 2/k_inf * exp(0.7 tanh z2);  M = 0.3 exp(tanh z0);
    ki = 30 exp(tanh z3) (<= 30);  u_ff = KPF * 1 N * tanh z4 (+ the dither);  K = 0;
  all clipped to the plant bounds. Zero network output reproduces EtR-A exactly (`check`:
  max force difference 2e-11 N). Decision period: every control tick (8 ms); ablation (e)
  every 12 ticks (0.096 s), prior and residual both held between decisions.

TEACHER (privileged, stage 1)
  baselines_ref.TruthK at 0 s latency with the same rule (80, 2, 30, 100), M 0.3: train J
  T3 0.136, T1 0.011. Because the teacher is a known function of k_true and its M, ki and KPF
  form equal the prior's, imitating its (often saturated) D/KPF is replaced by regressing its
  privileged state: target of the stiffness head = (log k_true - log k_hold)/2 (k held during
  contact loss, as in TruthK); the other heads are not distilled (stay 0 until stage 2).
  Loss: expectile (tau 0.95, under-estimates of k penalised 19x) because under-estimating the
  stiffness on bone destabilises contact (MSE students were worse than EtR-A; candidate2_log.md).
  The rule tuning grid also contained a stiffness-map feed-forward teacher (TruthKFF, best train
  J T3 0.135 / T1 0.007, candidate2_teacher.npz); its u_ff could not be distilled stably.

TRAINING (per seed; seeds 0, 1, 2; ONE shared policy for T1 and T3)
  Scenes: train split only (scene rng seeds 1000-1999); each batch half at the two cell doses,
  half randomised (A 0.5-5 mm, ratio 10-30, T3 ramp 0.2-0.6 s); noise seeds >= 1e6.
  Stage 1 DAgger: N_DAG = 8 iterations x EP_DAG = 32 episodes (16 T3 + 16 T1); iteration i
  executes the teacher (plus the student's dither) with per-episode probability 0.5^i, otherwise
  the student; every 3rd tick labelled; all layers trained by Adam (lr 1e-3, batch 512, 25 epochs
  per iteration on the aggregate, 50 in the last).
  Stage 2 task-cost ES: separable CMA-ES (baselines_learned.SepCMA) on the last layer (W3, b3:
  325 parameters), lambda 16, sigma0 0.01, 50 generations, minibatch 8 scenes (4 T3 + 4 T1) per
  generation with common random numbers; fitness mean(J / (J_B1 + 0.01)), J_B1 = the cell's B1 on
  the same scene and noise. Starts from the selected stage-1 network.
  Model selection (both stages): dev cohorts (seed 2000, 20 scenes per cell, the two cell
  doses), score = mean over cells of J_dev/J_dev(B1), after every DAgger iteration and every 10
  ES generations; only checkpoints with peak <= 7 N in >= 95 % of the train cohorts (seed 1000,
  20 scenes per cell, the baselines' training cohorts) are eligible (the B1 selection rule).
  Budget per seed (full): ~256 DAgger + ~6.8k ES episodes + ~1.0k selection episodes; the exact
  counts are in candidate2_results.json (E27-style baseline: 14k episodes per seed).

ABLATIONS (same budgets, same initialisation)
  a  no teacher: stage 2 only, from the same random initialisation (W3 = 0, i.e. from EtR-A)
  b  teacher only: the selected stage-1 checkpoint of `full`, no stage 2
  c  gains only: u_ff residual forced to 0 (dither kept), both stages
  d  no force history: force-derived features (0-13, 29-31) zeroed AND the prior's stiffness
     centre is the constant K_C = 7000 N/m instead of k_hold (k_hold is force-derived); both
     stages. Remaining inputs: encoder y increments/displacements and scan speed.
  e  decision every 12 ticks (0.096 s) instead of every tick, both stages

DEPARTURES FROM THE EARLIER CANDIDATE-2 DESIGN (research/separate_candidate_design_and_
iteration.md sec. 4 and 10; ../dual_candidate_task_simulation)
  * plant: the shared numba PC125 benchmark plant (1-D normal admittance, 125 Hz, delay + lag,
    F/T filter) instead of the 2-D B0 Hertz force/moment model with four K/D gains and an
    energy tank; no stored-energy projection or gain-rate constraint (gains bounded only).
  * no MHE / shared estimate / spatial preview / task-weight inputs: the student sees raw causal
    history plus the EtR-A lock-in estimate.
  * the network output is a bounded residual around a structured prior (the stiffness rule on
    the lock-in estimate) instead of absolute impedance targets; the main learned quantity is the
    contact stiffness fed to that rule.
  * gradients do not pass through the plant (numba, non-differentiable): the task-cost stage is
    ES on the last layer instead of 600 analytic task-gradient steps on the whole network.
  * a privileged teacher + DAgger stage is added (the earlier design used simulator truth only
    for the rollout cost); this is the YOPO-style division asked for here. The distilled target
    is the teacher's privileged state (k_true), not its saturated action.
  * outputs on the benchmark's complete controller family (M, D, KPF, ki, u_ff; K = 0) instead
    of four K/D gains plus compensation; no u_ff stiffness feed-forward (it failed, see log).
  * one shared policy for both tasks; the baselines were tuned per cell.

Usage:
  python candidate2_learned.py check                 zero-output student == EtR-A
  python candidate2_learned.py teacher               teacher rule grid -> candidate2_teacher.npz
  python candidate2_learned.py train VARIANT SEED    VARIANT in full,a,c,d,e -> candidate2_VARIANT_seedS.npz
  python candidate2_learned.py test                  -> candidate2_results.json
  python candidate2_learned.py test2 | kb            secondary cell | known-behaviour + timing
"""
import math, json, sys, time, itertools, os
import numpy as np
from numba import njit
import plant as P, scenarios as S

# ------------------------------------------------------------------ constants / budgets
TC = 0.008
NF, NH, NO = 35, 64, 5
FORCE_FEATS = list(range(0, 14)) + [29, 30, 31]
OFF_FEATS = [26, 33, 34]      # own-action channels (y_cmd - y, applied D/KPF): masked, DAgger copycat leak
USE_FF = False                # structured k_hat feed-forward c_ff*phi: disabled (harmful, see log)
P0 = (1.5, 120., 0., 10., 1.15e-4)          # fixed prior of ablation d (no force -> no k_hat)
RULE_P = (80., 2., 30., 100.)               # wI, cP, ki, Dmin: EtR-A's tuned rule (both cells) = teacher's rule
M_R = 0.3
SPAN = (1.0, 2.0, 0.7, 1.0)                 # residual e-folds: M, contact stiffness k, extra KPF, ki
EXPECTILE = 0.95                            # conservative stiffness regression (under-estimates destabilise bone contact)
TEACHER_RULE = RULE_P + (-1.0, 0.0)         # TruthK, 0 s latency, same rule as the student's (no FF)
K_C = 7000.                                 # constant stiffness centre of ablation d (no k_hat)
FU = 1.0                                    # u residual = KPF * FU * tanh(z4)  [force-equivalent, N]
ZTAU = 0.03                                 # first-order filter on the network output [s]
TREF = 0.25                                 # k_ref = running min of k_hold after TREF s
# lock-in front end (EtR-A, T3-tuned estimator parameters from adaptive_results.json)
F_DITHER, NCYC, NDIFF, RHO, TREL, FDITH, AMIN, AMAX = 20., 2, 3, 0.3, 0.1, 0.1, 2e-6, 30e-6
K0, KMIN, KMAX = 2e4, 300., 3e5
# training budgets
N_DAG, EP_DAG, EPOCHS = 8, 32, 25
ES_GENS, ES_LAM, ES_NB, ES_SIGMA = 50, 16, 8, 0.01
DEV_EVERY = 10
EVERY_E = 12
CELLS = [('T3', dict(ramp=0.5, A=0.001)), ('T1', dict(ratio=30., A=0.001))]
CELLS2 = [('T3', dict(ramp=0.25, A=0.001))]           # secondary cell (gate.py)
NTR, NDEV, NTE = 20, 20, 40
VARIANTS = dict(full=dict(teacher=True, es=True, gains_only=False, noforce=False, every=1),
                a=dict(teacher=False, es=True, gains_only=False, noforce=False, every=1),
                c=dict(teacher=True, es=True, gains_only=True, noforce=False, every=1),
                d=dict(teacher=True, es=True, gains_only=False, noforce=True, every=1),
                e=dict(teacher=True, es=True, gains_only=False, noforce=False, every=EVERY_E))


# ------------------------------------------------------------------ numba core
@njit(cache=True)
def _binom(nd, j):
    c = 1.
    for q in range(j): c = c * (nd - q) / (q + 1)
    return c


@njit(cache=True)
def li_update(i, F, y, a, sc, buf, prv):
    """LockInK.update (baselines_adaptive) in numba. sc scalars:
    0 SFr 1 SFi 2 SYr 3 SYi 4 SA 5 nbad 6 cnt 7 k 8 init 9 Hr 10 Hi 11 dgain 12 W 13 Nc"""
    W = int(sc[12]); Nc = int(sc[13])
    if sc[8] == 0.:
        for q in range(NDIFF + 1):
            prv[q, 0] = F; prv[q, 1] = y
        sc[8] = 1.
        return sc[7]
    for q in range(NDIFF):
        prv[q, 0] = prv[q + 1, 0]; prv[q, 1] = prv[q + 1, 1]
    prv[NDIFF, 0] = F; prv[NDIFF, 1] = y
    dF = 0.; dy = 0.
    for q in range(NDIFF + 1):
        c = _binom(NDIFF, q) * (1. if (NDIFF - q) % 2 == 0 else -1.)
        dF += c * prv[q, 0]; dy += c * prv[q, 1]
    th = 2 * math.pi * (i % Nc) / Nc
    rr = math.cos(th); ri = -math.sin(th)
    j = int(sc[6]) % W
    zFr = dF * rr; zFi = dF * ri; zYr = dy * rr; zYi = dy * ri
    ok = F > 0.5
    sc[0] += zFr - buf[0, j]; sc[1] += zFi - buf[1, j]; sc[2] += zYr - buf[2, j]; sc[3] += zYi - buf[3, j]
    sc[4] += a - buf[4, j]
    sc[5] += (0. if ok else 1.) - (0. if buf[5, j] > 0.5 else 1.)
    buf[0, j] = zFr; buf[1, j] = zFi; buf[2, j] = zYr; buf[3, j] = zYi; buf[4, j] = a; buf[5, j] = 1. if ok else 0.
    sc[6] += 1.
    if int(sc[6]) % (4 * W) == 0:
        for r in range(5):
            s = 0.
            for q in range(W): s += buf[r, q]
            sc[r] = s
    thr = max(1e-9 * W, RHO * 0.5 * sc[11] * sc[4])
    aY = math.sqrt(sc[2] ** 2 + sc[3] ** 2)
    if sc[6] >= W and sc[5] < 0.5 and aY > thr:
        hr = sc[9]; hi = sc[10]; den = hr * hr + hi * hi
        ZFr = (sc[0] * hr + sc[1] * hi) / den; ZFi = (sc[1] * hr - sc[0] * hi) / den
        kk = math.sqrt(ZFr ** 2 + ZFi ** 2) / aY
        sc[7] = min(max(kk, KMIN), KMAX)
    return sc[7]


@njit(cache=True)
def features(log, n, st, mu, sd, mask, out):
    """causal features at tick n (uses log rows < n only). st: 0 k_raw 1 k_hold 2 ef_int."""
    def ii(j):
        return max(n - 1 - j, 0)
    F = log[:, 5]; y = log[:, 2]
    out[0] = math.log(st[0] / 1e4); out[1] = math.log(st[1] / 1e4)
    for j in range(12):
        out[2 + j] = P.FD - F[ii(j)]
        out[14 + j] = (y[ii(j)] - y[ii(j + 1)]) / 1e-5
    out[26] = (log[ii(0), 11] - y[ii(0)]) / 1e-3
    out[27] = (y[ii(0)] - y[ii(12)]) / 1e-3
    out[28] = (y[ii(0)] - y[ii(50)]) / 1e-3
    out[29] = st[2]
    m = 0.; c = 0
    for j in range(50):
        m += F[ii(j)]; c += 1
    out[30] = P.FD - m / c
    m = 0.; v = 0.
    for j in range(25): m += F[ii(j)]
    m /= 25
    for j in range(25): v += (F[ii(j)] - m) ** 2
    out[31] = math.sqrt(v / 25) / 0.3
    out[32] = (log[ii(0), 1] - log[ii(1), 1]) / TC / 0.02
    out[33] = math.log(max(log[ii(0), 9], 1e-9))
    out[34] = math.log(max(log[ii(0), 10], 1e-9))
    for q in range(NF):
        v = (out[q] - mu[q]) / sd[q]
        v = min(max(v, -6.), 6.)
        out[q] = v * mask[q]


@njit(cache=True)
def mlp(x, W1, b1, W2, b2, W3, b3, z):
    h1 = np.empty(W1.shape[0]); h2 = np.empty(W2.shape[0])
    for a in range(W1.shape[0]):
        s = b1[a]
        for q in range(W1.shape[1]): s += W1[a, q] * x[q]
        h1[a] = math.tanh(s)
    for a in range(W2.shape[0]):
        s = b2[a]
        for q in range(W2.shape[1]): s += W2[a, q] * h1[q]
        h2[a] = math.tanh(s)
    for a in range(W3.shape[0]):
        s = b3[a]
        for q in range(W3.shape[1]): s += W3[a, q] * h2[q]
        z[a] = s


@njit(cache=True)
def prior(st, noforce, pr):
    """structured prior: the stiffness rule (EtR-A / teacher rule) on a stiffness centre, the
    causal lock-in k_hold (or the constant K_C without force, ablation d).
    pr: log M0, log D0, log KPF0, log ki0, phi (k_hat feed-forward basis, m; unused), log k centre."""
    kh = K_C if noforce else st[1]
    pr[0] = math.log(M_R); pr[1] = math.log(max(RULE_P[3], kh / RULE_P[0])); pr[2] = math.log(RULE_P[1] / kh)
    pr[3] = math.log(RULE_P[2]); pr[5] = math.log(kh)
    pr[4] = P.FD * (1.0 / kh - 1.0 / st[4]) if (st[4] > 0 and not noforce) else 0.


@njit(cache=True)
def z_to_gains(z, pr, cff, span, gains_only, g):
    kk = math.exp(pr[5] + span[1] * math.tanh(z[1]))          # inferred contact stiffness
    g[0] = min(max(math.exp(pr[0] + span[0] * math.tanh(z[0])), 0.25), 3.0)
    g[1] = min(max(max(RULE_P[3], kk / RULE_P[0]), 20.), 2500.)
    g[2] = 0.
    g[4] = min(max(RULE_P[1] / kk * math.exp(span[2] * math.tanh(z[2])), 1e-6), 1e-3)
    g[3] = min(max(math.exp(pr[3] + span[3] * math.tanh(z[3])), 0.), 30.)
    g[5] = 0. if gains_only else min(max(cff * pr[4] + g[4] * FU * math.tanh(z[4]), -0.02), 0.02)


@njit(cache=True)
def pstep(log, n, Tc, sc, buf, prv, st, mu, sd, mask, W1, b1, W2, b2, W3, b3, cff, span, every, gains_only, noforce,
          use_net, ztau, x, zr, z, pr, gcur, gout):
    """one control tick of the student: lock-in update + dither every tick, network every
    `every` ticks (prior + residual both held between decisions).
    st: 0 k_raw 1 k_hold 2 ef_int 3 a_prev 4 k_ref."""
    Nc = int(sc[13])
    if n > 0:
        F = log[n - 1, 5]
        ke = li_update(n - 1, F, log[n - 1, 2], st[3], sc, buf, prv)
        st[0] = ke
        if TREL > 0: ke = max(ke, st[1] * math.exp(-Tc / TREL))
        st[1] = ke
        al = math.exp(-Tc / 0.3)
        st[2] = al * st[2] + (1 - al) * (P.FD - F)
        if n * Tc >= TREF:
            st[4] = st[1] if st[4] <= 0 else min(st[4], st[1])
    prior(st, noforce, pr)
    a = min(max(FDITH / st[1], AMIN), AMAX)
    ud = a * math.sin(2 * math.pi * (n % Nc) / Nc)
    st[3] = a
    if n % every == 0:
        if n > 0 and use_net:
            features(log, n, st, mu, sd, mask, x)
            mlp(x, W1, b1, W2, b2, W3, b3, zr)
            az = math.exp(-every * Tc / ztau) if ztau > 0 else 0.
            for q in range(z.shape[0]): z[q] = az * z[q] + (1 - az) * zr[q]   # first-order output filter
        else:
            for q in range(z.shape[0]): z[q] = 0.
        z_to_gains(z, pr, cff[0], span, gains_only, gcur)
    for q in range(6): gout[q] = gcur[q]
    gout[5] += ud


def init_frontend(pl):
    Tc = pl['Tc']; Nc = max(4, int(round(1.0 / (Tc * F_DITHER)))); W = int(NCYC * Nc); f = 1.0 / (Nc * Tc)
    afc = 1.0 - math.exp(-2 * math.pi * pl['fc'] * 1e-3); zz = np.exp(-2j * np.pi * f * 1e-3)
    H = afc / (1.0 - (1.0 - afc) * zz)
    sc = np.zeros(14); sc[7] = K0; sc[9] = H.real; sc[10] = H.imag; sc[11] = (2 * math.sin(math.pi / Nc)) ** NDIFF
    sc[12] = W; sc[13] = Nc
    buf = np.zeros((6, W)); buf[5, :] = 1.
    prv = np.zeros((NDIFF + 1, 2)); st = np.array([K0, K0, 0., 0., 0.])
    return sc, buf, prv, st


# ------------------------------------------------------------------ policy
def init_net(seed):
    rng = np.random.default_rng(seed)
    return dict(W1=rng.standard_normal((NH, NF)) / math.sqrt(NF), b1=np.zeros(NH),
                W2=rng.standard_normal((NH, NH)) / math.sqrt(NH), b2=np.zeros(NH),
                W3=np.zeros((NO, NH)), b3=np.zeros(NO), mu=np.zeros(NF), sd=np.ones(NF), cff=np.zeros(1),
                span=np.array(SPAN, float))


class Cand2Policy:
    """causal student; never reads obs.k_true or scene parameters."""
    oracle = False

    def __init__(self, net, every=1, gains_only=False, noforce=False, use_net=True, record=False, ztau=None):
        self.ztau = ZTAU if ztau is None else ztau
        self.net = {k: np.ascontiguousarray(v, dtype=float) for k, v in net.items()}
        self.every = int(every); self.go = bool(gains_only); self.use_net = bool(use_net)
        self.mask = np.ones(NF); self.nf = bool(noforce); self.mask[OFF_FEATS] = 0.
        if noforce: self.mask[FORCE_FEATS] = 0.
        self.record = record

    def reset(self, pl):
        self.Tc = pl['Tc']; self.sc, self.buf, self.prv, self.st = init_frontend(pl)
        self.x = np.zeros(NF); self.z = np.zeros(NO); self.gcur = np.zeros(6); self.g = np.zeros(6); self.pr = np.zeros(6); self.zr = np.zeros(NO)
        self.raw = np.zeros(NF); self.khist = []

    def __call__(self, o):
        nt = self.net
        pstep(o.log, o.n, self.Tc, self.sc, self.buf, self.prv, self.st, nt['mu'], nt['sd'], self.mask,
              nt['W1'], nt['b1'], nt['W2'], nt['b2'], nt['W3'], nt['b3'], nt['cff'], nt['span'], self.every, self.go, self.nf,
              self.use_net, self.ztau, self.x, self.zr, self.z, self.pr, self.gcur, self.g)
        if self.record: self.khist.append(self.st[1])
        return self.g

    def raw_features(self, o):
        """unnormalised features (for normalisation statistics / DAgger records)."""
        features(o.log, o.n, self.st, np.zeros(NF), np.ones(NF), np.ones(NF), self.raw)
        return self.raw.copy()


# ------------------------------------------------------------------ teacher
def teacher_from(p, rule):
    """rule = (wI, cP, ki, Dmin, preview[, latency]); preview < 0: no feed-forward (plain TruthK)"""
    from baselines_ref import TruthKFF
    wI, cP, ki, Dmin, prev = rule[:5]; lat = rule[5] if len(rule) > 5 else 0.0
    if prev < 0:
        from baselines_ref import TruthK
        return TruthK(lat, wI, cP, ki, Dmin)
    return TruthKFF(p, prev, lat, wI, cP, ki, Dmin)


def gains_to_t(g, kT, pr, span):
    """teacher -> student targets in tanh space, relative to the student's prior at the same state.
    The teacher is the truth-stiffness rule, a known function of k_true; its M, ki equal the
    prior's and its KPF is cP/k, so imitating its unsaturated action is regressing log k_true:
    target z1 = (log k_T - log k_centre)/span_k; z0, z2, z3 targets 0 (teacher = prior form);
    u is not distilled. Last two entries: teacher u_ff and phi (diagnostics)."""
    t = np.array([0., (math.log(kT) - pr[5]) / span[1], 0., 0., 0.])
    return np.concatenate([np.clip(t, -0.97, 0.97), [g[5], pr[4]]])


class DaggerRollout:
    """executes teacher (+ student dither) or student; records student raw features and teacher
    labels every `stride` ticks."""
    def __init__(self, student, teacher, use_teacher, stride=3, offset=0):
        self.s = student; self.t = teacher; self.ut = use_teacher; self.stride = stride; self.off = offset

    def reset(self, pl):
        self.s.reset(pl); self.t.reset(pl); self.X = []; self.Y = []

    def __call__(self, o):
        gs = np.array(self.s(o))
        gt = np.clip(P._pad(self.t(o)), P.GMIN, P.GMAX)
        if o.n >= 2 and (o.n + self.off) % self.stride == 0:
            self.X.append(self.s.raw_features(o)); self.Y.append(gains_to_t(gt, self.t.k, self.s.pr, self.s.net['span']))
        if self.ut:
            d = self.s.g[5] - (self.s.gcur[5])        # the student's dither this tick
            g = gt.copy(); g[5] += d; return g
        return gs


# ------------------------------------------------------------------ scenes
def batch_scenes(seed, n_each):
    """train split scenes: half at the cell doses, half randomised."""
    assert S.SPLIT['train'] <= seed < S.SPLIT['dev']
    rng = np.random.default_rng(seed); out = []
    for task in ['T3'] * n_each + ['T1'] * n_each:
        if rng.random() < 0.5:
            dose = dict(CELLS[0][1]) if task == 'T3' else dict(CELLS[1][1])
        else:
            dose = dict(A=rng.uniform(0.0005, 0.005), ratio=rng.uniform(10, 30))
            if task == 'T3': dose['ramp'] = rng.uniform(0.2, 0.6)
        out.append((task,) + S.make(task, rng, dose))
    return out


def gate_cells():
    g = json.load(open('gate_PC125.json'))
    return {(c['task'], json.dumps(c['dose'], sort_keys=True)): c for c in g}


def b1_gains(task, dose=None):
    gc = gate_cells()
    dose = dose or (CELLS[0][1] if task == 'T3' else CELLS[1][1])
    return tuple(gc[(task, json.dumps(dose, sort_keys=True))]['B1']['gains'])


# ------------------------------------------------------------------ workers
def _ctrl(spec, p):
    kind = spec[0]
    if kind == 'gains': return spec[1]
    if kind == 'net':
        _, net, cfg = spec
        return Cand2Policy(net, every=cfg['every'], gains_only=cfg['gains_only'], noforce=cfg['noforce'])
    if kind == 'teacher': return teacher_from(p, spec[1])
    if kind == 'truthk':
        from baselines_ref import TruthK
        return TruthK(*spec[1])
    if kind == 'etr':
        from baselines_adaptive import EtR
        rp, ep = spec[1]
        return EtR(rp[0], rp[1], rp[2], rp[3], rule='truthk', **ep)
    raise ValueError(kind)


def _run(job):
    spec, p, T, nseed = job
    c0 = time.process_time()
    pol = _ctrl(spec, p)
    log = P.episode(p, pol, T=T, platform='PC125', seed=nseed)
    m = S.metrics(log, Tc=TC)
    if log is not None and spec[0] in ('net', 'teacher'):
        w = log[:, 0] >= 1.0; D = log[w, 9]; K = log[w, 10]
        m['D_med'] = float(np.median(D)); m['D_p5'] = float(np.percentile(D, 5)); m['D_p95'] = float(np.percentile(D, 95))
        m['KPF_med'] = float(np.median(K)); m['KPF_p5'] = float(np.percentile(K, 5)); m['KPF_p95'] = float(np.percentile(K, 95))
    m['cpu'] = time.process_time() - c0
    return m


def _run_dagger(job):
    net, cfg, rule, use_t, p, T, nseed, off = job
    st = Cand2Policy(net, every=cfg['every'], gains_only=cfg['gains_only'], noforce=cfg['noforce'],
                     use_net=net is not None and not cfg.get('prior', False))
    roll = DaggerRollout(st, teacher_from(p, rule) if rule is not None else _Null(), use_t, offset=off)
    log = P.episode(p, roll, T=T, platform='PC125', seed=nseed)
    m = S.metrics(log, Tc=TC)
    return np.array(roll.X), np.array(roll.Y), m['J']


class _Null:
    k = K0
    def reset(self, pl): pass
    def __call__(self, o): return P0


_POOL = None
def pool():
    global _POOL
    if _POOL is None:
        from multiprocessing import Pool
        _POOL = Pool(2)
    return _POOL


def run_many(jobs, chunk=None):
    return pool().map(_run, jobs, chunksize=chunk or max(1, len(jobs) // 16))


# ------------------------------------------------------------------ cohorts and scoring
def cohorts(split, n):
    return [(t, S.cohort(t, n, S.SPLIT[split], d)) for t, d in CELLS]


def eval_cohorts(spec, coh, split):
    jobs = [(spec, p, T, S.SPLIT[split] + i) for _, sc in coh for i, (p, T) in enumerate(sc)]
    r = run_many(jobs); out = {}; k = 0
    for t, sc in coh:
        rows = r[k:k + len(sc)]; k += len(sc)
        pk = np.array([m['peak'] if np.isfinite(m['peak']) else 1e9 for m in rows])
        out[t] = dict(J=float(np.mean([m['J'] for m in rows])), rmse=float(np.nanmean([m['rmse'] for m in rows])),
                      feas=bool(np.mean(pk <= P.FD + 3.0) >= 0.95), rows=rows)
    return out


_REF = {}
def refs():
    """B1 J on dev and train cohorts (normalisers for the selection score)."""
    if not _REF:
        for split, n in [('dev', NDEV), ('train', NTR)]:
            coh = cohorts(split, n)
            for t, sc in coh:
                rr = run_many([(('gains', b1_gains(t)), p, T, S.SPLIT[split] + i) for i, (p, T) in enumerate(sc)])
                _REF[(split, t)] = float(np.mean([m['J'] for m in rr]))
    return _REF


def select_score(spec, dev, train):
    R = refs()
    d = eval_cohorts(spec, dev, 'dev'); tr = eval_cohorts(spec, train, 'train')
    score = float(np.mean([d[t]['J'] / R[('dev', t)] for t in d]))
    feas = all(tr[t]['feas'] for t in tr)
    strip = lambda e: {t: {k: v for k, v in e[t].items() if k != 'rows'} for t in e}
    return score, feas, strip(d), strip(tr)


# ------------------------------------------------------------------ teacher tuning
TEACH_GRID = list(itertools.product([40., 80., 160.], [1., 2., 4.], [10., 30.], [20., 100.], [0.0, 0.05]))

def tune_teacher():
    t0 = time.time(); tr = cohorts('train', NTR); R = refs(); res = []
    for rule in TEACH_GRID:
        e = eval_cohorts(('teacher', rule), tr, 'train')
        sc = float(np.mean([e[t]['J'] / R[('train', t)] for t in e])); feas = all(e[t]['feas'] for t in e)
        res.append((sc + (0 if feas else 1e3), rule, {t: e[t]['J'] for t in e}))
        print('teacher', rule, round(sc, 3), feas, {t: round(e[t]['J'], 3) for t in e}, flush=True)
    res.sort(key=lambda r: r[0]); best = res[0]
    np.savez('candidate2_teacher.npz', rule=np.array(best[1]), info=json.dumps(dict(
        best_rule=best[1], best_score=best[0], best_train_J=best[2], n_configs=len(TEACH_GRID),
        top5=[(r[0], r[1], r[2]) for r in res[:5]], seconds=time.time() - t0, B1_train_J={str(k): v for k, v in R.items()})))
    print('best teacher', best, flush=True)
    return best[1]


def load_teacher():
    return tuple(float(x) for x in np.load('candidate2_teacher.npz')['rule'])


# ------------------------------------------------------------------ supervised training (numpy)
def fwd(net, X):
    a1 = np.tanh(X @ net['W1'].T + net['b1']); a2 = np.tanh(a1 @ net['W2'].T + net['b2'])
    z = a2 @ net['W3'].T + net['b3']; return a1, a2, z


def fit(net, X, Y, wout, epochs, rng, lr=1e-3, bs=512, tau=0.5):
    """Adam on the weighted expectile loss in tanh space: weight 2*tau where the output is below
    the target, 2*(1-tau) above (tau = 0.5: MSE; tau > 0.5 penalises under-estimates)."""
    keys = ['W1', 'b1', 'W2', 'b2', 'W3', 'b3']
    m = {k: np.zeros_like(net[k]) for k in keys}; v = {k: np.zeros_like(net[k]) for k in keys}; t = 0
    N = len(X); loss = 0.
    for ep in range(epochs):
        perm = rng.permutation(N); tot = 0.
        for s in range(0, N, bs):
            idx = perm[s:s + bs]; x = X[idx]; y = Y[idx]
            a1, a2, z = fwd(net, x); o = np.tanh(z); e = (o - y) * wout * np.where(o < y, 2 * tau, 2 * (1 - tau))
            tot += float(np.sum(e * (o - y)))
            dz = 2 * e * (1 - o * o) / len(idx)
            g = {'W3': dz.T @ a2, 'b3': dz.sum(0)}
            d2 = (dz @ net['W3']) * (1 - a2 * a2); g['W2'] = d2.T @ a1; g['b2'] = d2.sum(0)
            d1 = (d2 @ net['W2']) * (1 - a1 * a1); g['W1'] = d1.T @ x; g['b1'] = d1.sum(0)
            t += 1
            for k in keys:
                m[k] = 0.9 * m[k] + 0.1 * g[k]; v[k] = 0.999 * v[k] + 0.001 * g[k] ** 2
                net[k] -= lr * (m[k] / (1 - 0.9 ** t)) / (np.sqrt(v[k] / (1 - 0.999 ** t)) + 1e-8)
        loss = tot / N
    return loss


# ------------------------------------------------------------------ training
def cfg_of(variant):
    v = VARIANTS[variant]; return dict(every=v['every'], gains_only=v['gains_only'], noforce=v['noforce'])


def norm_stats(net, seed):
    """feature mean/sd from rollouts of the prior (zero-output) policy on train scenes."""
    sc = batch_scenes(S.SPLIT['train'] + 900 + seed, 4)
    jobs = [(net, dict(every=1, gains_only=False, noforce=False, prior=True), None, False, p, T, 10 ** 6 + 5000 + i, 0)
            for i, (_, p, T) in enumerate(sc)]
    r = pool().map(_run_dagger, jobs, chunksize=1)
    X = np.concatenate([a[0] for a in r]); mu = X.mean(0); sd = X.std(0) + 1e-6
    # angle-free robustness: scale-only features keep a floor
    return mu, np.maximum(sd, 1e-3 * (np.abs(mu) + 1e-3))


def netspec(net, cfg): return ('net', {k: v for k, v in net.items()}, cfg)


def train(variant, seed, log=print):
    V = VARIANTS[variant]; cfg = cfg_of(variant); rng = np.random.default_rng(100 + seed)
    t0 = time.time(); cpu0 = time.process_time(); neps = 0
    net = init_net(seed); net['mu'], net['sd'] = norm_stats(net, seed); neps += 8
    dev = cohorts('dev', NDEV); trc = cohorts('train', NTR); refs()
    hist = []; best = None
    wout = np.zeros(NO); wout[1] = 1.        # only the stiffness head is distilled (see gains_to_t)
    m = np.ones(NF); m[OFF_FEATS] = 0.
    if cfg['noforce']: m[FORCE_FEATS] = 0.
    # ---- stage 1: DAgger
    if V['teacher']:
        rule = TEACHER_RULE; XA = []; YA = []
        for it in range(N_DAG):
            beta = 0.5 ** it
            sc = batch_scenes(S.SPLIT['train'] + 300 + 20 * seed + it, EP_DAG // 2)
            ut = rng.random(len(sc)) < beta
            jobs = [({k: v for k, v in net.items()}, cfg, rule, bool(ut[i]), p, T, 10 ** 6 + 7000 + 100 * it + i, i % 3)
                    for i, (_, p, T) in enumerate(sc)]
            r = pool().map(_run_dagger, jobs, chunksize=2); neps += len(jobs)
            for X, Y, _ in r:
                if len(X): XA.append(X); YA.append(Y)
            Xr = np.concatenate(XA); YU = np.concatenate(YA); Y = YU[:, :NO]; uT = YU[:, NO]; ph = YU[:, NO + 1]
            Xn = np.clip((Xr - net['mu']) / net['sd'], -6, 6)
            Xn = Xn * m
            loss = fit(net, Xn, Y, wout, EPOCHS if it < N_DAG - 1 else 2 * EPOCHS, rng, tau=EXPECTILE)
            # stiffness feed-forward coefficient: least squares of the teacher's u_ff on phi
            if USE_FF and not (cfg['gains_only'] or cfg['noforce']) and np.sum(ph * ph) > 0:
                net['cff'] = np.array([float(np.sum(ph * uT) / np.sum(ph * ph))])
            score, feas, d, trr = select_score(netspec(net, cfg), dev, trc); neps += 2 * (NDEV + NTR)
            rec = dict(stage=1, it=it + 1, beta=beta, n_samples=int(len(Y)), loss=loss,
                       student_J_rollouts=float(np.mean([q[2] for q, u in zip(r, ut) if not u])) if (~ut).any() else None,
                       dev=d, train_feas=feas, dev_score=score, wall=time.time() - t0)
            hist.append(rec)
            log(f'[{variant} s{seed}] DAgger {it+1} beta {beta:.2f} N {len(Y)} loss {loss:.4f} dev T3 {d["T3"]["J"]:.3f} '
                f'T1 {d["T1"]["J"]:.3f} score {score:.3f} feas {feas} cff {net["cff"][0]:.3f} wall {time.time()-t0:.0f}s')
            if feas and (best is None or score < best[0]):
                best = (score, {k: v.copy() for k, v in net.items()}, ('dagger', it + 1), d)
        if best is None:
            best = (1e9, {k: v.copy() for k, v in net.items()}, ('dagger', N_DAG), None)
        stage1 = dict(best[1]); stage1_sel = dict(score=best[0], at=best[2], dev=best[3])
        net = {k: v.copy() for k, v in best[1].items()}
    else:
        stage1 = None; stage1_sel = None
        score, feas, d, trr = select_score(netspec(net, cfg), dev, trc); neps += 2 * (NDEV + NTR)
        best = (score if feas else 1e9, {k: v.copy() for k, v in net.items()}, ('init', 0), d)
        log(f'[{variant} s{seed}] init dev score {score:.3f} feas {feas}')
    # ---- stage 2: ES on the last layer
    from baselines_learned import SepCMA
    best2 = (best[0], {k: v.copy() for k, v in net.items()}, ('es', 0) if not V['teacher'] else best[2], best[3])
    if V['es']:
        th0 = np.concatenate([net['W3'].ravel(), net['b3'], net['cff']])
        es = SepCMA(th0, ES_SIGMA, ES_LAM, 1000 + seed)
        def unpack(th):
            q = {k: v for k, v in net.items()}; q = dict(q)
            q['W3'] = th[:NO * NH].reshape(NO, NH).copy(); q['b3'] = th[NO * NH:NO * NH + NO].copy()
            q['cff'] = np.zeros(1) if (cfg['gains_only'] or cfg['noforce'] or not USE_FF) else th[-1:].copy(); return q
        for g in range(ES_GENS):
            X = es.ask(); assert g < 100; sc = batch_scenes(S.SPLIT['train'] + 100 * seed + g, ES_NB // 2)
            jobs = [(netspec(unpack(X[k]), cfg), p, T, 10 ** 6 + 100 * g + i) for k in range(ES_LAM) for i, (_, p, T) in enumerate(sc)]
            jobs += [(('gains', b1_gains(t)), p, T, 10 ** 6 + 100 * g + i) for i, (t, p, T) in enumerate(sc)]
            r = run_many(jobs, chunk=4); neps += len(jobs)
            Jm = np.array([q['J'] for q in r[:ES_LAM * len(sc)]]).reshape(ES_LAM, len(sc))
            Jref = np.array([q['J'] for q in r[ES_LAM * len(sc):]])
            f = (Jm / (Jref + 0.01)).mean(1); es.tell(X, f)
            rec = dict(stage=2, gen=g + 1, f_best=float(f.min()), f_mean=float(f.mean()), sigma=float(es.sigma))
            if (g + 1) % DEV_EVERY == 0 or g + 1 == ES_GENS:
                cand = unpack(es.m)
                score, feas, d, trr = select_score(netspec(cand, cfg), dev, trc); neps += 2 * (NDEV + NTR)
                rec.update(dev=d, train_feas=feas, dev_score=score)
                if feas and score < best2[0]: best2 = (score, cand, ('es', g + 1), d)
                log(f'[{variant} s{seed}] ES {g+1} f {f.mean():.3f} sigma {es.sigma:.4f} dev T3 {d["T3"]["J"]:.3f} T1 {d["T1"]["J"]:.3f} '
                    f'score {score:.3f} feas {feas} best@{best2[2]} wall {time.time()-t0:.0f}s')
            hist.append(rec)
    out = dict(variant=variant, seed=seed, cfg=cfg, selected_at=best2[2], dev_score=best2[0], dev=best2[3],
               stage1_selected=stage1_sel, episodes=neps, wall_s=time.time() - t0,
               budgets=dict(N_DAG=N_DAG if V['teacher'] else 0, EP_DAG=EP_DAG, EPOCHS=EPOCHS, ES_GENS=ES_GENS if V['es'] else 0,
                            ES_LAM=ES_LAM, ES_NB=ES_NB, ES_SIGMA=ES_SIGMA, DEV_EVERY=DEV_EVERY), hist=hist)
    sv = {k: v for k, v in best2[1].items()}
    if stage1 is not None: sv.update({'s1_' + k: v for k, v in stage1.items()})
    np.savez(f'candidate2_{variant}_seed{seed}.npz', **sv, meta=json.dumps(out, default=float))
    return out


def load_net(variant, seed, stage1=False):
    z = np.load(f'candidate2_{variant}_seed{seed}.npz')
    pre = 's1_' if stage1 else ''
    return {k: z[pre + k] for k in ['W1', 'b1', 'W2', 'b2', 'W3', 'b3', 'mu', 'sd', 'cff', 'span']}, json.loads(str(z['meta']))


# ------------------------------------------------------------------ checks
def check_equivalence():
    """zero-output student (numba lock-in + dither + rule prior, c_ff = 0) vs baselines_adaptive.EtR-A
    (same estimator settings and rule) on one train scene per cell."""
    from baselines_adaptive import EtR
    out = []
    for task, dose in CELLS:
        p, T = S.cohort(task, 1, S.SPLIT['train'], dose)[0]
        ep = dict(est='lockin', ncyc=NCYC, f_dither=F_DITHER, ndiff=NDIFF, rho=RHO, trel=TREL, Fdith=FDITH)
        ref = EtR(*RULE_P, rule='truthk', **ep)       # zero network output must reproduce EtR-A
        lr = P.episode(p, ref, T=T, seed=1)
        net = init_net(0); pol = Cand2Policy(net, use_net=False, record=True)
        lc = P.episode(p, pol, T=T, seed=1)
        kr = np.array(ref.khist[:len(lr)]); kc = np.array(pol.khist[:len(lc)])
        out.append(dict(task=task, max_rel_k_diff=float(np.max(np.abs(kr - kc) / kr)),
                        max_force_diff=float(np.max(np.abs(lr[:, 3] - lc[:, 3]))),
                        max_ycmd_diff=float(np.max(np.abs(lr[:, 11] - lc[:, 11])))))
    return out


def inference_time(net, cfg, n_rep=3):
    """wall time per decision of the full policy call (front end + features + MLP) and of the
    MLP alone, measured inside an episode replay on one T3 scene."""
    p, T = S.cohort('T3', 1, S.SPLIT['test'], CELLS[0][1])[0]
    pol = Cand2Policy(net, **{k: cfg[k] for k in ['every', 'gains_only', 'noforce']})
    log = P.episode(p, pol, T=T, seed=1)
    pl = P.PLATFORMS['PC125']; ts = []
    for _ in range(n_rep):
        pol.reset(pl); t0 = time.perf_counter()
        for n in range(len(log)): pol(P.Obs(log, n))
        ts.append((time.perf_counter() - t0) / len(log))
    x = np.random.default_rng(0).standard_normal(NF); z = np.zeros(NO); nt = pol.net
    t0 = time.perf_counter()
    for _ in range(20000): mlp(x, nt['W1'], nt['b1'], nt['W2'], nt['b2'], nt['W3'], nt['b3'], z)
    tm = (time.perf_counter() - t0) / 20000
    t0 = time.perf_counter()
    for _ in range(2000):
        a1 = np.tanh(nt['W1'] @ x + nt['b1']); a2 = np.tanh(nt['W2'] @ a1 + nt['b2']); zz = nt['W3'] @ a2 + nt['b3']
    tn = (time.perf_counter() - t0) / 2000
    return dict(policy_call_us_per_tick=float(np.min(ts) * 1e6), mlp_numba_us=tm * 1e6, mlp_numpy_us=tn * 1e6,
                note='single core, includes Python call overhead and Obs construction; control period 8000 us')


# ------------------------------------------------------------------ test
SUMK = ['J', 'rmse', 'p95', 'maxerr', 'inband', 'peak', 't_over', 'loss', 'loss_max', 'osc', 'shear', 'diverged', 'frac_peak_over']

def test(seeds=(0, 1, 2), variants=('full', 'a', 'b', 'c', 'd', 'e'), cells=None):
    import tune as U
    gc = gate_cells(); ad = json.load(open('adaptive_results.json'))
    etr = {c['task']: (c['methods']['EtR-A']['tuning']['rule_params'], c['methods']['EtR-A']['tuning']['est_params']) for c in ad['cells']}
    res = json.load(open('candidate2_results.json')) if os.path.exists('candidate2_results.json') else {}
    res.setdefault('cells', {}); res.setdefault('training', {})
    for v in ('full', 'a', 'c', 'd', 'e'):
        for s_ in seeds:
            fn = f'candidate2_{v}_seed{s_}.npz'
            if os.path.exists(fn):
                meta = json.loads(str(np.load(fn)['meta']))
                res['training'][f'{v}_seed{s_}'] = {k: meta[k] for k in ['selected_at', 'dev_score', 'dev', 'stage1_selected',
                                                                          'episodes', 'wall_s', 'budgets']}
    for task, dose in (cells or CELLS):
        key = f'{task} {json.dumps(dose, sort_keys=True)}'; g = gc[(task, json.dumps(dose, sort_keys=True))]
        sc = S.cohort(task, NTE, S.SPLIT['test'], dose)
        cell = res['cells'].get(key, dict(task=task, dose=dose, methods={}))
        refsp = {'B1': ('gains', tuple(g['B1']['gains'])), 'B3b_0.1': ('truthk', (0.1,) + tuple(g['B3b']['0.1']['rule'])),
                 'EtR-A': ('etr', etr[task] if dose.get('ramp', 0.5) == 0.5 else etr[task])}
        rr = {}
        for name, spec in refsp.items():
            rows = run_many([(spec, p, T, S.SPLIT['test'] + i) for i, (p, T) in enumerate(sc)])
            rr[name] = rows; cell['methods'][name] = dict(test=U.summary(rows))
        jr = {k: [m['J'] for m in v] for k, v in rr.items()}
        for var in variants:
            per = {}; allrows = []
            for s in seeds:
                vv = 'full' if var == 'b' else var
                fn = f'candidate2_{vv}_seed{s}.npz'
                if not os.path.exists(fn): continue
                net, meta = load_net(vv, s, stage1=(var == 'b'))
                rows = run_many([(netspec(net, meta['cfg']), p, T, S.SPLIT['test'] + i) for i, (p, T) in enumerate(sc)])
                jj = [m['J'] for m in rows]; allrows.append(rows)
                per[str(s)] = dict(test=U.summary(rows),
                                   dJ_vs_B1=U.paired_ci(jj, jr['B1']), dJ_vs_EtRA=U.paired_ci(jj, jr['EtR-A']),
                                   dJ_vs_B3b01=U.paired_ci(jj, jr['B3b_0.1']),
                                   gains={k: float(np.mean([m[k] for m in rows if k in m])) for k in
                                          ['D_med', 'D_p5', 'D_p95', 'KPF_med', 'KPF_p5', 'KPF_p95']})
            if not per: continue
            # seed-averaged per-scene J (for a pooled paired CI over scenes)
            jm = np.mean([[m['J'] for m in rows] for rows in allrows], 0)
            e = dict(seeds=per,
                     mean={k: float(np.mean([per[s]['test'][k] for s in per])) for k in SUMK},
                     sd={k: float(np.std([per[s]['test'][k] for s in per])) for k in SUMK},
                     seedavg_dJ_vs_B1=U.paired_ci(jm, jr['B1']), seedavg_dJ_vs_EtRA=U.paired_ci(jm, jr['EtR-A']),
                     seedavg_dJ_vs_B3b01=U.paired_ci(jm, jr['B3b_0.1']))
            cell['methods']['cand2_' + var] = e
            print(key, var, {k: round(e['mean'][k], 3) for k in ['J', 'rmse', 'frac_peak_over', 'inband']},
                  'sd J', round(e['sd']['J'], 3), 'dJ vs B1', [round(x, 3) for x in e['seedavg_dJ_vs_B1']],
                  'vs EtR-A', [round(x, 3) for x in e['seedavg_dJ_vs_EtRA']], flush=True)
        for name in refsp:
            s = cell['methods'][name]['test']; print(key, name, {k: round(s[k], 3) for k in ['J', 'rmse', 'frac_peak_over', 'inband']})
        res['cells'][key] = cell
        json.dump(res, open('candidate2_results.json', 'w'), indent=1, default=float)
    return res


def known_behaviour(seeds=(0, 1, 2), n=4):
    """homogeneous N1 tissue, soft 15 kPa and stiff 1 MPa, breathing 1 mm: chosen D/KPF and J
    vs B1 (both cells' B1 gains)."""
    import tune as U
    out = {}
    for E in [1.5e4, 1e6]:
        sc = S.cohort('N1', n, 5000, dict(E=E, A=0.001)); row = {}
        for name, spec in [('B1_T3', ('gains', b1_gains('T3'))), ('B1_T1', ('gains', b1_gains('T1')))]:
            rows = run_many([(spec, p, T, 5000 + i) for i, (p, T) in enumerate(sc)])
            row[name] = U.summary(rows)
        for s in seeds:
            if not os.path.exists(f'candidate2_full_seed{s}.npz'): continue
            net, meta = load_net('full', s)
            rows = run_many([(netspec(net, meta['cfg']), p, T, 5000 + i) for i, (p, T) in enumerate(sc)])
            e = U.summary(rows); e.update({k: float(np.mean([m[k] for m in rows if k in m])) for k in ['D_med', 'D_p5', 'D_p95', 'KPF_med', 'KPF_p5', 'KPF_p95']})
            row[f'cand2_full_seed{s}'] = e
        out['N1_E%.1e' % E] = row
        print('N1 E', E, {k: (round(v['J'], 3), round(v.get('D_med', 0), 1), v.get('KPF_med', 0)) for k, v in row.items()}, flush=True)
    return out


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'test'
    if cmd == 'check':
        print(json.dumps(check_equivalence(), indent=1))
    elif cmd == 'teacher':
        tune_teacher()
    elif cmd == 'train':
        variant = sys.argv[2]; seed = int(sys.argv[3])
        o = train(variant, seed, log=lambda s: print(s, flush=True))
        print('done', variant, seed, o['selected_at'], o['dev_score'], o['episodes'], round(o['wall_s']), flush=True)
    elif cmd == 'test':
        test()
    elif cmd == 'test2':
        test(cells=CELLS2, variants=('full', 'a', 'b'))
    elif cmd == 'kb':
        r = json.load(open('candidate2_results.json')); r['known_behaviour_N1'] = known_behaviour()
        net, meta = load_net('full', 0); r['inference_time'] = inference_time(net, meta['cfg'])
        net, meta = load_net('e', 0) if os.path.exists('candidate2_e_seed0.npz') else (net, meta)
        r['inference_time_check_equivalence'] = check_equivalence()
        json.dump(r, open('candidate2_results.json', 'w'), indent=1, default=float)
        print(json.dumps(r['inference_time'], indent=1))
