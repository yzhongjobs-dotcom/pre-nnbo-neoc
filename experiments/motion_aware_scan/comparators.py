"""Literature-based normal force comparators; SI units, 8 ms direct commands.

Only WijType reads privileged k_eff. All other observations are columns
0 (time), 2 (encoder), 5 (filtered force), 11 (actually sent command).
YuenType additionally reads encoder x (column 1) and a delayed privileged
surface-position sensor; its class documents this information advantage.
These are task adaptations, not reproductions of the papers' full robots.
"""
import math
import numpy as np
from scipy.linalg import expm
from scipy.signal import iirpeak

try:
    from . import common
except ImportError:
    import common


def _positive(**values):
    for key, value in values.items():
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f'{key} must be finite and positive')


class _MatchedDisplacement:
    """Causal 1 kHz displacement LPF in the force measurement domain.

    As in common._build, reconstruct interval j from encoder y[j] and the
    delayed sent command log[j-ndel,11], using first-order servo interpolation.
    Carry the LPF state across ticks instead of inverting the force filter.
    Only columns 2 and 11 are read; window_fit is not used. Initial delay
    slots contain the initial encoder position, as in plant.init_state.
    """
    def __init__(self, policy):
        self.ns, self.eb = policy.ns, policy.eb
        self.afc, self.ndel = policy.afc, policy.ndel
        self.last = -1

    def update(self, o):
        if self.last < 0:
            self.initial = self.filtered = float(o.log[0, 2])
            self.last = 0
        for j in range(self.last, o.n-1):
            y = float(o.log[j, 2])
            target = float(o.log[j-self.ndel, 11]) if j >= self.ndel else self.initial
            for _ in range(self.ns):
                # same order as the plant: sample with the current position, then advance the servo
                self.filtered += self.afc*(y-self.filtered)
                y = target+(y-target)*self.eb
        self.last = o.n-1
        return self.filtered


def nominal_mpc(T, Np, q_over_r):
    """S03 (24)-(27), R=1: return first-move gain and poles of A-BK.

    The poles concern only the ideally compensated double integrator, with
    exact state feedback. They do not include the ESO, servo, LPF or delay.
    """
    _positive(T=T, q_over_r=q_over_r)
    if isinstance(Np, bool) or int(Np) != Np or Np not in (8, 15, 25, 40):
        raise ValueError('Np must be one of {8, 15, 25, 40}')
    Np = int(Np)
    A, B = np.array([[1., T], [0., 1.]]), np.array([T*T/2, T])
    H = np.column_stack((np.ones(Np), T*np.arange(1, Np+1)))
    i, j = np.indices((Np, Np))
    G = np.where(i >= j, T*T*(i-j+0.5), 0.)
    gain = np.linalg.solve(q_over_r*G.T@G+np.eye(Np), q_over_r*G.T@H)[0]
    poles = np.linalg.eigvals(A-np.outer(B, gain))
    return gain, poles


class _RossiEstimator:
    """Rossi Eq. (5)-(8), III-A.2: second-order BP and scalar EWRLS.

    iirpeak implements a prewarped bilinear bandpass of unity peak gain.
    Input y must already have the same 30 Hz LPF as measured F (see
    _MatchedDisplacement). Both channels then use identical causal BP filters.
    Subtract first observations
    before filtering to avoid initial static-offset transients. Normalize
    theta=k/K_MAX, phi=K_MAX*y_bp for numerical conditioning. No inverse
    sensor filter or motion feed-forward is added. Servo interpolation is
    measurement preprocessing, not compensation in the control command.
    Readiness uses a warm-up and discounted regressor energy, not a claim
    of statistical convergence. Before readiness the exposed k stays K_MAX.
    """
    def __init__(self, Tc, f_exc, Q, lam, P0, ready_cycles, info_min):
        _positive(Tc=Tc, f_exc=f_exc, Q=Q, P0=P0,
                  ready_cycles=ready_cycles, info_min=info_min)
        if not 0 < lam <= 1 or f_exc >= .5/Tc:
            raise ValueError('Require 0 < lam <= 1 and f_exc below Nyquist')
        self.b, self.a = iirpeak(f_exc, Q, fs=1/Tc)
        self.state = np.zeros((2, 2))
        self.offset = None
        self.theta, self.P, self.info = 1., P0, 0.
        self.lam, self.info_min = lam, info_min
        self.warm_ticks = math.ceil(ready_cycles/(f_exc*Tc))
        self.n, self.k_hat = 0, common.K_MAX

    def update(self, F, y):
        sample = np.array([F, y])
        if self.offset is None:
            self.offset = sample.copy()
        x = sample-self.offset
        bp = self.b[0]*x+self.state[:, 0]
        self.state[:, 0] = self.b[1]*x-self.a[1]*bp+self.state[:, 1]
        self.state[:, 1] = self.b[2]*x-self.a[2]*bp
        phi = common.K_MAX*bp[1]
        gain = self.P*phi/(self.lam+phi*self.P*phi)
        self.theta += gain*(bp[0]-phi*self.theta)
        self.theta = float(np.clip(self.theta, common.K_MIN/common.K_MAX, 1.))
        self.P = (self.P-gain*phi*self.P)/self.lam
        self.info = self.lam*self.info+phi*phi
        self.n += 1
        self.k_hat = (self.theta*common.K_MAX if self.n >= self.warm_ticks
                      and self.info >= self.info_min else common.K_MAX)
        return self.k_hat


class RossiType(common.DirectPolicy):
    """Rossi ICRA2016 Eq. (5)-(10), III-A.2 estimator; scalar Eq. (11).

    y points inward, so F_bp=+k*y_bp (paper z has the opposite sign).
    Delta y=(K_I/k_hat)*(FD-F)*Tc plus successive excitation differences.
    K_I=30 s^-1 is this task's ideal normalized integral-loop pole scale;
    the original experiment used 5 s^-1. Neither is the actual bandwidth
    of this delayed servo. Task adaptation: use fixed displacement amplitude
    A_E from Eq. (5) and tunable f_exc in [2,30] Hz, rather than the paper's
    experimental 0.5 Hz, 3 N force-reference sinusoid (the theory uses
    high-frequency displacement excitation). Preprocess encoder displacement
    with the servo-interpolated 1 ms, 30 Hz LPF used by this task's force
    chain (afc=1-exp(-2*pi*30*1e-3) on PC125), then identical BP and EWRLS.
    Actual ripple depends on servo lag and stiffness; amplitude is not scaled
    by k_hat. Robot joint/link compliance and motion FF are omitted.
    P0 is the normalized RLS covariance; info_min is in normalized phi^2
    units (N^2); ready_cycles is the operational warm-up in dither periods.
    """
    name = 'RossiType'
    params = dict(K_I=30., lam=.99, Q=3., f_exc=20.8,
                  P0=100., ready_cycles=3., info_min=.01, A_E=2e-5)

    def __init__(self, *, K_I=30., lam=.99, Q=3., f_exc=20.8,
                 P0=100., ready_cycles=3., info_min=.01, A_E=2e-5):
        self.params = dict(K_I=K_I, lam=lam, Q=Q, f_exc=f_exc,
                           P0=P0, ready_cycles=ready_cycles, info_min=info_min, A_E=A_E)
        _positive(K_I=K_I, A_E=A_E)
        if not 2 <= f_exc <= 30:
            raise ValueError('f_exc must be in [2, 30] Hz')
        self.__dict__.update(self.params)

    def reset(self, platform):
        super().reset(platform)
        self.displacement = _MatchedDisplacement(self)
        self.estimator = _RossiEstimator(self.Tc, self.f_exc, self.Q, self.lam,
                                         self.P0, self.ready_cycles, self.info_min)
        self.k_hat, self.exc_prev = common.K_MAX, 0.

    def step(self, o):
        row = o.log[o.n-1]
        t, F, sent = (float(row[i]) for i in (0, 5, 11))
        self.k_hat = self.estimator.update(F, self.displacement.update(o))
        exc = self.A_E*math.sin(2*math.pi*self.f_exc*t)
        dy = self.K_I/self.k_hat*(common.FD-F)*self.Tc
        cmd = sent+dy+exc-self.exc_prev
        self.exc_prev = exc
        self.trace.append((t, self.k_hat))
        return float(cmd)


class WijType(common.DirectPolicy):
    """Noisy privileged-stiffness PI surrogate (带噪特权刚度 PI 代理).

    Task adaptation inspired by Wijayarathne ICRA2020 IV, Fig. 4.

    Replace the dedicated probe's spring/damper identification (10)-(13)
    with k_eff*(1+noise_std*N(0,1)), independent of robot excitation and
    surface motion. This privileged advantage is unavailable to causal
    force/encoder-only methods. Measure at probe_hz, holding between samples.
    Scalar unit-regressor RLS substitutes for the two-parameter (14)-(16).
    In (17), normalize innovation by current stiffness for dimensionless
    tanh; clip to [lambda_min,lambda_max], then use it as the conventional
    RLS denominator forgetting factor when forget_mode='tanh'. Alternatively
    forget_mode='constant' uses standard constant forgetting lam in [0.8,0.999].
    probe_hz is tunable in [20,125] Hz; samples occur on control ticks (thus
    noninteger period ratios have tick-quantized intervals). This is an adaptation rather than
    the original covariance map f. PI uses incremental form
    C*(K_P*Delta e+K_I*e*Tc). No motion FF or robot dither. P0 is scalar
    covariance in the unit-regressor, unit-measurement-variance convention.
    """
    name = 'WijType'
    params = dict(K_P=.2, K_I=30., probe_hz=20., noise_std=.05,
                  lambda_min=.05, lambda_max=.99, innovation_scale=1.,
                  P0=100., seed=0, forget_mode='tanh', lam=.99)

    def __init__(self, *, K_P=.2, K_I=30., probe_hz=20., noise_std=.05,
                 lambda_min=.05, lambda_max=.99, innovation_scale=1.,
                 P0=100., seed=0, forget_mode='tanh', lam=.99):
        self.params = dict(K_P=K_P, K_I=K_I, probe_hz=probe_hz,
                           noise_std=noise_std, lambda_min=lambda_min,
                           lambda_max=lambda_max, innovation_scale=innovation_scale,
                           P0=P0, seed=seed, forget_mode=forget_mode, lam=lam)
        _positive(K_I=K_I, probe_hz=probe_hz, innovation_scale=innovation_scale, P0=P0)
        if not 0 < lambda_min <= lambda_max <= 1:
            raise ValueError('Require 0 < lambda_min <= lambda_max <= 1')
        if not 20 <= probe_hz <= 125:
            raise ValueError('probe_hz must be in [20, 125] Hz')
        if forget_mode not in ('tanh', 'constant') or not .8 <= lam <= .999:
            raise ValueError("Require forget_mode in {'tanh','constant'} and lam in [0.8,0.999]")
        if not np.isfinite(K_P) or K_P < 0 or not np.isfinite(noise_std) or noise_std < 0:
            raise ValueError('K_P and noise_std must be finite and nonnegative')
        self.__dict__.update(self.params)

    def reset(self, platform):
        super().reset(platform)
        self.rng = np.random.default_rng(self.seed)
        self.k_hat, self.cov, self.next_probe = common.K_MAX, self.P0, 0.
        self.e_prev = None

    def step(self, o):
        row = o.log[o.n-1]
        t, F, sent = (float(row[i]) for i in (0, 5, 11))
        if t+1e-10 >= self.next_probe:
            measurement = float(row[8])*(1+self.noise_std*self.rng.standard_normal())
            innovation = measurement-self.k_hat
            lam = self.lam if self.forget_mode == 'constant' else float(np.clip(math.tanh(abs(innovation)/
                        (max(self.k_hat, common.K_MIN)*self.innovation_scale)),
                        self.lambda_min, self.lambda_max))
            gain = self.cov/(lam+self.cov)
            self.k_hat = float(np.clip(self.k_hat+gain*innovation, common.K_MIN, common.K_MAX))
            self.cov = (1-gain)*self.cov/lam
            self.next_probe += (math.floor((t-self.next_probe)*self.probe_hz+1e-9)+1)/self.probe_hz
        e = common.FD-F
        de = 0. if self.e_prev is None else e-self.e_prev
        self.e_prev = e
        self.trace.append((t, self.k_hat))
        return float(sent+(self.K_P*de+self.K_I*e*self.Tc)/self.k_hat)


class DOMPAC(common.DirectPolicy):
    """Shen CEP2027 S03 Eq. (17)-(29), scalar normal-axis task adaptation.

    Paper coordinates: e=F-FD, outward x=-y, stationary reference x_d.
    (17)-(18): e''=b0*u+f+D; b0=-k_h0/m, f=-(k_h0/m)*(e+b*v), v=y_cmd'.
    (19)-(21): ESO z=[e_hat,e_dot_hat,D_hat], gains [3w,3w^2,w^3].
      Exact ZOH of continuous linear ESO, newest measured e held over the
      update; known input is previous interval's average f+b0*u. Initialize
      z1 at first e, z2=z3=0 to avoid initial peaking. D includes mismatch,
      surface motion, servo/delay and omitted nonlinear effects.
    (22)-(23): u=(u0-f-z3)/b0, giving e''=u0+D-z3.
    (24): double integrator A=[[1,T],[0,1]], B=[T^2/2,T] (exact ZOH).
      The text says forward Euler but prints T^2/2: use the latter.
    (25)-(27): sum Q*e[i]^2+R*u0[i]^2 over i=1..Np / 0..Np-1,
      with no terminal cost. Eq. (25) penalizes u0, not Delta u0.
      Condense prediction and solve once in reset; online take first row.
    (28)-(29): substitute u in m*v'+b*v=-e-u. Integrate held e,u
      exactly for velocity and position; add increment to actually sent
      command. Reconcile velocity only when plant rate clipping changes the
      requested increment, removing optional dither before reconciliation.
      This prevents hidden command windup without bypassing plant VMAX.

    Defaults m=3.5,b=1000,k_h0=1500,Np=15,w=8 use the hardware setting.
    Tune only q_over_r=Q/R, fixing R=1 (a common cost scaling does not
    change the gain), and Np in {8,15,25,40}. Default q_over_r=20 preserves
    the former Q=10,R=.5 gain for comparison. This legacy default is slower
    than space_patch's task-duration-based pole screen permits; use that
    screened space for equal-budget retuning, not these legacy weights.
    Spring k=0. reset verifies Schur stability of A-BK at the actual Tc.
    There is no reference-motion feed-forward. Remark 4 mentions execution
    bounds but gives no numeric values: no extra u0/displacement thresholds
    are invented; the plant's speed bound remains in effect. MPC neglects
    servo and delay, so this is no stability certificate for this platform.
    """
    name = 'DOMPAC'
    params = dict(m=3.5, b=1000., k_h0=1500., eso_bandwidth=8., q_over_r=20., Np=15)

    def __init__(self, *, m=3.5, b=1000., k_h0=1500., eso_bandwidth=8.,
                 q_over_r=20., Np=15):
        self.params = dict(m=m, b=b, k_h0=k_h0, eso_bandwidth=eso_bandwidth,
                           q_over_r=q_over_r, Np=Np)
        _positive(m=m, b=b, k_h0=k_h0, eso_bandwidth=eso_bandwidth, q_over_r=q_over_r)
        if isinstance(Np, bool) or int(Np) != Np or Np not in (8, 15, 25, 40):
            raise ValueError('Np must be one of {8, 15, 25, 40}')
        self.__dict__.update(self.params)
        self.Np = int(Np)

    def reset(self, platform):
        super().reset(platform)
        T, w = self.Tc, self.eso_bandwidth
        self.Q, self.R = self.q_over_r, 1.
        self.mpc_gain, self.nominal_poles = nominal_mpc(T, self.Np, self.q_over_r)
        if np.max(np.abs(self.nominal_poles)) >= 1.:
            raise ValueError('Nominal compensated MPC loop is not Schur-stable')
        L = np.array([3*w, 3*w*w, w**3])
        Ac = np.array([[0., 1., 0.], [0., 0., 1.], [0., 0., 0.]])
        Ac[:, 0] -= L
        aug = np.zeros((5, 5))
        aug[:3, :3], aug[:3, 3], aug[1, 4] = Ac, L, 1.
        disc = expm(aug*T)
        self.eso_A, self.eso_B = disc[:3, :3], disc[:3, 3:]
        self.z = None
        self.k_hat = self.k_h0
        self.velocity, self.held_input = 0., 0.
        self.sent_prev = None
        self.exc_prev, self.exc_delta_prev = 0., 0.
        self.decay = math.exp(-self.b*T/self.m)

    def _stiffness(self, o):
        return self.k_h0

    def _excitation(self, t):
        return 0.

    def step(self, o):
        row = o.log[o.n-1]
        t, F, sent = (float(row[i]) for i in (0, 5, 11))
        e = F-common.FD
        self.k_hat = self._stiffness(o)
        if self.sent_prev is not None:
            realized = sent-self.sent_prev-self.exc_delta_prev
            if abs(realized-self.dy_prev) > 1e-12:
                self.velocity = realized/self.Tc
        if self.z is None:
            self.z = np.array([e, 0., 0.])
        else:
            self.z = self.eso_A@self.z+self.eso_B@np.array([e, self.held_input])
        u0 = -float(self.mpc_gain@self.z[:2])
        b0 = -self.k_hat/self.m
        f = b0*(e+self.b*self.velocity)
        u = (u0-f-self.z[2])/b0
        v_inf = (-e-u)/self.b
        dy = v_inf*self.Tc+(self.velocity-v_inf)*self.m/self.b*(1-self.decay)
        v_next = v_inf+(self.velocity-v_inf)*self.decay
        self.held_input = self.k_hat*(v_next-self.velocity)/self.Tc
        self.velocity = v_next
        exc = self._excitation(t)
        self.exc_delta_prev = exc-self.exc_prev
        self.exc_prev = exc
        self.sent_prev, self.dy_prev = sent, dy
        self.trace.append((t, self.k_hat))
        return float(sent+dy+self.exc_delta_prev)


class RossiDOMPAC(DOMPAC):
    """Combination baseline, not a controller claimed by either paper.

    Same estimator/dither as RossiType. Replace k_h0 in S03 (17),(22),
    (28),(29) each tick, recomputing b0 and f by scalar scaling. Preserve
    previous held input in force units when stiffness changes. ESO poles
    depend only on bandwidth; compensated prediction (24)-(27) is a double
    integrator independent of stiffness, so both matrices remain unchanged.
    Add dither after admittance, leaving its effect to the ESO; do not add
    an excitation-cancellation model. This is a task adaptation, with the
    same matched 30 Hz displacement preprocessing and fixed A_E, tunable
    2-30 Hz dither as RossiType; the Rossi experiment instead used a 0.5 Hz
    force-reference sinusoid. bp_Q is the bandpass quality factor, distinct
    from q_over_r=Q/R in the predictive cost (R=1).
    k_h0 is retained for constructor compatibility but the conservative
    Rossi estimate replaces it from the first observation onward.
    """
    name = 'RossiDOMPAC'
    params = dict(DOMPAC.params, lam=.99, bp_Q=3., f_exc=20.8,
                  P0=100., ready_cycles=3., info_min=.01, A_E=2e-5)

    def __init__(self, *, m=3.5, b=1000., k_h0=1500., eso_bandwidth=8.,
                 q_over_r=20., Np=15, lam=.99, bp_Q=3.,
                 f_exc=20.8, P0=100., ready_cycles=3., info_min=.01, A_E=2e-5):
        super().__init__(m=m, b=b, k_h0=k_h0, eso_bandwidth=eso_bandwidth,
                         q_over_r=q_over_r, Np=Np)
        self.params.update(lam=lam, bp_Q=bp_Q, f_exc=f_exc,
                           P0=P0, ready_cycles=ready_cycles, info_min=info_min, A_E=A_E)
        _positive(A_E=A_E)
        if not 2 <= f_exc <= 30:
            raise ValueError('f_exc must be in [2, 30] Hz')
        self.__dict__.update(self.params)

    def reset(self, platform):
        super().reset(platform)
        self.displacement = _MatchedDisplacement(self)
        self.estimator = _RossiEstimator(self.Tc, self.f_exc, self.bp_Q, self.lam,
                                         self.P0, self.ready_cycles, self.info_min)
        self.k_hat = common.K_MAX

    def _stiffness(self, o):
        return self.estimator.update(float(o.log[o.n-1, 5]), self.displacement.update(o))

    def _excitation(self, t):
        return self.A_E*math.sin(2*math.pi*self.f_exc*t)


class YuenType(common.DirectPolicy):
    """Route-3 privileged motion sensor, adapted to the position interface.

    Yuen T-RO2010 (2),(12),(13), V supplies Fourier EKF and motion FF;
    Cao TIM2025 (8)-(11) supplies residual-triggered refitting. This is a
    task adaptation, not reproduction of the torque-controlled MCI. Here
    velocity/acceleration FF becomes predicted-position increments and PI:
      dy=(K_P*Delta e+K_I*e*Tc)/k_h0-Delta h_pred, e=FD-F.
    Minus height is necessary because penetration=y+surface_height.
    K_P is dimensionless, K_I is s^-1, k_h0 is N/m. They are not Yuen's
    torque gains Kf,Kv; no exact gain equivalence is claimed. The space
    includes his nominal ke=133 N/m and task-appropriate feedback scales.

    Favorable assumptions: undeformed local surface is visible throughout
    contact (Yuen stops updating), exact registration/latency, no occlusion
    or contact deformation. p_scene is copied and used ONLY by the sensor's
    surface() call; no stiffness, derivative, frequency or future truth is
    read. Fixed sensor: 30 Hz, 30 ms, 0.1 mm Gaussian noise; seed 0 gives
    common noise realizations across candidates. Acquisition k/fs is used
    only after k/fs+tau_m<=current tick. Acquisition x is interpolated from
    the available encoder prefix (not a future robot path).

    EKF state [offset,a1,b1,...,am,bm,omega] uses rotating harmonic pairs,
    equivalent to (13)'s amplitude/phase form. Internal positions are mm.
    Q=q_scale*diag([.05 mm]^2,...,[.04 rad/s]^2)*dt, R=(.1 mm)^2;
    random walks allow amplitude/phase/frequency changes. Frequency bounds
    .1-2 Hz cover respiration to heart; .3 Hz is a fixed initial prior.
    After 2 s, causal least-squares frequency-grid initialization uses up
    to 5 s of arrived samples. Innovations above refit_threshold [m] reset
    the EKF. Grid refits are spaced >=.5 s; intervening resets reanchor the
    offset to the arrived measurement (Cao's measured-position fallback).
    No pre-episode data or scene-truth initialization is supplied.

    Target=t+(ndel+.5)*Tc+tau_a. From the latest acquisition stamp the
    horizon is (ndel+.5)*Tc+tau_a+tau_m plus sample/tick age. Latency is
    compensated once. Successive predictions include EKF corrections;
    the first establishes an offset without a command jump. Incremental
    PI starts from actually sent command, avoiding hidden position windup.
    No stiffness adaptation, extra excitation or contact supervisor.
    """
    name = 'YuenType'
    params = dict(k_h0=1500., K_P=.2, K_I=30., harmonics=2,
                  q_scale=1., refit_threshold=1e-3)
    f_s, tau_m, sigma_m = 30., .030, .0001

    def __init__(self, p_scene=None, *, k_h0=1500., K_P=.2, K_I=30.,
                 harmonics=2, q_scale=1., refit_threshold=1e-3):
        _positive(k_h0=k_h0, K_I=K_I, q_scale=q_scale, refit_threshold=refit_threshold)
        if not np.isfinite(K_P) or K_P < 0:
            raise ValueError('K_P must be finite and nonnegative')
        if isinstance(harmonics, bool) or harmonics not in (1, 2, 3):
            raise ValueError('harmonics must be one of {1,2,3}')
        self.params = dict(k_h0=k_h0, K_P=K_P, K_I=K_I, harmonics=int(harmonics),
                           q_scale=q_scale, refit_threshold=refit_threshold)
        self.__dict__.update(self.params)
        self.p_scene = None if p_scene is None else np.array(p_scene, float, copy=True)

    def reset(self, platform):
        super().reset(platform)
        if self.p_scene is None:
            raise ValueError('YuenType requires p_scene for its simulated motion sensor')
        self.rng = np.random.default_rng(0)
        self.sample_index = 0; self.sensor_trace = []
        self.state = None; self.sample_t = None
        self.e_prev = self.pred_prev = None
        self.last_grid = -np.inf; self.initial_fit = False; self.refits = 0
        self.innovation = np.nan
        d = 2*self.harmonics+2
        self.H = np.zeros(d); self.H[0] = 1.; self.H[1:-1:2] = 1.
        self.Qdiag = np.full(d, .05**2); self.Qdiag[-1] = .04**2
        self.R = (self.sigma_m*1000)**2

    def _covariance(self):
        diag = np.full(len(self.H), 5.**2); diag[0] = self.R
        diag[-1] = (2*math.pi*.3)**2
        return np.diag(diag)

    def _transition(self, dt):
        z = self.state.copy(); A = np.eye(len(z))
        for i in range(self.harmonics):
            j = 1+2*i; angle = (i+1)*z[-1]*dt
            c, s = math.cos(angle), math.sin(angle)
            a, b = self.state[j:j+2]
            z[j:j+2] = [c*a+s*b, -s*a+c*b]
            A[j:j+2, j:j+2] = [[c, s], [-s, c]]
            A[j:j+2, -1] = [(i+1)*dt*z[j+1], -(i+1)*dt*z[j]]
        return z, A

    def _refit(self, ts, measurement):
        arrived = np.array([(r['sample_t'], r['z']*1000) for r in self.sensor_trace
                            if r['sample_t'] >= ts-5.])
        if ts-self.last_grid >= .5 and len(arrived) >= max(20, 4*self.harmonics+4):
            tt, yy = arrived[:, 0]-ts, arrived[:, 1]; best = None
            for freq in np.linspace(.1, 2., 48):
                angles = 2*math.pi*freq*tt[:, None]*np.arange(1, self.harmonics+1)
                X = np.column_stack((np.ones(len(tt)), np.cos(angles), np.sin(angles)))
                reg = np.eye(X.shape[1])*.01; reg[0, 0] = 1e-8
                coef = np.linalg.solve(X.T@X+reg, X.T@yy)
                score = np.mean((X@coef-yy)**2)+.01*np.sum(coef[1:]**2)/len(tt)
                if best is None or score < best[0]:
                    best = score, freq, coef
            _, freq, coef = best
            self.state[0], self.state[-1] = coef[0], 2*math.pi*freq
            self.state[1:-1:2] = coef[1:1+self.harmonics]
            self.state[2:-1:2] = coef[1+self.harmonics:]
            self.last_grid = ts
        self.state[0] += measurement-float(self.H@self.state)
        self.cov = self._covariance(); self.refits += 1

    def _measure(self, ts, z):
        z *= 1000.
        if self.state is None:
            self.state = np.zeros(len(self.H)); self.state[0] = z
            self.state[-1] = 2*math.pi*.3
            self.cov = self._covariance(); self.sample_t = ts
            return
        dt = ts-self.sample_t
        self.state, A = self._transition(dt)
        self.cov = A@self.cov@A.T+np.diag(self.Qdiag*self.q_scale*dt)
        self.sample_t = ts
        self.innovation = z-float(self.H@self.state)
        if abs(self.innovation) > self.refit_threshold*1000 or (ts >= 2. and not self.initial_fit):
            self._refit(ts, z)
            if ts >= 2.: self.initial_fit = True
        else:
            ph = self.cov@self.H
            gain = ph/(float(self.H@ph)+self.R)
            self.state += gain*self.innovation
            B = np.eye(len(gain))-np.outer(gain, self.H)
            self.cov = B@self.cov@B.T+self.R*np.outer(gain, gain)
        self.state[-1] = float(np.clip(self.state[-1], 2*math.pi*.1, 2*math.pi*2.))

    def step(self, o):
        row = o.log[o.n-1]
        t, F, sent = (float(row[i]) for i in (0, 5, 11))
        while self.sample_index/self.f_s+self.tau_m <= t+1e-12:
            ts = self.sample_index/self.f_s
            prefix = o.log[:o.n]
            x = float(np.interp(ts, prefix[:, 0], prefix[:, 1]))
            z = float(common.P.surface(x, ts, self.p_scene)[0])+self.sigma_m*self.rng.normal()
            self.sensor_trace.append(dict(sample_t=ts, available_t=ts+self.tau_m,
                                          used_t=t, x=x, z=z))
            self._measure(ts, z)
            self.sample_index += 1
        target = t+(self.ndel+.5)*self.Tc+self.tau
        pred = None if self.state is None else float(self.H@self._transition(target-self.sample_t)[0])/1000.
        ff = 0. if pred is None or self.pred_prev is None else -(pred-self.pred_prev)
        self.pred_prev = pred
        e = common.FD-F
        de = 0. if self.e_prev is None else e-self.e_prev
        self.e_prev = e
        self.trace.append(dict(t=t, sample_t=self.sample_t, target_t=target, prediction=pred,
                               innovation_mm=self.innovation,
                               frequency_hz=None if self.state is None else self.state[-1]/(2*math.pi),
                               refits=self.refits))
        return float(sent+(self.K_P*de+self.K_I*e*self.Tc)/self.k_h0+ff)
