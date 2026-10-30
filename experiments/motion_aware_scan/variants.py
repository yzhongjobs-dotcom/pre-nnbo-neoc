"""Motion tracking variants; import with motion_aware_scan on sys.path.

build('Proposed-IMM', tuned_callback, q_smooth=...) or build(name, params_dict).
Only causal columns t, y and filtered force are used. The base stiffness window,
admittance, excitation and displacement feed-forward sign remain unchanged.
New defaults are starting points, not closed-loop tuned performance claims.
"""
import json
import math
from pathlib import Path

import numpy as np

import common
import proposed
from imm import IMMTracker


DEFAULTS = {
    'Proposed-IMM': dict(q_smooth=1e-5, q_maneuver=0.1, p_stay=0.98, r_scale=4.),
    'Proposed-KF': dict(q_kf=1e-3, r_scale=4.),
    'Proposed-IMMg': dict(q_smooth=1e-5, q_maneuver=0.1, p_stay=0.98, r_scale=4.),
    'Proposed-KFg': dict(q_kf=1e-3, r_scale=4.),
}
GATED = ('Proposed-IMMg', 'Proposed-KFg')
SINGLE = ('Proposed-KF', 'Proposed-KFg')
# Continuous white-jerk PSD, m^2/s^5; smooth and maneuver ranges do not overlap.
BOUNDS = {'q_smooth': (1e-9, 1e-3), 'q_maneuver': (1e-2, 10.),
          'q_kf': (1e-9, 10.), 'p_stay': (.9, .999), 'r_scale': (.1, 100.)}


class MotionAwareIMM(proposed.MotionAware):
    name = 'Proposed-IMM'

    def __init__(self, q_smooth=1e-5, q_maneuver=0.1, p_stay=.98,
                 r_scale=4., q_kf=1e-3, single_model=False, gate=False, **kw):
        values = (q_smooth, q_maneuver, p_stay, r_scale, q_kf)
        if not all(np.isfinite(v) for v in values):
            raise ValueError('tracker parameters must be finite')
        if min(q_smooth, q_maneuver, q_kf) < 0 or r_scale <= 0 or not 0 <= p_stay <= 1:
            raise ValueError('invalid tracker parameters')
        super().__init__(**kw)
        self.single_model = bool(single_model)
        self.motion_params = (dict(q_kf=q_kf, r_scale=r_scale) if single_model else
                              dict(q_smooth=q_smooth, q_maneuver=q_maneuver,
                                   p_stay=p_stay, r_scale=r_scale))
        self.params.update(self.motion_params)
        self.gate = bool(gate)
        self.name = ('Proposed-KF' if single_model else 'Proposed-IMM')+('g' if gate else '')

    def reset(self, platform):
        super().reset(platform)
        p = self.motion_params
        q = p['q_kf'] if self.single_model else p['q_smooth']
        self.tracker = IMMTracker(Tc=self.Tc, q_smooth=q,
                                  q_maneuver=p.get('q_maneuver', q),
                                  p_stay=p.get('p_stay', 1.),
                                  single_model=self.single_model)
        # Force samples in the plant's last interval are at m=0,...,ns-1 ms.
        # Match their exponential weights using a servo-shaped interpolation
        # between the TWO observed encoder endpoints; no future/true state read.
        # This is more accurate than treating the endpoint as the interval mean,
        # especially at the deliberate 20.8 Hz excitation frequency.
        m = np.arange(self.ns)
        if self.tau > 0:
            fraction = -np.expm1(-m*.001/self.tau)/(-math.expm1(-self.Tc/self.tau))
        else:
            fraction = (m > 0).astype(float)
        weights = self.afc*(1-self.afc)**(self.ns-1-m)
        self._encoder_fraction = float(weights @ fraction/weights.sum())
        self._clear_motion()

    def _clear_motion(self):
        self.Z = 0.
        self._motion_n = None
        self._fin_prev = None
        self._y_prev = None
        self._motion_prediction = None
        self.tracker.initialized = False
        # Base step also clears this on recovery; explicit here avoids a reset offset.
        self.q_prev = None

    def step(self, o):
        # drift_pred is skipped during approach, so supervise resets here.
        if self.mode == 1 or o.log[o.n-1, 5] < .3:
            self._clear_motion()
        return super().step(o)

    def drift_pred(self, o, n):
        if self._motion_n == n:
            return self._motion_prediction
        if self._motion_n is not None and n < self._motion_n:
            raise ValueError('motion observations must arrive in increasing tick order')
        previous_n = self._motion_n
        self._motion_n = n
        self._motion_prediction = None
        if n < 2:
            return None
        j = n-1
        F, Fprev = o.log[j, 5], o.log[j-1, 5]
        y0, y1 = o.log[j-1, 2], o.log[j, 2]
        y = float(y0+self._encoder_fraction*(y1-y0))
        decay = (1-self.afc)**self.ns
        denom = 1-decay
        if denom <= 0:
            raise ValueError('positive force filter bandwidth is required')
        fin = (F-decay*Fprev)/denom
        kref = self.k_good
        if not np.isfinite(kref) or kref <= 0:
            kref = self.kbar
        if not np.isfinite(kref) or kref <= 0:
            raise ValueError('a finite positive reference stiffness is required')
        # noise is the inherited force-noise scale [N]; r_scale absorbs inversion,
        # encoder approximation and measurement correlations, not modeled exactly here.
        R = (self.params['noise']/kref)**2*self.motion_params['r_scale']
        if self.gate and self.fit is not None and np.isfinite(self.fit.get('chi2', np.nan)):
            # consistency gate: when the stiffness window is inconsistent (local stiffness changing
            # inside the window), force changes are not attributable to tissue motion; inflate the
            # motion-measurement variance by the normalised residual (innovation-based adaptation)
            R *= max(1., self.fit['chi2']/self.params['chi2_max'])
        self.tracker.R = max(R, np.finfo(float).tiny)
        if self._fin_prev is None or previous_n != n-1:
            self.Z = 0.
            self.tracker.reset(0.)
            self.q_prev = None
        else:
            # Use the SAME current kref for both forces. Never redivide past Z by
            # a new stiffness: stiffness changes alone must not create level jumps.
            self.Z += (fin-self._fin_prev)/kref-(y-self._y_prev)
            self.tracker.update(self.Z)
        self._fin_prev, self._y_prev = float(fin), y
        L = (self.ndel+0.5)*self.Tc
        s, sd = self.tracker.predict(L)
        self._motion_prediction = (None, s+self.tau*sd)
        # Before a first accepted stiffness, the base tuple path multiplies dq
        # by kbar but divides the final command by k_seek (different scales).
        # Track with the requested kbar fallback, but anchor q_prev every tick
        # until k_good exists so this startup mismatch cannot amplify motion FF.
        if not np.isfinite(self.k_good) or self.k_good <= 0:
            self.q_prev = None
        return self._motion_prediction


def tuned(method):
    """Read existing best.params without writes; new methods inherit Proposed.

    If a new method's JSON exists after later tuning, overlay its best.params.
    Absent such a file, the new tracker parameters are documented defaults.
    """
    if method not in ('Proposed', *DEFAULTS):
        raise ValueError(method)
    root = Path(common.HERE)/'tuned'
    with (root/'Proposed.json').open(encoding='utf-8') as f:
        base = dict(json.load(f)['best']['params'])
    if method in DEFAULTS:
        base.update(DEFAULTS[method])
        path = root/(method+'.json')
        if path.exists():
            with path.open(encoding='utf-8') as f:
                base.update(json.load(f)['best']['params'])
    return base


def build(name, tuned=tuned, **overrides):
    """Accept a tuned(method) callback or a parameter mapping, then overrides.

    Callback always supplies inherited Proposed parameters; a separate variant
    JSON is optional, supporting run_v2's callback before variant tuning exists.
    """
    if name not in ('Proposed', *DEFAULTS):
        raise ValueError(name)
    gate = name in GATED
    if gate:
        overrides = dict(overrides); overrides.setdefault('edge_rule', 'raise')
    if callable(tuned):
        params = dict(tuned('Proposed'))
        if name in DEFAULTS:
            params.update(DEFAULTS[name])
            try:
                extra = tuned(name)
            except (FileNotFoundError, KeyError):
                extra = {}
            params.update(extra)
    else:
        params = dict(DEFAULTS.get(name, {}))
        params.update(tuned)
    params.update(overrides)
    if name == 'Proposed':
        return proposed.MotionAware(**params)
    return MotionAwareIMM(single_model=(name in SINGLE), gate=gate, **params)


def space(name, r):
    """Sample only new parameters; retain tuned Proposed parameters in build()."""
    if name not in DEFAULTS:
        raise ValueError(name)
    return {key: float(math.exp(r.uniform(math.log(BOUNDS[key][0]),
                                          math.log(BOUNDS[key][1]))))
            for key in DEFAULTS[name]}


def perturb(name, params, r):
    """Log-local perturbation of new parameters; inherited parameters are fixed."""
    if name not in DEFAULTS:
        raise ValueError(name)
    out = dict(params)
    for key, default in DEFAULTS[name].items():
        lo, hi = BOUNDS[key]
        value = out.get(key, default)
        if not np.isfinite(value) or value <= 0:
            raise ValueError('perturb requires positive finite parameters')
        out[key] = float(np.clip(value*math.exp(r.normal(0., .35)), lo, hi))
    return out
