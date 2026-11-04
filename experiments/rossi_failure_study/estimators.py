"""Causal stiffness estimators for the Rossi/projection comparison."""
import math
import numpy as np
from scipy import signal

import common
from comparators import _MatchedDisplacement

DT = 0.008
KMIN, KMAX = common.K_MIN, common.K_MAX


def _clip_k(k):
    return float(np.clip(k, KMIN, KMAX)) if np.isfinite(k) else np.nan


def _fit_record(k=np.nan, sigma=np.inf, chi=np.inf, *, degree=3, tw=DT):
    info = 0. if not np.isfinite(sigma) or sigma <= 0 else 1./(sigma*sigma)
    return dict(k=k, beta=np.zeros(degree+1), info=info, chi2=chi,
                valid=0, tw=tw, te=0.)


class ProjectionEstimator:
    """The existing common.window_fit estimator, wrapped unchanged."""
    name = 'Projection'

    def __init__(self, nw=16, degree=3, sigma_scale=1., chi_scale=1.):
        self.nw = int(nw); self.degree = int(degree)
        self.dof = max(self.nw-self.degree-2, 1)
        self.sigma_scale = float(sigma_scale); self.chi_scale = float(chi_scale)
        if self.nw < 1 or self.degree < 0:
            raise ValueError('invalid projection window')

    def reset(self, policy):
        self.fit = _fit_record(degree=self.degree)

    def estimate(self, o, n, policy):
        p = policy.params
        fit = common.window_fit(o.log, n, self.nw, self.degree, p['noise'],
                                policy.ns, policy.eb, policy.afc, policy.ndel,
                                p['enc_y'])
        fit = dict(fit)
        if np.isfinite(fit.get('k', np.nan)) and fit.get('info', 0.) > 0:
            sigma = self.sigma_scale/math.sqrt(fit['info'])
            chi = self.chi_scale*float(fit['chi2'])
            k = _clip_k(fit['k'])
        else:
            sigma, chi, k = np.inf, np.inf, np.nan
        fit['k'] = k; fit['info'] = 0. if not np.isfinite(sigma) else 1./(sigma*sigma)
        fit['chi2'] = chi; self.fit = fit
        return k, sigma, chi


class _StreamingBandPass:
    """Identical-channel causal BP followed by discounted scalar EWRLS."""
    def __init__(self, kind='butterworth', order=2, half=4., Q=3., memory=.4,
                 sigma_scale=1., chi_scale=1.):
        self.kind = kind; self.order = int(order); self.half = float(half)
        self.Q = float(Q); self.memory = float(memory)
        self.sigma_scale = float(sigma_scale); self.chi_scale = float(chi_scale)
        if kind not in ('butterworth', 'r2'):
            raise ValueError(kind)
        if self.order < 1 or self.half <= 0 or self.Q <= 0 or self.memory <= 0:
            raise ValueError('invalid band-pass parameters')

    def reset(self, policy):
        f0 = float(policy.params.get('f_exc', 20.8)); fs = 1./policy.Tc
        if self.kind == 'butterworth':
            lo, hi = f0-self.half, f0+self.half
            if not 0 < lo < hi < .5*fs:
                raise ValueError(f'band outside Nyquist: {lo}, {hi}, fs={fs}')
            self.sos = signal.butter(self.order, [lo, hi], btype='band', fs=fs, output='sos')
            self.zi_f = np.zeros((len(self.sos), 2)); self.zi_y = np.zeros_like(self.zi_f)
            imp = np.zeros(8192); imp[0] = 1.; h = signal.sosfilt(self.sos, imp)
        else:
            self.b, self.a = signal.iirpeak(f0, self.Q, fs=fs)
            self.zi_f = np.zeros(max(len(self.a), len(self.b))-1); self.zi_y = self.zi_f.copy()
            imp = np.zeros(8192); imp[0] = 1.; h = signal.lfilter(self.b, self.a, imp)
        self.filter_noise_gain = float(np.sqrt(np.sum(h*h)))
        self.displacement = _MatchedDisplacement(policy)
        self.settle = int(np.searchsorted(np.cumsum(h*h)/np.sum(h*h), .95))+1   # ticks for the filter to forget its start
        self.dof = max(int(round(self.memory/policy.Tc))-1, 1)
        self.lam = math.exp(-policy.Tc/self.memory)
        self.f_prev = np.nan
        self._restart()

    def _restart(self):
        """Discard filter and regression state; samples out of contact carry no stiffness information."""
        self.zi_f = np.zeros_like(self.zi_f); self.zi_y = np.zeros_like(self.zi_y)
        self.f0 = None; self.y0 = None; self.age = 0
        self.sxy = self.syy = self.sse = self.weight = 0.
        self.fit = _fit_record(tw=self.memory)

    def _filter(self, f, y):
        if self.f0 is None: self.f0, self.y0 = float(f), float(y)
        xf, xy = float(f)-self.f0, float(y)-self.y0
        if self.kind == 'butterworth':
            fb, self.zi_f = signal.sosfilt(self.sos, [xf], zi=self.zi_f)
            yb, self.zi_y = signal.sosfilt(self.sos, [xy], zi=self.zi_y)
        else:
            fb, self.zi_f = signal.lfilter(self.b, self.a, [xf], zi=self.zi_f)
            yb, self.zi_y = signal.lfilter(self.b, self.a, [xy], zi=self.zi_y)
        return float(fb[0]), float(yb[0])

    def estimate(self, o, n, policy):
        row = o.log[n-1]; f = float(row[5]); y = float(self.displacement.update(o))
        in_contact = min(f, self.f_prev) > .5      # same contact test as common._build
        self.f_prev = f
        if not in_contact:
            self._restart()
            return np.nan, np.inf, np.inf
        fb, yb = self._filter(f, y)
        self.age += 1
        if self.age <= self.settle:
            return np.nan, np.inf, np.inf
        self.sxy = self.lam*self.sxy + yb*fb; self.syy = self.lam*self.syy + yb*yb
        k = _clip_k(self.sxy/self.syy) if self.syy > 1e-24 else np.nan
        residual = fb-k*yb if np.isfinite(k) else fb
        self.sse = self.lam*self.sse + residual*residual; self.weight = self.lam*self.weight + 1.
        noise_eff = float(policy.params['noise'])*self.filter_noise_gain
        sigma0 = noise_eff/math.sqrt(max(self.syy, 1e-30)); sigma = self.sigma_scale*sigma0
        chi = self.chi_scale*self.sse/(max(noise_eff*noise_eff, 1e-30)*max(self.weight-1., 1.))
        self.fit = _fit_record(k, sigma, chi, tw=self.memory); self.fit['te'] = float(row[0])
        return k, sigma, chi


class BandPassEstimator(_StreamingBandPass):
    name = 'BandPass'
    def __init__(self, order=2, half=4., memory=.4, sigma_scale=1., chi_scale=1.):
        super().__init__('butterworth', order=order, half=half, memory=memory,
                         sigma_scale=sigma_scale, chi_scale=chi_scale)


class R2Estimator(_StreamingBandPass):
    name = 'R2'
    def __init__(self, Q=3., memory=.4, sigma_scale=1., chi_scale=1.):
        super().__init__('r2', Q=Q, memory=memory, sigma_scale=sigma_scale,
                         chi_scale=chi_scale)


# Offline attribution helpers.  Only R2, Butterworth, Kaiser and Projection
# families are kept here; local-linear and closed-loop code is intentionally absent.
def matched_displacement_series(log, policy_like):
    d = _MatchedDisplacement(policy_like); out = np.empty(len(log))
    class O: pass
    for n in range(1, len(log)+1):
        o = O(); o.log = log; o.n = n; out[n-1] = d.update(o)
    return out


def _offline_bp(x, f0, kind, order=2, half=4., Q=3.):
    x = np.asarray(x, float); x0 = x[0] if len(x) else 0.
    if kind == 'R2':
        b, a = signal.iirpeak(f0, Q, fs=125.); return signal.lfilter(b, a, x-x0)
    if kind.startswith('B'):
        sos = signal.butter(order, [f0-half, f0+half], btype='band', fs=125., output='sos')
        return signal.sosfilt(sos, x-x0)
    if kind.startswith('K'):
        n, adb = kind[1:].split('_'); n, adb = int(n), float(adb)
        beta = signal.kaiser_beta(adb); w = signal.windows.kaiser(n, beta)
        h = w*np.cos(2*np.pi*f0*DT*(np.arange(n)-(n-1)/2)); h -= w*(h.sum()/w.sum())
        h /= abs(np.sum(h*np.exp(-2j*np.pi*f0*DT*np.arange(n))))
        return signal.lfilter(h, [1.], x-x0)
    raise ValueError(kind)


def offline_rossi(log, name, memory=.4, f0=20.8, order=2, half=4., Q=3.):
    F = np.asarray(log[:, 5]); y = np.asarray(log[:, 2])
    class P:
        ns = 8; eb = math.exp(-.001/.019); afc = 1-math.exp(-2*math.pi*30*.001); ndel = 1
    yf = matched_displacement_series(log, P())
    fb = _offline_bp(F, f0, name, order, half, Q); yb = _offline_bp(yf, f0, name, order, half, Q)
    lam = math.exp(-DT/memory); sxy = syy = 0.; out = np.full(len(log), np.nan)
    for i, (a, b) in enumerate(zip(fb, yb)):
        sxy = lam*sxy+b*a; syy = lam*syy+b*b
        if syy > 1e-24: out[i] = _clip_k(sxy/syy)
    return out


def offline_projection(log, degree=3, nw=16):
    class P:
        ns = 8; eb = math.exp(-.001/.019); afc = 1-math.exp(-2*math.pi*30*.001); ndel = 1
    out = np.full(len(log), np.nan)
    for n in range(nw+3, len(log)):
        fit = common.window_fit(log, n, nw, degree, .04, P.ns, P.eb, P.afc, P.ndel, True)
        if np.isfinite(fit['k']): out[n-1] = _clip_k(fit['k'])
    return out


def offline_named(log, name):
    if name.startswith('P'):
        d, w = name[1:].split('_'); return offline_projection(log, int(d), int(w))
    parts = name.split('|'); kind = parts[0]
    if kind == 'R2': return offline_rossi(log, kind, float(parts[1]))
    if kind.startswith('B'):
        _, o, h = kind[0], kind[1:].split('_')[0], kind[1:].split('_')[1]
        return offline_rossi(log, kind, float(parts[1]), order=int(o), half=float(h))
    if kind.startswith('K'): return offline_rossi(log, kind, float(parts[1]))
    raise ValueError(name)


OFFLINE_GRID = (['R2|0.8', 'R2|0.2'] +
                [f'B{o}_{h}|{m}' for o in (1, 2, 3, 4) for h in (2, 4, 8, 12) for m in (.1, .4, 1.)] +
                [f'K{n}_{a}|{m}' for n in (32, 48, 64) for a in (40, 60, 80) for m in (.1, .4, 1.)] +
                [f'P{d}_{w}' for d in (2, 3, 4, 5) for w in (12, 16, 24, 32, 40)])
