"""Paired, deterministic scan scenes; SI units throughout.

Quantities without direct evidence are drawn over their physical range per scene, so no unproven value stands for the real
scene. Discrete factors are balanced within each tissue and regime, and results can be read per level. overrides fix any of them.
"""
import math
import zlib
import json
from pathlib import Path
from functools import lru_cache

import numpy as np
from scipy.optimize import nnls

DX = .0005
X0 = .42
CONDITIONS = ('rest', 'breath', 'shift_slow', 'shift_fast', 'flinch', 'cough', 'tremor', 'rib', 'combo')
# Additive scenes use the original randomization strata. Do not append them to
# CONDITIONS: its nine-position balance and historical seeds must stay intact.
OVERLAY_CONDITIONS = {'breath_rib': 'breath', 'breath_fast': 'shift_fast',
                      'breath_rib_fast': 'shift_fast'}
TISSUES = ('relaxed', 'active', 'rib', 'mixed')
CELLS3 = [(c, s) for c in CONDITIONS for s in TISSUES]
REGIMES = ('nominal', 'slow', 'stop_go', 'reversal', 'fast', 'curved', 'curved_untracked', 'strong', 'external')
PROBES = {
    'convex': dict(length=.069, width=.017, radius=.060),
    'linear': dict(length=.04843, width=.01156, radius=math.inf),
}
REMUS_S1 = json.loads((Path(__file__).with_name('remus_s1.json')).read_text(encoding='utf-8'))
REMUS = {key: value['rangeMean']['coeff'] for key, value in REMUS_S1.items()}
DYNAMIC = json.loads((Path(__file__).with_name('dynamic_truth.json')).read_text(encoding='utf-8'))
SHAPES = ('neural_white', 'neural_grey', 'viscous')
PHASES = (8., 14., 21., 28.)


def sample_curves(rng, tissue, overrides):
    region = str(overrides.get('region', 'R'+str(int(rng.integers(1, 6)))))
    number = region[1:]
    if region not in REMUS_S1:
        raise ValueError(region)
    paired = sorted(set(REMUS_S1['R'+number]) & set(REMUS_S1['A'+number]) - {'rangeMean'})
    subject = str(overrides.get('subject', rng.choice(paired)))
    if subject not in paired:
        raise ValueError('A relaxed/active subject pair is required: '+subject)
    state = 'A' if tissue == 'active' or region.startswith('A') else 'R'
    return state+number, subject, REMUS_S1[state+number][subject], REMUS_S1['A'+number][subject]


def _measured_response(kind, phase_deg, rng, slow_scale=None):
    """The default truth, identified rather than chosen: shape from the Neumann records, phase from the stratum, amount of slow memory varying as it does between their specimens. See derive_dynamic_truth.py."""
    table = DYNAMIC['shapes'][kind]; low = int(math.floor(phase_deg)); top = min(low+1, 40); frac = phase_deg-low
    alpha = (1.-frac)*np.array(table[str(low)]['alpha'])+frac*np.array(table[str(top)]['alpha'])
    tau = np.array(DYNAMIC['tau']); w = 2.*np.pi*20.8
    slow = tau >= .4
    # Across the Neumann specimens the force falls by 12-26 percent from 1 s to 10 s, which this factor spans; the fit sits at the upper end. No abdominal measurement exists.
    alpha[slow] *= float(np.exp(rng.uniform(math.log(.4), math.log(1.2)))) if slow_scale is None else slow_scale
    h = 1j*w*tau/(1.+1j*w*tau)
    lo, hi = .1, 10.
    for _ in range(50):
        m = .5*(lo+hi); a = alpha.copy(); a[~slow] *= m
        if math.degrees(np.angle(1.+np.sum(a*h))) < phase_deg:
            lo = m
        else:
            hi = m
    alpha[~slow] *= m
    return alpha, tau, 0.


def _dynamic_response(kind, phase_deg, rng, slow_scale=None):
    """Returns alpha, tau, phi_hc in degrees for one mechanism at a given total phase at the probe frequency.

    neural_white, neural_grey, viscous   identified or evidenced shapes, the default
    hc, power, kb   earlier hand-chosen laws, each the twin of one candidate (H, W, S). The Neumann records reject hc and kb; kept only for explicit sensitivity runs.
    """
    if kind in DYNAMIC['shapes']:
        return _measured_response(kind, phase_deg, rng, slow_scale)
    w = 2.*np.pi*20.8; target = math.tan(math.radians(phase_deg))
    if kind == 'hc':
        # Exploratory parameters, not identified abdominal relaxation constants.
        alpha = np.array([.35, .35, .35])*np.array([1., 1., 1.])*np.array([1., float(rng.uniform(.1, 1.)), float(rng.uniform(.1, 1.))]); tau = np.array([.008, .08, .8])
        m = 1.+np.sum(alpha*1j*w*tau/(1.+1j*w*tau))
        return alpha, tau, math.degrees(math.atan(max(target*m.real-m.imag, 0.)))
    if kind == 'kb':
        tau = float(rng.uniform(.005, .012)); h = 1j*w*tau/(1.+1j*w*tau)
        alpha = target/(h.imag-target*h.real)
        return np.array([alpha, 0., 0.]), np.array([tau, .08, .8]), 0.
    return _power_response(float(phase_deg))


@lru_cache(maxsize=32)
def _power_response(phase_deg):
    w = 2.*np.pi*20.8
    tau = np.array([.003, .03, .3])
    f = np.geomspace(.3, 100., 60); s = 2j*np.pi*f
    A = np.array([s*t/(1.+s*t) for t in tau]).T; h = 1j*w*tau/(1.+1j*w*tau)
    best = None
    for beta in np.arange(.02, .46, .01):
        G = (f/.3)**beta*np.exp(.5j*beta*np.pi*(1.-np.exp(-f/.3))); wt = 1./np.abs(G)
        a, _ = nnls(np.r_[(A.real.T*wt).T, (A.imag.T*wt).T], np.r_[(G.real-1.)*wt, G.imag*wt])
        gap = abs(float(np.angle(1.+np.sum(a*h)))-math.radians(phase_deg))
        if best is None or gap < best[0]:
            best = (gap, a)
    return best[1], tau, 0.


def _smooth(x):
    x = np.clip(x, 0., 1.)
    return x**3 * (10. - 15.*x + 6.*x*x)


def _level(name, levels, cond, tissue, code, i):
    """Balanced over the 36 positions of nine conditions and four repetitions, within one tissue and regime."""
    perm = np.random.default_rng([zlib.crc32(name.encode()), TISSUES.index(tissue), code]).permutation(36)
    return levels[int(perm[CONDITIONS.index(cond)+9*(i % 4)]) % len(levels)]


def _log_uniform(rng, lo, hi):
    return math.exp(rng.uniform(math.log(lo), math.log(hi)))


def _field(rng, grid, corr, sigma):
    n = int(4*corr/DX)
    white = rng.standard_normal(len(grid)+2*n)
    kernel = np.exp(-.5*(np.arange(-n, n+1)*DX/corr)**2)
    return sigma*np.convolve(white, kernel/np.sqrt(np.sum(kernel**2)), 'same')[n:n+len(grid)]


def _profile(t, regime, rng, o):
    limits = (.001, .005) if regime == 'slow' else ((.04, .06) if regime == 'fast' else (.005, .03))
    speed = float(o.get('speed', rng.uniform(*limits)))
    accel = float(o.get('accel', _log_uniform(rng, .05, .5)))
    ramp = 1.875*speed/accel
    vel = speed*_smooth(t/ramp)
    pauses = []
    if regime == 'stop_go':
        for start in np.linspace(.25, .7, int(rng.integers(2, 4)))*t[-1]:
            hold = float(o.get('pause', rng.uniform(.2, .8)))
            vel *= 1.-_smooth((t-start)/ramp)+_smooth((t-start-ramp-hold)/ramp)
            pauses.append((start+ramp, start+ramp+hold))
    if regime == 'reversal':
        for start in np.linspace(.3, .7, int(rng.integers(1, 3)))*t[-1]:
            vel *= 1.-2.*_smooth((t-start)/(2.*ramp))
    pos = np.r_[0., np.cumsum((vel[1:]+vel[:-1])*.0005)]
    pos += X0+.05-pos.min()
    return pos, speed, accel, pauses


def _cardiac(t, T, rng, amp):
    hr = rng.uniform(1., 1.4); now = rng.uniform(0., 1./hr)
    pulse = np.zeros(len(t))
    while now < T:
        x = t-now
        m = (x >= 0.) & (x < .25)
        pulse[m] += np.sin(2.*np.pi*x[m]/.25)*np.exp(-x[m]/.08)
        now += 1./hr*rng.uniform(.95, 1.05)
    return amp*pulse


def cell_scene3(cond, tissue, regime='nominal', i=0, base=420000, T=10., overrides=None):
    requested_cond = cond
    cond = OVERLAY_CONDITIONS.get(cond, cond)
    if cond not in CONDITIONS or tissue not in TISSUES or regime not in REGIMES:
        raise ValueError((cond, tissue, regime))
    o = dict(overrides or {})
    cell = CONDITIONS.index(cond)*len(TISSUES)+TISSUES.index(tissue)
    # Tracking changes only pose, preserving the entire physical scene and noise seed.
    code = REGIMES.index('curved' if regime == 'curved_untracked' else regime)
    rng = np.random.default_rng([base, cell, code, i])
    t = np.arange(int(round(T*1000.))+2)*.001
    xcmd, speed, accel, pauses = _profile(t, regime, rng, o)
    probe_name = str(o.get('probe', _level('probe', ('convex', 'linear'), cond, tissue, code, i)))
    probe = dict(PROBES[probe_name]); probe.update(o.get('geometry', {}))
    force = float(o.get('force', _level('force', (4., 8., 12.), cond, tissue, code, i)))
    slide = str(o.get('slide', _level('slide', ('long', 'short'), cond, tissue, code, i)))
    law = str(o.get('contact_law', _level('law', ('foundation', 'halfspace'), cond, tissue, code, i)))
    # Bare-material phase at the probe frequency. No abdominal contact measurement exists: 8 is near skin and pancreas, 14-21 human liver by MR elastography, 28 spinal cord by indentation.
    total_phase = float(o.get('phase', _level('phi', PHASES, cond, tissue, code, i)))
    dynamics = str(o.get('dynamics', _level('dyn', SHAPES, cond, tissue, code, i)))
    cover = float(o.get('cover', _level('cover', (.007, .018, .032), cond, tissue, code, i)))
    region, subject, curve, active_curve = sample_curves(rng, tissue, o)
    q0, E0 = 1.4e6, 23.e3     # reference units; absolute forces come from S1
    length = float(xcmd.max()-X0+.08)
    grid = X0+np.arange(int(math.ceil(length/DX))+1)*DX
    hetero = np.exp(_field(rng, grid, float(rng.uniform(.01, .03)), float(o.get('field_sigma', .2))))
    hetero /= hetero.mean()         # the field holds relative variation only; absolute force comes from the measured curve
    if tissue == 'mixed':
        for center in np.linspace(grid[0]+.04, grid[-1]-.04, 4):
            hetero *= np.exp(math.log(float(rng.uniform(2., 8.)))*np.exp(-.5*((grid-center)/.006)**2))      # inclusion in soft tissue; bone is the rib cells
    qscale = float(o.get('q_scale', 1.))
    qfield, efield = hetero.copy(), hetero.copy()
    ribs = []
    rib_gap = float(o.get('rib_gap', rng.uniform(.006, .019)))
    start_rib = float(rng.uniform(.04, .09))
    if tissue == 'rib' or cond in ('rib', 'combo') or requested_cond in ('breath_rib', 'breath_rib_fast'):
        center = grid[0]+start_rib
        while center < grid[-1]:
            # the cover spreads the bone's influence at the surface over about its own thickness
            width = float(rng.uniform(.010, .020)); edge = cover*float(rng.uniform(.5, 1.))
            ribs.append((center, width, edge, cover))
            center += width+rib_gap
    radius_s = float(o.get('surface_radius', rng.uniform(.15, .5))) if regime.startswith('curved') else math.inf
    center_s = float(xcmd[0])
    # Undulation of the body surface: the scan itself drives normal motion through its slope.
    h_s = _field(rng, grid, float(rng.uniform(.03, .06)), float(o.get('surface_amp', rng.uniform(.0003, .002))))
    h_s = h_s-np.interp(center_s, grid, h_s)-(grid-center_s)**2/(2.*radius_s)
    u = np.zeros(len(t)); act = np.zeros(len(t)); events = []; event_info = []; components = []
    breath_pp = float(o.get('breath_pp', rng.uniform(.003, .012))); breath_period = float(o.get('breath_period', rng.uniform(3., 6.)))
    breathing = np.zeros(len(t)); cardiac = np.zeros(len(t))
    if cond not in ('rest', 'rib'):
        # cos^2 cycles have zero slope at both joins; amp is peak-to-peak.
        breathing = .5*breath_pp*(1.-np.cos(2.*np.pi*t/breath_period))
        u += breathing
    if cond != 'rest':
        cardiac = _cardiac(t, T, rng, float(rng.uniform(5.e-5, 2.e-4)))
        u += cardiac
    event_kinds = {'shift_slow': ['slow'], 'shift_fast': ['fast'], 'flinch': ['flinch'],
                   'cough': ['cough', 'cough'], 'combo': ['slow', 'fast', 'flinch', 'cough']}.get(cond, [])
    strong = regime == 'strong'
    for j, kind in enumerate(event_kinds):
        start = .2*T+j*.17*T
        sign = float(rng.choice([-1., 1.]))
        if kind in ('slow', 'fast'):
            amp = float(rng.uniform(.002, .015 if kind == 'slow' else .012))
            dur = float(rng.uniform(.5, 1.) if kind == 'slow' else rng.uniform(.15, .4))
            peak = float(o.get('motion_peak', 1.875*amp/dur))
            dur = 1.875*amp/peak
            hold = .15
            wave = amp*(_smooth((t-start)/dur)-_smooth((t-start-dur-hold)/dur))
            end = start+2.*dur+hold
        else:
            amp = float(rng.uniform(.001, .005) if kind == 'flinch' else rng.uniform(.0005, .003))
            peak = float(o.get('motion_peak', rng.uniform(.11, .14) if strong else rng.uniform(.02, .035)))
            # sin^2 pulse has zero velocity at onset and return; peak speed = pi A / duration.
            dur = np.pi*amp/peak
            phase = np.clip((t-start)/dur, 0., 1.)
            wave = amp*np.sin(np.pi*phase)**2
            end = start+dur
            # The muscle tenses with the movement, stiffening the tissue under the probe for its duration.
            level = float(rng.uniform(.3, 1.) if kind == 'flinch' else rng.uniform(.2, .6))
            act = np.maximum(act, level*np.sin(np.pi*np.clip((t-start)/(dur+.1), 0., 1.))**2)
        wave *= sign
        u += wave; components.append(wave)
        events.append((kind, start, end))
        event_info.append(dict(kind=kind, start=start, end=end, amplitude=amp, peak_limit=peak))
    if cond == 'tremor':
        f0 = float(o.get('tremor_f', rng.uniform(8., 12.)))
        phase = np.cumsum(2.*np.pi*(f0+rng.normal(0., .3, len(t)))*.001)
        u += float(o.get('tremor_amp', rng.uniform(5.e-6, 5.e-5)))*_smooth((t-2.)/.5)*np.sin(phase)
    fext = np.zeros(len(t))
    if regime == 'external':
        # Separate drive interpretation: no duplicate displacement or activation.
        fext = u*float(o.get('external_gain', 5000.)); u = np.zeros_like(u)
    m = float(o.get('m', rng.uniform(15., 40.))); freq = float(rng.uniform(7., 10.))
    k = m*(2.*np.pi*freq)**2
    mu = float(o.get('mu', _log_uniform(rng, .03, .30)))
    mag = rng.normal(2.9, 1.6)
    pose = float(o.get('pose_error', float(rng.choice([-1., 1.]))*max(mag, 0.)))       # in vivo normal error 2.9 +/- 1.6 degrees (mean absolute)
    sc = dict(cond=cond, tissue=tissue, cell=(cond, tissue), regime=regime, T=T,
              probe=probe_name, geometry=probe, FD=force, contact_law=law,
              slide=slide, xcmd=xcmd, v=speed, accel=accel, pauses=pauses,
              field_x0=X0, dx=DX, q_field=qfield, E_field=efield, h_s=h_s,
              q_base=q0, E_base=E0, q_scale=qscale, ribs=np.array(ribs, float).reshape(-1, 4),
              cover=cover, rib_gap=rib_gap, full_cover=bool(o.get('full_cover', False)),
              surface_radius=radius_s, surface_center=center_s,
              tracking=regime != 'curved_untracked', theta0=math.radians(pose),
              u=u, fext=fext, act=np.zeros(len(t)) if regime == 'external' else act,
              rip=np.zeros(len(t)),
              k_rib=float(o.get('k_rib', _log_uniform(rng, 1.e4, 4.5e4))),    # whole cage about 9 N/mm to a single rib about 45 N/mm
              events=events, event_info=event_info,
              motion_components=components, breath_pp=breath_pp, breath_period=breath_period,
              phase=total_phase, dynamics=dynamics,
              m=m, c=2.*float(o.get('zeta', rng.uniform(.25, .5)))*m*2.*np.pi*freq, k=k,
              kr=float(o.get('kr', _log_uniform(rng, 5.e4, 2.e5))), K_t=float(o.get('K_t', _log_uniform(rng, 5.e3, 5.e4))),
              rho=float(o.get('rho', rng.uniform(-.3, .3))),
              m_t=float(o.get('m_t', rng.uniform(1., 3.))), zeta_t=float(o.get('zeta_t', rng.uniform(.3, .7))),
              mu_k=mu, mu_s=mu*float(o.get('mu_ratio', rng.uniform(1., 1.3))), v_s=float(o.get('v_s', rng.uniform(.002, .01))),
              v_e=.0005, k_d=float(o.get('k_d', rng.uniform(0., 5.))), k1=float(o.get('k1', .001)), k2=float(o.get('k2', .003)),
              asymmetry=bool(o.get('asymmetry', True)), l_tool=float(o.get('l_tool', .15)),
              noise_scale=float(o.get('noise_scale', 1.)), sensor_lp=o.get('sensor_lp', None),
              sensor_sigma=float(o.get('sensor_sigma', .1 if o.get('sensor_lp') else .08)))
    sc['alpha'], sc['tau'], sc['phi_hc'] = _dynamic_response(dynamics, total_phase, rng, o.get('slow_scale'))
    if dynamics in DYNAMIC['shapes']:
        # The measured curve is the response at its own loading rate, near .3 Hz, so the equilibrium lies below it. At the start of a scene the probe arrived some time ago.
        sc['equilibrium'] = 1./abs(1.+np.sum(sc['alpha']*2j*np.pi*.3*sc['tau']/(1.+2j*np.pi*.3*sc['tau'])))
        ramp, wait = float(o.get('load_ramp', rng.uniform(1., 5.))), float(o.get('load_wait', _log_uniform(rng, 2., 60.)))
        sc['memory_start'] = sc['alpha']*sc['tau']/ramp*(1.-np.exp(-ramp/sc['tau']))*np.exp(-wait/sc['tau'])
    sc['region'] = region
    sc['subject'] = subject
    # The true system, and what calibration gives the Method. Ranges are assumptions; no measurement for this robot is available.
    sc['tool_mass'] = float(o.get('tool_mass', rng.uniform(.15, .45)))        # probe, holder and cable beyond the sensor
    sc['tool_mass_hat'] = sc['tool_mass']*(1.+float(o.get('mass_error', rng.uniform(-.05, .05))))
    sc['tau_a_true'] = float(o.get('tau_a_true', rng.uniform(.015, .025)))
    sc['tau_a_hat'] = sc['tau_a_true']*(1.+float(o.get('tau_error', rng.uniform(-.05, .05))))
    # A sine in free air calibrates mass and force time offset together: the reading is m a, its phase the offset.
    sc['sensor_delay_ms'] = float(o.get('sensor_delay_ms', rng.uniform(0., 2.)))
    sc['sensor_delay_hat'] = sc['sensor_delay_ms']+float(o.get('delay_error', rng.uniform(-.25, .25)))
    if tissue == 'active':
        sc['act'][:] = 0.
    if not o.get('no_curve', False):
        sc['remus_poly'] = np.array(curve['coeff']); sc['fit_max_mm'] = curve['fitMaxDispl']
        sc['active_poly'] = np.array(active_curve['coeff']); sc['active_fit_max_mm'] = active_curve['fitMaxDispl']
        sc['transfer_gamma'] = .5 if law == 'halfspace' else 1.
    if 'phi_hc' in o:         # direct setting for checks and sensitivity: nominal memory, given rate phase
        sc['alpha'] = np.array([.35]*3)*float(o.get('alpha_scale', 1.)); sc['tau'] = np.array([.008, .08, .8])*float(o.get('tau_scale', 1.)); sc['phi_hc'] = float(o['phi_hc']); sc['dynamics'] = 'hc'
        sc.pop('equilibrium', None); sc.pop('memory_start', None)
    if not abs(sc['rho']) < 1. or min(sc['kr'], sc['K_t'], *sc['tau']) <= 0.:
        raise ValueError('Compliance must be positive definite and time constants positive')
    if 'remus_poly' in sc:
        import physics_plant3 as P
        P.attach_curve(sc)
    if requested_cond in OVERLAY_CONDITIONS:
        sc.update(cond=requested_cond, cell=(requested_cond, tissue),
                  source_cond=cond, motion_breath=breathing, motion_cardiac=cardiac)
    return sc


def make3(cond, tissue, regime, i=0, base=420000, T=10., overrides=None):
    return cell_scene3(cond, tissue, regime, i, base, T, overrides)
