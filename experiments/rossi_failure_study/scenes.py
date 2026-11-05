"""Scan scenes with drives chosen for physical reasons, not to fit any estimator.

Forms with a source:
  breathing     cos^(2n) cycles with cycle-to-cycle variation (Lujan 1999 form; variability reported in Med Phys 2005)
  voluntary     minimum-jerk set-point (Flash and Hogan 1985), a quintic passed through body dynamics
  body          second-order support with resonance 7-10 Hz (supine apparent mass, Huang and Griffin)
  tissue memory three Prony branches giving a weak power law, roughly a factor 2 over 0.1-100 Hz (liver data)
Assumptions, no direct measurement found, always swept: flinch and cough force and timing (only startle EMG
latency and duration are sourced), reflex stiffening, damping ratio, robot tool stiffness and ripple, rib geometry.
"""
import math
import numpy as np

DX = .001
DT_DRIVE = .001


def _field(rng, length, corr, sigma):
    n = int(length/DX)+1
    w = rng.standard_normal(n+200)
    kern = np.exp(-.5*(np.arange(-100, 101)*DX/corr)**2)
    f = np.convolve(w, kern/np.sqrt((kern**2).sum()), mode='same')[100:100+n]
    return sigma*f


def _breathing(rng, t, amp, period, n_pow):
    out = np.zeros(len(t)); now = -rng.uniform(0, period); a0 = None
    while now < t[-1]:
        tau = period*rng.uniform(.85, 1.15); a = amp*rng.uniform(.8, 1.2)
        m = (t >= now) & (t < now+tau)
        # Exhale pause at zero, rise to amp: this is the even cosine power form, one cycle at a time.
        out[m] = a*(1-np.cos(np.pi*(t[m]-now)/tau)**(2*n_pow))
        now += tau
    return out


def _minjerk(t, t0, dur, amp, hold):
    a = np.clip((t-t0)/dur, 0, 1); b = np.clip((t-t0-dur-hold)/dur, 0, 1)
    s = lambda x: x**3*(10-15*x+6*x**2)
    return amp*(s(a)-s(b))


def _twitch(t, t0, rise, decay):
    """Muscle activation: fast rise, slower decay, non-smooth at onset."""
    x = np.maximum(t-t0, 0.)
    y = (1-np.exp(-x/rise))*np.exp(-x/decay)
    return y/max(y.max(), 1e-9)


def make(cond, tissue, rng, T=10., v=None, breath=None, overrides=None):
    """cond: rest, breath, shift_slow, shift_fast, flinch, cough, tremor, rib, combo; tissue: soft, mid, hard, mixed."""
    o = dict(overrides or {})
    t = np.arange(int(T/DT_DRIVE)+2)*DT_DRIVE
    speed = float(v if v is not None else rng.uniform(.01, .03))
    length = speed*T+.1
    base = {'soft': 1.0e4, 'mid': 8e4, 'hard': 5e5, 'mixed': 1.5e4, 'relaxed': 7e3}[tissue]*math.exp(rng.normal(0, .2))
    ln_e = math.log(base)+_field(rng, length, rng.uniform(.01, .03), o.get('field_sigma', .2))
    h_s = _field(rng, length, rng.uniform(.03, .06), rng.uniform(.0003, .002))
    h_s = h_s-h_s[0]
    ribs = []
    if cond in ('rib', 'combo') or o.get('ribs'):
        ratio = o.get('rib_ratio', rng.uniform(5, 30))
        x = rng.uniform(.012, .03)
        grid = np.arange(len(ln_e))*DX
        while x < length-.02:
            hw = rng.uniform(.005, .012); edge = rng.uniform(.001, .004); cover = rng.uniform(.5, 1.)
            w = 1/(1+np.exp(-(grid-x+hw)/edge))*(1-1/(1+np.exp(-(grid-x-hw)/edge)))
            ln_e = ln_e+w*cover*math.log(ratio)
            h_s = h_s+w*rng.uniform(0, .001)
            ribs.append((x, hw, edge))
            x += rng.uniform(.022, .04)
    u = np.zeros(len(t)); fext = np.zeros(len(t)); act = np.zeros(len(t)); rip = np.zeros(len(t))
    events = []
    sev = o.get('severity', 1.)
    def add_cardiac(a=None):
        amp = a if a is not None else rng.uniform(5e-5, 2e-4)
        hr = rng.uniform(1., 1.4); now = rng.uniform(0, 1/hr)
        pulse = np.zeros(len(t))
        while now < T:
            x = t-now
            m = (x >= 0) & (x < .25)
            pulse[m] += np.sin(2*np.pi*x[m]/.25)*np.exp(-x[m]/.08)
            now += 1/hr*rng.uniform(.95, 1.05)
        return amp*pulse
    if cond != 'rest':
        u = u+add_cardiac()
    if cond in ('breath', 'combo') or o.get('breath') or cond in ('shift_slow', 'shift_fast', 'flinch', 'cough', 'tremor'):
        amp = breath if breath is not None else o.get('breath_amp', rng.uniform(.001, .008))
        u = u+_breathing(rng, t, amp, rng.uniform(3., 6.), int(rng.integers(1, 4)))
    t_ev = rng.uniform(2.5, 4.)
    n_ev = int(rng.integers(1, 3))
    event_kinds = {'shift_slow': ['slow'], 'shift_fast': ['fast'], 'flinch': ['flinch'], 'cough': ['cough'],
                   'combo': ['slow', 'fast', 'flinch', 'cough']}.get(cond, [])
    if event_kinds:
        for i in range(n_ev):
            kind = str(rng.choice(event_kinds))
            sign = rng.choice([-1., 1.])
            if kind == 'slow':
                dur = rng.uniform(.5, 1.); amp = sign*rng.uniform(.002, .015)
                u = u+_minjerk(t, t_ev, dur, amp, rng.uniform(.3, 1.))
                events.append((kind, t_ev, t_ev+dur*2+1.))
            elif kind == 'fast':
                dur = rng.uniform(.15, .4); amp = sign*rng.uniform(.002, .012)
                u = u+_minjerk(t, t_ev, dur, amp, rng.uniform(.3, 1.))
                events.append((kind, t_ev, t_ev+dur*2+1.))
            elif kind == 'flinch':
                rise = rng.uniform(.03, .1); decay = rng.uniform(.1, .4); amp = sign*sev*rng.uniform(.001, .005)
                y = _twitch(t, t_ev, rise, decay)
                u = u+amp*y
                fext = fext+sign*sev*rng.uniform(5., 40.)*y
                act = np.maximum(act, y*rng.uniform(.3, 1.))
                events.append((kind, t_ev, t_ev+decay*3))
            else:
                n_b = int(rng.integers(2, 4)); tb = t_ev
                for _ in range(n_b):
                    y = _twitch(t, tb, rng.uniform(.02, .05), rng.uniform(.05, .12))
                    fext = fext+sev*rng.uniform(15., 80.)*y
                    act = np.maximum(act, y*rng.uniform(.2, .6))
                    u = u+sev*rng.uniform(.0005, .003)*y
                    tb += rng.uniform(.25, .5)
                events.append((kind, t_ev, tb+.3))
            t_ev += rng.uniform(1.5, 3.)
    if cond == 'tremor' or o.get('tremor'):
        # Muscle tremor while the patient holds tension; narrow band close to the probe frequency range.
        f0 = o.get('tremor_f', rng.uniform(8., 12.))
        phase = np.cumsum(2*np.pi*(f0+rng.normal(0, .3, len(t)))*DT_DRIVE)
        env = _minjerk(t, 2., .5, 1., 5.)
        u = u+o.get('tremor_amp', rng.uniform(5e-6, 5e-5))*env*np.sin(phase)
        events.append(('tremor', 2.5, 7.5))
    if o.get('ripple_amp'):
        # Robot transmission ripple: unmeasured by the encoder, frequency set by joint speed; magnitude is an assumption.
        rip = o['ripple_amp']*np.sin(2*np.pi*o.get('ripple_f', 20.8)*t+rng.uniform(0, 6.28))
    fn = rng.uniform(7., 10.)
    m = o.get('m', rng.uniform(15., 40.))
    k = m*(2*np.pi*fn)**2
    zeta = o.get('zeta', rng.uniform(.25, .5))
    return dict(cond=cond, tissue=tissue, T=T, v=speed, ln_e=ln_e, h_s=h_s, dx=DX, u=u, fext=fext, act=act, rip=rip,
                alpha=np.array(o.get('alpha', [.35, .35, .35])), tau=np.array(o.get('tau', [.008, .08, .8])),
                m=m, c=2*zeta*m*(2*np.pi*fn), k=k, kr=o.get('kr', 10**rng.uniform(4.7, 5.3)), lam=rng.uniform(.1, .4), beta=o.get('beta', rng.uniform(.3, 1.)),
                events=events, ribs=ribs, noise_scale=o.get('noise_scale', 1.))
