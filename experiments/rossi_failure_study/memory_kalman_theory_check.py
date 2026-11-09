"""Numerical check of the theory in kalman优化idea.txt, before any plant or closed loop.

Part 1, real windows: does the noise model the note proposes match the recorded windows, and is the non-sensor error negligible?
Part 2, a controlled world where only the memory question is alive: tick noise at the measured level, probe force at the measured amplitude,
no d, true gain from four shapes. Same window D4 for every window-based filter. Reports each filter family's frontier, not one tuned point:
for a given plateau noise, how well it follows a change. No truth enters any filter.

Families (all two-mode IMM unless noted):
  LV  level state, window input, R inflated by the long-run factor of the measured correlation
  V2  level + velocity, the note's improvement 1 alone
  C2  level, noise states MA(3) from the window kernel, the note's improvement 2 alone
  VC  both
  CK  C2 plus the same kernel on the signal: z_t is a weighted mean of the last four block gains (the note's "later" item 3)
  T2  tick-level: ell plus the D4 prior as a 4-state disturbance, measurement z_i = exp(ell) Y_i + d_i + e_i, no windows
Writes results/memory_kalman_theory_check.json.
"""
import sys
sys.dont_write_bytecode = True
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import huber_prior_check as H

HERE = Path(__file__).resolve().parent
TC, F0, SIGMA, RHO = .008, 20.8, H.SIGMA, H.RHO
HOP = 6; DT = HOP*TC; NW = 24
FORCE = .037      # probe force in z, measured median
N = 750; K0 = 2000.; BURN = 125
D = np.diff(np.eye(NW), 4, axis=0)
MW = np.eye(NW)-np.linalg.inv(np.eye(NW)+RHO*D.T@D)


def ell_shape(name):
    t = np.arange(N)*TC
    if name == 'plateau':
        return np.zeros(N)
    if name.startswith('step'):
        return float(name[4:])*np.clip((t-3.)/.024, 0., 1.)
    if name == 'ramp':       # A: .27 over .56 s, then flat
        return .27*np.clip((t-3.)/.56, 0., 1.)
    if name == 'dome':       # A_large: round tops, sharp valleys, .7 over .52 s
        return .7*np.abs(np.sin(np.pi*t/1.04))
    raise ValueError(name)


def world(name, seed):
    rng = np.random.default_rng(seed)
    ell = ell_shape(name); k = K0*np.exp(ell); t = np.arange(N)*TC
    Y = FORCE/k*np.sin(2*np.pi*F0*t+rng.uniform(0, 2*np.pi))
    z = np.r_[0., np.cumsum(.5*(k[1:]+k[:-1])*np.diff(Y))]+SIGMA*rng.standard_normal(N)
    return ell, Y, z


def windows(Y, z):
    yw = np.lib.stride_tricks.sliding_window_view(Y, NW); zw = np.lib.stride_tricks.sliding_window_view(z, NW)
    my = yw@MW
    k = np.sum(my*zw, axis=1)/np.sum(my*yw, axis=1)
    out = np.full(N, np.nan); out[NW-1:] = np.log(np.maximum(k, 200.))
    return out


def kernel():
    """Noise and signal weights of one ln-window estimate on its ticks, phase locked since 6 ticks = one probe period."""
    t = np.arange(NW)*TC
    G, S = [], []
    for ph in np.linspace(0, 2*np.pi, 12, endpoint=False):
        Y = FORCE/K0*np.sin(2*np.pi*F0*t+ph); my = MW@Y; info = float(Y@my)
        G.append(my/(info*K0)); S.append(Y*my/info)
    return G, S


def ma_model():
    G, S = kernel()
    c = np.zeros(4)
    for L in range(4):
        c[L] = SIGMA**2*np.mean([g[:NW-HOP*L]@g[HOP*L:] for g in G])
    roots = np.roots(np.r_[c[3], c[2], c[1], c[0], c[1], c[2], c[3]])
    inside = roots[np.abs(roots) < 1.]
    b = np.real(np.poly(inside))[:4]
    su2 = c[0]/float(b@b)
    long_run = float((c[0]+2*c[1:].sum())/c[0])
    # Signal weight of the four 6-tick blocks, newest first.
    s = np.mean(S, axis=0)
    blocks = np.array([s[NW-6*(q+1):NW-6*q].sum() for q in range(4)]); blocks /= blocks.sum()
    return dict(c=c, b=b, su2=su2, long_run=long_run, rho=(c[1:]/c[0]).tolist(), root_mod=np.abs(inside).tolist(), su_over_R=su2/c[0], signal_blocks=blocks)


MA = ma_model()


def imm(modes, h, R, zs, x0, P0, stay=.98, temper=1., gate_h=None):
    """Generic IMM. modes: list of (A, Q). h: observation row (or function of x). Returns the output of the first state per update."""
    n = len(modes); T = np.full((n, n), (1.-stay)/(n-1)); np.fill_diagonal(T, stay)
    xs = [x0.copy() for _ in modes]; Ps = [P0.copy() for _ in modes]; w = np.ones(n)/n
    out = []
    for z in zs:
        joint = w[:, None]*T; pw = joint.sum(axis=0); mix = joint/pw[None, :]
        mx = [sum(mix[i, j]*xs[i] for i in range(n)) for j in range(n)]
        mP = [sum(mix[i, j]*(Ps[i]+np.outer(xs[i]-mx[j], xs[i]-mx[j])) for i in range(n)) for j in range(n)]
        ll = np.zeros(n)
        for j, (A, Q) in enumerate(modes):
            x = A@mx[j]; P = A@mP[j]@A.T+Q
            e = z-h@x; Phx = P@h; S = h@Phx+R
            K = Phx/S
            xs[j] = x+K*e; Ps[j] = P-np.outer(K, Phx)
            ll[j] = -.5*(math.log(2*math.pi*S)+e*e/S)
        lw = np.log(pw)+temper*ll; lw -= lw.max(); w = np.exp(lw); w /= w.sum()
        out.append(sum(w[j]*xs[j][0] for j in range(n)))
    return np.array(out)


def q_velocity(lam, q, h):
    u = np.linspace(0., h, 41)
    if lam > 0:
        a = (1-np.exp(-lam*u))/lam; b = np.exp(-lam*u)
    else:
        a = u; b = np.ones_like(u)
    return q*np.trapz(np.array([[a*a, a*b], [a*b, b*b]]), u, axis=2)


def a_velocity(lam, h):
    if lam > 0:
        return np.array([[1., (1-math.exp(-lam*h))/lam], [0., math.exp(-lam*h)]])
    return np.array([[1., h], [0., 1.]])


def block(A, Q, extra, Qe):
    n, m = A.shape[0], extra.shape[0]
    AA = np.zeros((n+m, n+m)); QQ = np.zeros((n+m, n+m))
    AA[:n, :n] = A; AA[n:, n:] = extra; QQ[:n, :n] = Q; QQ[n:, n:] = Qe
    return AA, QQ


def noise_states():
    E = np.zeros((4, 4)); E[1:, :-1] = np.eye(3)
    Qe = np.zeros((4, 4)); Qe[0, 0] = MA['su2']
    return E, Qe


def f_lv(zs, p):
    R = MA['c'][0]*MA['long_run']
    modes = [(np.eye(1), np.array([[p[0]*DT]])), (np.eye(1), np.array([[p[1]*DT]]))]
    return imm(modes, np.array([1.]), R, zs, np.array([zs[0]]), np.array([[R]]))


def f_current(zs, p):
    # As in rib_transition_eval.make_filter: R = 4 (sig/k)^2 with sig on the nominal noise, gain on R, likelihood on prior + R/4 to the power 1/4.
    r_true = MA['c'][0]; R = 4.*(.796*.04/SIGMA)**2*r_true
    w = np.array([.5, .5]); T = np.array([[.98, .02], [.02, .98]]); q = np.array(p)
    m = np.full(2, zs[0]); V = np.full(2, R); out = []
    for y in zs:
        joint = w[:, None]*T; pw = joint.sum(axis=0); mix = joint/pw[None, :]
        mm = mix.T@m; mp = np.sum(mix*(V[:, None]+(m[:, None]-mm)**2), axis=0)
        prior = mp+q*DT; e = y-mm; iv = prior+R/4.
        lw = np.log(pw)-.5*.25*(np.log(2*np.pi*iv)+e*e/iv); lw -= lw.max(); w = np.exp(lw); w /= w.sum()
        K = prior/(prior+R); m = mm+K*e; V = (1-K)*prior
        out.append(float(w@m))
    return np.array(out)


def f_v2(zs, p):
    R = MA['c'][0]*MA['long_run']
    modes = [(a_velocity(p[0], DT), q_velocity(p[0], p[1], DT)), (a_velocity(0., DT), q_velocity(0., p[2], DT))]
    return imm(modes, np.array([1., 0.]), R, zs, np.array([zs[0], 0.]), np.diag([R, 1.]))


def f_c2(zs, p):
    E, Qe = noise_states(); r = p[2]*MA['c'][0]
    modes = [block(np.eye(1), np.array([[p[0]*DT]]), E, Qe), block(np.eye(1), np.array([[p[1]*DT]]), E, Qe)]
    h = np.r_[1., MA['b']]
    return imm(modes, h, max(r, 1e-9), zs, np.r_[zs[0], np.zeros(4)], np.diag(np.r_[MA['c'][0], np.full(4, MA['su2'])]))


def f_vc(zs, p):
    E, Qe = noise_states(); r = p[3]*MA['c'][0]
    modes = [block(a_velocity(p[0], DT), q_velocity(p[0], p[1], DT), E, Qe), block(a_velocity(0., DT), q_velocity(0., p[2], DT), E, Qe)]
    h = np.r_[1., 0., MA['b']]
    return imm(modes, h, max(r, 1e-9), zs, np.r_[zs[0], 0., np.zeros(4)], np.diag(np.r_[MA['c'][0], 1., np.full(4, MA['su2'])]))


def f_ck(zs, p):
    # Signal states: gain of the last four blocks; noise states as C2. Output is the newest block's gain.
    E, Qe = noise_states(); s = MA['signal_blocks']; r = p[2]*MA['c'][0]
    def mode(q):
        A = np.zeros((4, 4)); A[0, 0] = 1.; A[1:, :-1] = np.eye(3)
        Q = np.zeros((4, 4)); Q[0, 0] = q*DT
        return block(A, Q, E, Qe)
    h = np.r_[s, MA['b']]
    return imm([mode(p[0]), mode(p[1])], h, max(r, 1e-9), zs, np.r_[np.full(4, zs[0]), np.zeros(4)], np.diag(np.r_[np.full(4, MA['c'][0]), np.full(4, MA['su2'])]))


def f_t2(Y, z, ell0, p):
    """EKF per tick; state ell and the D4 disturbance as a companion form, prior innovation variance sigma^2 / rho."""
    A0 = np.zeros((5, 5)); A0[0, 0] = 1.; A0[1, 1:] = [4., -6., 4., -1.]; A0[2:, 1:4] = np.eye(3)
    Qd = SIGMA**2/RHO
    mk = lambda q: np.diag([q*TC, Qd, 0., 0., 0.])
    Qs = [mk(p[0]), mk(p[1])]
    stay = .98**(1./HOP); T = np.array([[stay, 1-stay], [1-stay, stay]])
    xs = [np.r_[ell0, np.zeros(4)] for _ in Qs]; Ps = [np.diag([.25, 1., 1., 1., 1.])*1. for _ in Qs]; w = np.array([.5, .5])
    out = np.full(N, np.nan)
    for i in range(N):
        joint = w[:, None]*T; pw = joint.sum(axis=0); mix = joint/pw[None, :]
        mx = [mix[0, j]*xs[0]+mix[1, j]*xs[1] for j in range(2)]
        mP = [sum(mix[a, j]*(Ps[a]+np.outer(xs[a]-mx[j], xs[a]-mx[j])) for a in range(2)) for j in range(2)]
        ll = np.zeros(2)
        for j in range(2):
            x = A0@mx[j]; P = A0@mP[j]@A0.T+Qs[j]
            kk = math.exp(x[0]+math.log(K0)); h = np.array([kk*Y[i], 1., 0., 0., 0.])
            e = z[i]-(kk*Y[i]+x[1]); Ph = P@h; S = h@Ph+SIGMA**2
            G = Ph/S
            xs[j] = x+G*e; Ps[j] = P-np.outer(G, Ph)
            ll[j] = -.5*(math.log(2*math.pi*S)+e*e/S)
        lw = np.log(pw)+ll; lw -= lw.max(); w = np.exp(lw); w /= w.sum()
        out[i] = w[0]*xs[0][0]+w[1]*xs[1][0]
    return out


GRIDS = dict(
    LV=[(a, b) for a in (.002, .02, .1) for b in (.5, 2., 8., 30.)],
    V2=[(10., a, b) for a in (.01, .1, 1.) for b in (1., 10., 100., 1000.)],
    C2=[(a, b, r) for a in (.002, .02, .1) for b in (.5, 2., 8., 30.) for r in (0., .05)],
    VC=[(10., a, b, r) for a in (.01, .1, 1.) for b in (1., 10., 100., 1000.) for r in (0., .05)],
    CK=[(a, b, r) for a in (.002, .02, .1) for b in (.5, 2., 8., 30.) for r in (0., .05)],
    T2=[(a, b) for a in (.002, .02, .1) for b in (.5, 2., 8., 30.)],
)
FN = dict(LV=f_lv, V2=f_v2, C2=f_c2, VC=f_vc, CK=f_ck)
SHAPES = ('plateau', 'step0.3', 'step0.7', 'ramp', 'dome')
SEEDS = 16


def score(name, ell, est):
    e = est-ell; t = np.arange(N)*TC; ok = np.isfinite(e) & (t >= BURN*TC)
    if name == 'plateau':
        return dict(sd=float(np.sqrt(np.mean(e[ok]**2))))
    if name == 'dome':
        cusp = np.abs((t+.52) % 1.04-.52) < .1
        return dict(rms=float(np.sqrt(np.mean(e[ok]**2))), cusp=float(np.sqrt(np.mean(e[ok & cusp]**2))))
    span = (t >= 3.) & (t < (3.56 if name == 'ramp' else 3.5))
    after = (t >= 3.) & (t < 4.)
    return dict(mae=float(np.mean(np.abs(e[span & ok]))), under=float(max(0., -np.min(e[after & ok]))))


def run(task):
    fam, p = task
    acc = {}
    for s in SHAPES:
        rows = []
        for seed in range(SEEDS):
            ell, Y, z = world(s, [seed, SHAPES.index(s)])
            if fam == 'WIN':
                est = windows(Y, z)
            elif fam == 'T2':
                est = f_t2(Y, z, 0., p)
            else:
                w = windows(Y, z)
                idx = np.arange(NW-1+HOP, N, HOP)
                out = (f_current if fam == 'CUR' else FN[fam])(w[idx], p)
                est = np.full(N, np.nan)
                est[idx] = out
                for i in range(idx[0], N):     # zero-order hold between updates
                    if not np.isfinite(est[i]):
                        est[i] = est[i-1]
                est = est-math.log(K0)
            if fam == 'WIN':
                est = est-math.log(K0)
            if fam == 'T2':
                pass
            rows.append(score(s, ell, est))
        acc[s] = {k: float(np.mean([r[k] for r in rows])) for k in rows[0]}
    return dict(fam=fam, p=p, **acc)


def real_windows():
    """Part 1: measured against predicted correlation on recorded windows, and the size and memory of the non-sensor error."""
    import memory_scene_facts as F
    pred, res = [], [F.case(p) for p in sorted(H.WIN.glob('case_*.npz'))]
    S = H.Solver()
    for path in sorted(H.WIN.glob('case_*.npz'))[:12]:
        data = H.load(path)
        for j in range(100, len(data['Y'])-20, 25):
            row = []
            base = []
            for L in (0, 6, 12, 18):
                Y = H.tail(data['Y'][j+L], 24)
                if len(Y) != 24:
                    row = None; break
                a = S.fit(Y, H.tail(data['z'][j+L], 24), iters=1)
                if not a['k'] > 0:
                    row = None; break
                base.append(a['my']/(a['info']*a['k']))
            if row is None:
                continue
            v0 = base[0]@base[0]
            pred.append([base[0][:24-L]@base[q][L:]/v0 for q, L in ((1, 6), (2, 12), (3, 18))])
    lags = [6, 12, 18, 24]
    out = dict(pred_rho=np.median(pred, axis=0).tolist(), meas_rho=F.acf([r['e_n'] for r in res], lags), e0_rho=F.acf([r['e_0']-np.nanmedian(r['e_0']) for r in res], lags))
    for g in ('A', 'A_large'):
        rr = [r for r in res if r['meta']['group'] == g]
        out[g] = dict(e_n_rms=float(np.sqrt(np.nanmean(np.concatenate([r['e_n'] for r in rr])**2))), e0_rms=float(np.sqrt(np.nanmean(np.concatenate([r['e_0']-np.nanmedian(r['e_0']) for r in rr])**2))), e0_rho=F.acf([r['e_0']-np.nanmedian(r['e_0']) for r in rr], lags))
    return out


def main():
    print('MA model', {k: (np.round(v, 4).tolist() if not isinstance(v, float) else round(v, 4)) for k, v in MA.items()})
    I = (FORCE/SIGMA)**2/2*.77
    print('tick info', round(I, 3))
    print('slope accuracy: endpoint sd of level and sd of velocity by linear fit over T; level-only sd over T')
    for T in (.1, .2, .3, .5, .8):
        n = T/TC
        print(f'  T {T:.1f} s: level-only {1/math.sqrt(I*n):.3f}, linear-fit endpoint {2/math.sqrt(I*n):.3f}, velocity sd {math.sqrt(12/(I*n))/T:.2f} /s')
    print('reversal +v to -v: time for the gap 2 v t to reach log-likelihood 7 (10 s between false alarms), and the gap then')
    for v in (1., 2.4, 5., 8.):
        t = (7.*6.*TC/(I*4*v*v))**(1/3)
        print(f'  v {v}: {t:.3f} s, gap {2*v*t:.2f}')
    rw = real_windows()
    print('real windows', json.dumps(rw, indent=0))
    tasks = [('WIN', ()), ('CUR', (.0, .0))]
    from memory_filters2 import Q_RATE_HI, Q_RATE_LO
    tasks[1] = ('CUR', (Q_RATE_LO, Q_RATE_HI))
    for fam, grid in GRIDS.items():
        tasks += [(fam, p) for p in grid]
    with ProcessPoolExecutor(8) as ex:
        rows = list(ex.map(run, tasks, chunksize=2))
    (HERE/'results/memory_kalman_theory_check.json').write_text(json.dumps(dict(ma={k: np.asarray(v).tolist() for k, v in MA.items()}, tick_info=I, real_windows=rw, rows=rows), indent=1), encoding='utf-8')
    fmt = lambda r: f"sd {r['plateau']['sd']:.3f} | step.3 {r['step0.3']['mae']:.3f}/{r['step0.3']['under']:.2f} step.7 {r['step0.7']['mae']:.3f}/{r['step0.7']['under']:.2f} ramp {r['ramp']['mae']:.3f}/{r['ramp']['under']:.2f} dome {r['dome']['rms']:.3f} cusp {r['dome']['cusp']:.3f}"
    print('raw window  ', fmt(next(r for r in rows if r['fam'] == 'WIN')))
    print('current IMM ', fmt(next(r for r in rows if r['fam'] == 'CUR')))
    for target in (.06, .08, .10, .13):
        print(f'--- best dome RMS with plateau sd <= {target}')
        for fam in GRIDS:
            q = [r for r in rows if r['fam'] == fam and r['plateau']['sd'] <= target]
            if q:
                b = min(q, key=lambda r: r['dome']['rms'])
                print(f'  {fam} {fmt(b)}  p={b["p"]}')
            else:
                print(f'  {fam} none; lowest sd', round(min(r['plateau']['sd'] for r in rows if r['fam'] == fam), 3))


if __name__ == '__main__':
    main()
