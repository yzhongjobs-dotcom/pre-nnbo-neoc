"""Offline design of the unit-stiffness admittance A1 (paper draft §三.4).

Force-domain controller with unit stiffness: M1 w'' + D1 w' = e + ki I, I' = e, and the
command increment Δy = (Δw + K_PF1 Δe)/k_bar. With k = κ k_bar the loop is the same
for every stiffness up to the gain ratio κ. The design picks (M1, D1, K_PF1) on a grid
(ki fixed) to minimise a nominal force-error cost plus a noise-to-command term, subject
to frozen-gain stability for every κ in [a, b(1+δ)] and the circle-criterion disk test
(Nyquist of L1 outside the disk spanned by -1/a and -1/(b(1+δ))). Platform: 8 ms
period, one-period command delay, 19 ms first-order servo, 30 Hz force filter (sampled
at the control period here; the plant filters at 1 kHz).
Run: python -X utf8 design_A1.py  -> writes A1_design.json
"""
import json, math, os
import numpy as np

Tc, tau_s, fc = 0.008, 0.019, 30.0
a_s = math.exp(-Tc/tau_s)
a_f = math.exp(-2*math.pi*fc*Tc)
KI = 0.0     # with ki > 0 the loop is type 2 and no design passes the circle test (see main)
SECTOR = (0.5, 2.0)          # design sector [a, b]
DELTA = 0.07                 # commutation-error allowance on the upper end (rho_d = 0.02)
HERE = os.path.dirname(os.path.abspath(__file__))


BETA = 1-math.exp(-1e-3/tau_s)              # servo per 1 ms sub-step (plant)
ALPHA = 1-math.exp(-2*math.pi*fc*1e-3)      # force filter per 1 ms sub-step (plant)
NSUB = int(round(Tc/1e-3))
NX = 7   # y, cmd_prev, Ff, I, wd, cmd, m_prev


def _tick(x, s, m, M1, D1, KPF1, kappa):
    """one control period of the linearised plant + force-domain controller, unit k_bar.
    m: force measurement fed to the controller (closed loop: m = Ff)."""
    y, cp, Ff, I, wd, cmd, mp = x
    e, ep = -m, -mp
    eD = math.exp(-D1*Tc/M1)
    I2 = I+e*Tc
    wd2 = wd*eD+(1-eD)*(e+KI*I2)/D1
    cp2 = cmd
    cmd2 = cmd+0.5*Tc*(wd+wd2)+KPF1*(e-ep)
    for _ in range(NSUB):
        Ff = Ff+ALPHA*(kappa*(y-s)-Ff)
        y = y+BETA*(cp2-y)
    return np.array([y, cp2, Ff, I2, wd2, cmd2, m])


def loop_matrices(M1, D1, KPF1, kappa):
    """closed loop x+ = A x + Bs s (+ Bd d, force noise added to the measurement)."""
    A = np.zeros((NX, NX))
    for i in range(NX):
        x = np.zeros(NX); x[i] = 1.
        A[:, i] = _tick(x, 0., x[2], M1, D1, KPF1, kappa)
    Bs = _tick(np.zeros(NX), 1., 0., M1, D1, KPF1, kappa)
    Bd = _tick(np.zeros(NX), 0., 1., M1, D1, KPF1, kappa)
    return A, Bs, Bd


def _open_loop(M1, D1, KPF1):
    """open loop at the measurement: x+ = Ao x + Bo m, output Ff = Co x (unit kappa)."""
    Ao = np.zeros((NX, NX))
    for i in range(NX):
        x = np.zeros(NX); x[i] = 1.
        Ao[:, i] = _tick(x, 0., 0., M1, D1, KPF1, 1.0)
    Bo = _tick(np.zeros(NX), 0., 1., M1, D1, KPF1, 1.0)
    Co = np.zeros(NX); Co[2] = 1.
    return Ao, Bo, Co


def open_loop_freq(M1, D1, KPF1, w):
    """L1(e^{jwTc}) with negative feedback convention (closed loop: m = Ff)."""
    Ao, Bo, Co = _open_loop(M1, D1, KPF1)
    out = []
    for wi in w:
        z = np.exp(1j*wi*Tc)
        out.append(-(Co@np.linalg.solve(z*np.eye(NX)-Ao, Bo)))
    return np.array(out)


def stable(M1, D1, KPF1, kappas):
    for k in kappas:
        A, _, _ = loop_matrices(M1, D1, KPF1, k)
        if KI == 0:                      # integral state unused: drop its unit eigenvalue
            keep = [i for i in range(NX) if i != 3]
            A = A[np.ix_(keep, keep)]
        if max(abs(np.linalg.eigvals(A))) >= 1-1e-9:
            return False
    return True


def circle_ok(M1, D1, KPF1, a, b):
    w = np.logspace(-1, math.log10(math.pi/Tc*0.999), 1500)
    L = open_loop_freq(M1, D1, KPF1, w)
    c = -(1/a+1/b)/2; r = (1/a-1/b)/2
    return bool(np.min(np.abs(L-c)) > r*1.02), float(np.min(np.abs(L-c))/r)


def cost(M1, D1, KPF1, lam_n=0.02, T=4.0):
    """nominal (kappa=1) cost: force error (unit stiffness, in mm-equivalent) for a 0.5 mm
    surface step plus 2.5 mm, 0.25 Hz breathing; plus noise-to-command gain."""
    A, Bs, Bd = loop_matrices(M1, D1, KPF1, 1.0)
    n = int(T/Tc); x = np.zeros(NX); err = 0.
    for i in range(n):
        t = i*Tc
        s = 0.5e-3*(t > 0.5)+2.5e-3*math.sin(2*math.pi*0.25*t)
        x = A@x+Bs*s
        f = x[0]-s
        err += f*f
    rms = math.sqrt(err/n)*1e3
    # noise: H2 gain from force noise to command increment
    x = Bd.copy(); g = 0.; xp = np.zeros(NX)
    for i in range(400):
        g += (x[5]-xp[5])**2; xp = x.copy(); x = A@x
    return rms+lam_n*math.sqrt(g)*1e3, rms, math.sqrt(g)


def design():
    a, b = SECTOR
    kap = np.concatenate([np.linspace(a, 1, 6), np.linspace(1, b*(1+DELTA), 8)])
    best = None; rows = []
    for M1 in np.logspace(-6, -4, 11):
        for D1 in np.logspace(-3, -1, 13):
            for KPF1 in np.linspace(0, 3, 7):
                if not stable(M1, D1, KPF1, kap):
                    continue
                ok, marg = circle_ok(M1, D1, KPF1, a, b*(1+DELTA))
                if not ok:
                    continue
                J, rms, ng = cost(M1, D1, KPF1)
                rows.append((J, M1, D1, KPF1, rms, ng, marg))
                if best is None or J < best[0]:
                    best = (J, M1, D1, KPF1, rms, ng, marg)
    return best, rows


def feasible_set():
    """the circle-screened design set actually used for tuning (A1_feasible.json)."""
    a, b = SECTOR
    kap = np.concatenate([np.linspace(a, 1, 4), np.linspace(1, b*(1+DELTA), 6)])
    feas = []
    for M1 in np.logspace(-6, -1, 16):
        for D1 in np.logspace(-3, 0, 25):
            for K in (0.0, 0.5, 1.0):
                if not stable(M1, D1, K, kap):
                    continue
                ok, m = circle_ok(M1, D1, K, a, b*(1+DELTA))
                if ok:
                    feas.append(dict(M1=float(M1), D1=float(D1), KPF1=float(K), ki=KI, margin=m, bw=1/D1))
    feas.sort(key=lambda x: -x['bw'])
    return feas


if __name__ == '__main__':
    # 1) type-2 check: with ki = 30 no design passes the disk test
    global_ki = KI
    KI = 30.0
    n30 = sum(circle_ok(M1, D1, K, SECTOR[0], SECTOR[1]*(1+DELTA))[0]
              for M1 in np.logspace(-6, -1, 6) for D1 in np.logspace(-3, 0, 8) for K in (0.0, 1.0))
    KI = global_ki
    feas = feasible_set()
    json.dump(feas, open(os.path.join(HERE, 'A1_feasible.json'), 'w'), indent=1)
    print(json.dumps(dict(n_feasible_ki0=len(feas), n_feasible_ki30_coarse=int(n30),
                          fastest=feas[:3]), indent=1))
