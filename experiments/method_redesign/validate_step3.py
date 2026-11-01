"""Targeted checks: prediction chronology, QP, inherited modules, robustness map.
No final cohort or tuning/evaluation scene outcomes are used here.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
from types import SimpleNamespace
import numpy as np
from mpc_controller import *
from mpc_robustness import base_map, screen, kgain, AT, Q1, FRAC, SHIFT


def main():
    rng = np.random.default_rng(40901)
    # Independent 1 ms rollout of the prediction model; catches sample/servo
    # order, one-period command delay, CA sign and endpoint raw-force constraint.
    for H in (10, 15, 25):
        BF, BR, GF, GR = prediction_matrices(H)
        st = rng.normal(size=5); z = rng.normal(size=H)
        u, ff, f0, v, a = st
        y = 0.; fps, rps = [], []
        for i in range(H):
            for m in range(NS):
                dt = i*TC+m*.001
                raw = f0+y+v*dt+.5*a*dt*dt
                ff += AF*(raw-ff)
                y = np.exp(-.001/TAU)*y+(1-np.exp(-.001/TAU))*u
            u += z[i]
            dt = (i+1)*TC
            fps.append(ff); rps.append(f0+y+v*dt+.5*a*dt*dt)
        np.testing.assert_allclose(BF@st+GF@z, fps, atol=1e-12)
        np.testing.assert_allclose(BR@st+GR@z, rps, atol=1e-12)
        np.testing.assert_allclose(GF[0], 0., atol=0.)
        np.testing.assert_allclose(GR[0], 0., atol=0.)
    print('PASS condensed prediction versus independent substep rollout')
    u = ForceMPC(H=15, qr=.1, terminal=2., constrained=False)
    c = ForceMPC(H=15, qr=.1, terminal=2., constrained=True)
    st = np.array([0., 4., 4., 0., 0.])
    assert abs(c.solve(st, 0., 1e4, 1e4)-u.solve(st, 0., 1e4, 1e4)) < 1e-8
    st[1:3] = 5.
    assert abs(c.solve(st, 0., 1e4, 1e4)-u.solve(st, 0., 1e4, 1e4)) < 1e-8
    c.solve(np.array([0., 30., 30., 0., 0.]), 0., 1e4, 1e4)
    assert c.infeasible == 1
    print('PASS unconstrained/QP agreement with inactive bounds; delayed infeasibility is explicit')
    # Replay the SAME observation buffers into new and original policies:
    # estimation, excitation state, IMM and recovery must match exactly despite
    # different output commands. Synthetic forcing includes lost/recontact ticks.
    old = CE.build('Proposed')
    new = EstimationAwareMPC()
    pl = common.P.PLATFORMS['PC125']
    old.reset(pl); new.reset(pl)
    n = 180; log = np.zeros((n, 16))
    for j in range(n):
        t = j*TC
        log[j, 0] = t
        log[j, 2] = .001+1e-5*np.sin(2*np.pi*20.8*t)
        log[j, 11] = log[j, 2]+1e-5
        log[j, 5] = 4.+.1*np.sin(2*np.pi*20.8*t)
        if 95 <= j < 110: log[j, 5] = .1
        log[j, 8] = np.nan  # neither deployment policy may read true stiffness
        log[j, 12:] = np.nan
        ob = SimpleNamespace(n=j+1, log=log)
        old_command = old.step(ob); new_command = new.step(ob)
        if not new._controlled:
            np.testing.assert_equal(old_command, new_command, err_msg='recovery command')
        for key in ('kbar', 'k_good', 'k_hat', 'first', 'mode', 'k_seek', 'x_prev', 'Z'):
            np.testing.assert_equal(getattr(old, key), getattr(new, key), err_msg=key)
        np.testing.assert_array_equal(old.tracker.x, new.tracker.x)
        np.testing.assert_array_equal(old.tracker.mu, new.tracker.mu)
    print('PASS exact estimation/IMM/dither/supervisor replay; deployment truth columns poisoned')
    # Check DARE gains against the original iterative covariance formula.
    for q in (CE.FINAL['q_smooth'], CE.FINAL['q_maneuver']):
        kr = 1e4; R = CE.FINAL['r_scale']*(.04/kr)**2
        P = np.eye(3)*1e-6
        for _ in range(10000):
            pp = AT@P@AT.T+q*Q1
            K = pp[:, 0]/(pp[0, 0]+R)
            P = (np.eye(3)-np.outer(K, [1, 0, 0]))@pp
        np.testing.assert_allclose(kgain(q, kr), K, rtol=1e-8)
    print('PASS steady-state gains against iterative FINAL-R covariance update')
    # Direct one-tick construction versus condensed screening map for a frozen
    # mode, with offset-free predictor state and mismatched stiffness channels.
    kap, kg, k, mu = 1.4, 1.09, 3e4, .5
    A, S, B = base_map(kap, kg, k, mu)
    par = dict(H=15, qr=.1, terminal=1.)
    m = ForceMPC(**par, constrained=False)
    r = kap/kg
    for _ in range(10):
        x = rng.normal(size=12)*1e-3
        du = m.solve(S@x, B@x, r, 1.)-m.solve(np.zeros(5), 0., r, 1.)
        nx = A@x; nx[1] += du; nx[11] = m.BF[0]@(S@x)
        gain = np.linalg.solve(m.GTG+np.eye(m.H)/(m.qr*r*r), m.GF.T*m.weights)[0]
        am = A.copy(); am[1] -= gain@(m.BF@S+np.ones((m.H, 1))*B)/r; am[11] = m.BF[0]@S
        np.testing.assert_allclose(nx, am@x, atol=1e-12)
    print('PASS unconstrained solver equals closed-loop screening map')
    good = screen(par); bad = screen(dict(H=15, qr=100., terminal=1.))
    assert good['accepted'] and not bad['accepted']
    assert good['unit_count_range'] == [0, 0]
    print('PASS full 525-point screen separates stable/aggressive configurations')
    print('ALL CHECKS PASSED')


if __name__ == '__main__': main()
