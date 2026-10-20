"""Nonlinear sliding contact, four actual gains and a common energy filter.

An explicitly reduced, known-model test, NOT the cylindrical Hertz/MHE/robot B0.
The contact law is a unilateral 3/2-power spring; no image or clinical claims.
"""
import numpy as np
from numba import njit

HORIZON = .30
DT = .001
STEPS = 300
FD = 3.
SPEED = .008
MASS = np.array([.8, .5])
ROBOT_DAMPING = .5
SCALE = np.array([500., 500., 40., 40.])
LOW = np.array([5., 5., 2., 2.]) / SCALE
HIGH = np.array([1600., 1600., 120., 120.]) / SCALE
RATE = np.array([800., 800., 80., 80.]) / SCALE
E_MIN = 1e-5
CHI = .5
ALPHA = 40.


def sample(n, seed):
    """State is position/velocity errors, normalized actual gains, tank energy.

    Old/new local slopes are revised at t=0, not advance knowledge of tissue.
    Both models and exact states are shared by all methods (an oracle layer).
    """
    rng = np.random.default_rng(seed)
    k = np.exp(rng.uniform(np.log(300.), np.log(1500.), n))
    old_k = k * rng.uniform(.75, 1.25, n)
    friction = rng.uniform(.1, .6, n)
    state = np.zeros((n, 9))
    state[:, 0] = rng.uniform(-.002, .002, n)
    state[:, 1] = 1.5 * FD / old_k - 1.5 * FD / k
    state[:, 2:4] = rng.uniform(-.02, .02, (n, 2))
    state[:, 4:8] = np.exp(rng.uniform(np.log([40., 20., 6., 6.]),
                                     np.log([600., 500., 50., 50.]), (n, 4))) / SCALE
    state[:, 8] = E_MIN + np.exp(rng.uniform(np.log(5e-5), np.log(.005), n))
    return np.column_stack((k, friction)), state


@njit(cache=True)
def contact(y, p):
    penetration = y[1] + 1.5 * FD / p[0]
    fn = FD * (max(penetration, 0.) / (1.5 * FD / p[0])) ** 1.5
    ft = p[1] * fn * np.tanh((y[2] + SPEED) / .002)
    return fn, ft


@njit(cache=True)
def action(y, target):
    g = y[4:8]
    raw = (target - g) / .025
    lo = np.maximum(-RATE, -ALPHA * (g - LOW))
    hi = np.minimum(RATE, ALPHA * (HIGH - g))
    u = np.minimum(np.maximum(raw, lo), hi)
    # Projection in normalized gain-rate coordinates onto box + one halfspace.
    a = np.array([.5*y[0]**2*SCALE[0], .5*y[1]**2*SCALE[1], 0., 0.])
    pd = g[2]*SCALE[2]*y[2]**2 + g[3]*SCALE[3]*y[3]**2
    pd += ROBOT_DAMPING*((y[2]+SPEED)**2 + y[3]**2)
    budget = CHI*pd + ALPHA*(y[8]-E_MIN)
    active = np.dot(a, u) > budget
    if active:
        left, right = 0., 1.
        for _ in range(60):
            candidate = np.minimum(np.maximum(raw-right*a, lo), hi)
            if np.dot(a, candidate) <= budget:
                break
            right *= 2.
        for _ in range(40):
            middle = (left+right)/2
            candidate = np.minimum(np.maximum(raw-middle*a, lo), hi)
            if np.dot(a, candidate) > budget:
                left = middle
            else:
                right = middle
        u = np.minimum(np.maximum(raw-right*a, lo), hi)
    return u, CHI*pd-np.dot(a, u), active


@njit(cache=True)
def flow(y, target, p):
    fn, ft = contact(y, p)
    u, edot, active = action(y, target)
    g = y[4:8]*SCALE
    ut_raw = p[1]*FD - g[0]*y[0] - g[2]*y[2]
    un_raw = FD - g[1]*y[1] - g[3]*y[3]
    ut = min(max(ut_raw, -12.), 12.)
    un = min(max(un_raw, -12.), 12.)
    dy = np.empty(9)
    dy[:2] = y[2:4]
    dy[2] = (ut-ft-ROBOT_DAMPING*(y[2]+SPEED))/MASS[0]
    dy[3] = (un-fn-ROBOT_DAMPING*y[3])/MASS[1]
    dy[4:8] = u
    dy[8] = edot
    # Tracking, contact-force error and executed actuator variation; no gain goal.
    cost = (y[0]/.002)**2 + (y[1]/.006)**2
    cost += (y[2]/.03)**2 + (y[3]/.03)**2 + ((fn-FD)/.7)**2
    cost += ((ut-p[1]*FD)/4.)**2 + ((un-FD)/4.)**2
    return dy, cost, fn, ut_raw, un_raw, active


@njit(cache=True)
def rollout(targets, state, p, delay=0., dt=DT):
    y = state.copy()
    total = 0.
    diagnostics = np.array([1e10, 0., 0., 0., 1e10, 0., 0., 0.])
    # Integrate exact release/block times, including sub-millisecond inference.
    bounds=np.array([0.,min(delay,HORIZON),min(delay+HORIZON/2,HORIZON),HORIZON])
    for phase in range(3):
        duration=bounds[phase+1]-bounds[phase]
        count=int(np.ceil(duration/dt-1e-12))
        target=state[4:8] if phase==0 else targets[phase-1]
        for i in range(count):
            step=min(dt,duration-i*dt)
            a, ja, fn, ut, un, active = flow(y, target, p)
            b, jb, _, _, _, _ = flow(y+step*a/2, target, p)
            c, jc, _, _, _, _ = flow(y+step*b/2, target, p)
            d, jd, _, _, _, _ = flow(y+step*c, target, p)
            diagnostics[0] = min(diagnostics[0], fn)
            diagnostics[1] = max(diagnostics[1], fn)
            diagnostics[2] += (1. if fn < .2 else 0.)*step
            diagnostics[3] += (1. if fn > 6. else 0.)*step
            diagnostics[4] = min(diagnostics[4], y[8]-E_MIN)
            diagnostics[5] += (1. if active else 0.)*step
            diagnostics[6] = max(diagnostics[6], max(abs(ut), abs(un)))
            diagnostics[7] = max(diagnostics[7], max(np.max(LOW-y[4:8]), np.max(y[4:8]-HIGH)))
            y += step*(a+2*b+2*c+d)/6
            total += step*(ja+2*jb+2*jc+jd)/6
    diagnostics[4] = min(diagnostics[4], y[8]-E_MIN)
    return total, diagnostics, y


@njit(cache=True)
def batch_cost(targets, states, params):
    out = np.empty(len(states))
    for i in range(len(states)):
        out[i] = rollout(targets[i], states[i], params[i])[0]
    return out


@njit(cache=True)
def candidate_cost(targets, state, param):
    out = np.empty(len(targets))
    for i in range(len(targets)):
        out[i] = rollout(targets[i], state, param)[0]
    return out


def normalized_targets(z):
    return LOW + np.asarray(z).reshape(-1, 4)*(HIGH-LOW)


def target_coordinates(g):
    return (np.asarray(g)-LOW)/(HIGH-LOW)


def features(params, states):
    return np.column_stack((np.log(params[:, 0]/600.), params[:, 1],
        states[:, :2]/.003, states[:, 2:4]/.02,
        np.log(states[:, 4:8]), np.log((states[:, 8]-E_MIN)/.001)))
