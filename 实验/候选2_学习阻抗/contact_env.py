"""Reduced Q2 carrier: causal observations, actual gains and shared tank shield.

True material parameters are used only by the plant and cost evaluation.
The surface-motion program is an exogenous, commonly available task input.
"""
from pathlib import Path
import sys
import numpy as np
from numba import njit
from gymnasium.spaces import Box
from stable_baselines3.common.vec_env import VecEnv

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent / '动态接触强对照' / 'Q2'))
from contact_stages import Problem, environment, margins
from projected_rival import step_projected, action

BASE = Problem()
LO, HI = BASE.lo, BASE.hi
H, DT, SUB, N = .125, .005, 25, 32


def cases(seed, n, kind='train'):
    r = np.random.default_rng(seed)
    p = np.tile(BASE.p, (n, 1))
    p[:, 4] = r.uniform(200, 400, n)
    p[:, 5] = r.uniform(100, 400, n)
    p[:, 1] = r.uniform(.8, 1.2, n)
    pulses = np.tile(BASE.pulses, (n, 1, 1))
    pulses[:, :, 0] += r.uniform(-.1, .1, (n, 2))
    pulses[:, :, 1] *= r.uniform(.9, 1.1, (n, 2))
    pulses[:, :, 2] *= r.uniform(.8, 1.2, (n, 2))
    noise = np.tile([1e-5, 1e-4, .005], (n, 1))
    if kind == 'material_ood':
        p[:, 4] = r.uniform(450, 600, n)
        p[:, 5] = r.uniform(450, 650, n)
    elif kind == 'noise_mass_ood':
        p[:, 0] = r.uniform(.24, .28, n)
        noise *= 3
    return p, pulses, noise


@njit(cache=True)
def advance(states, actions, index, ps, pulses, sub=SUB):
    """Continuous shield at every RK stage, held desired rate each decision."""
    batch = len(states)
    costs = np.zeros(batch)
    metrics = np.zeros((batch, 5))  # min margin, peak Fn, integral Fe², integral e², intervention
    metrics[:, 0] = 1e99
    dt = H / sub
    for b in range(batch):
        y = states[b].copy()
        target = np.exp(np.log(LO) + (.5 + .5 * np.minimum(1., np.maximum(-1., actions[b]))) * np.log(HI / LO))
        desired = np.minimum(1., np.maximum(-1., (target - y[2:4]) / (H * ps[b, 6:8])))
        for s in range(sub):
            t = index[b] * H + s * dt
            yn, ell, ya, yb, yc, u, ua, ub, uc = step_projected(y, desired, t, dt, ps[b], pulses[b], LO, HI)
            costs[b] += ell
            un = action(yn, desired, ps[b], LO, HI)
            for state, rate, phase in ((y, u, t), (ya, ua, t+dt/2), (yb, ub, t+dt/2), (yc, uc, t+dt), (yn, un, t+dt)):
                mar = margins(state, rate, phase, ps[b], pulses[b])
                metrics[b, 0] = min(metrics[b, 0], np.min(mar))
            th, gam = environment(t, ps[b], pulses[b])
            fe = th * (gam-y[0])
            metrics[b, 1] = max(metrics[b, 1], ps[b, 2]+fe)
            metrics[b, 2] += fe*fe*dt
            metrics[b, 3] += y[0]*y[0]*dt
            if np.max(np.abs(u-desired)) > 1e-8:
                metrics[b, 4] += 1. / sub
            y = yn
        if index[b] == N-1:
            th, gam = environment(4., ps[b], pulses[b])
            costs[b] += ps[b, 29]*((y[0]/ps[b, 14])**2+ps[b, 20]*(th*(gam-y[0])/ps[b, 16])**2)
        states[b] = y
    return costs, metrics


def action_from_target(target):
    return np.clip(2*(np.log(np.clip(target, LO, HI)/LO)/np.log(HI/LO)) - 1, -1, 1)


class ContactVec(VecEnv):
    def __init__(self, n=64, history=4, seed=0, fixed=None, sub=SUB):
        self.rng = np.random.default_rng(seed)
        self.history = history
        self.fixed = fixed
        self.sub = sub
        self.states = np.zeros((n, 5))
        self.index = np.zeros(n, dtype=np.int64)
        self.ps, self.pulses, self.noise = cases(seed, n)
        if fixed is not None:
            self.ps, self.pulses, self.noise = [a.copy() for a in fixed]
        self.frames = np.zeros((n, history, 6))
        self.sums = np.zeros((n, 6))
        self.sums[:, 1] = 1e99
        self.sums[:, 2] = 3.
        self.actions = None
        self.last_metrics = None
        self.render_mode = None
        super().__init__(n, Box(-np.inf, np.inf, (6*history+7,), dtype=np.float32), Box(-1., 1., (2,), dtype=np.float32))

    def _frame(self):
        y = self.states
        th_gam = np.array([environment(self.index[i]*H, self.ps[i], self.pulses[i]) for i in range(self.num_envs)])
        measured = np.column_stack([y[:, 0], y[:, 1], th_gam[:, 0]*(th_gam[:, 1]-y[:, 0])])
        measured += self.rng.normal(size=measured.shape)*self.noise
        return np.column_stack([measured[:, 0]/.003, measured[:, 1]/.03, measured[:, 2],
                                (y[:, 2]-50)/400, (y[:, 3]-3)/12, (y[:, 4]-.001)/.0001])

    def _obs(self):
        program = self.pulses / np.array([4., 1., .003])
        return np.column_stack([self.frames.reshape(self.num_envs, -1), self.index/N,
                                program.reshape(self.num_envs, -1)]).astype(np.float32)

    def reset(self):
        self.states[:] = BASE.y0
        self.index[:] = 0
        self.sums[:] = 0
        self.sums[:, 1] = 1e99
        self.sums[:, 2] = 3.
        self.frames[:] = self._frame()[:, None, :]
        return self._obs()

    def step_async(self, actions):
        self.actions = np.asarray(actions, dtype=float)

    def step_wait(self):
        cost, met = advance(self.states, self.actions, self.index, self.ps, self.pulses, self.sub)
        self.index += 1
        self.sums[:, 0] += cost
        self.sums[:, 1] = np.minimum(self.sums[:, 1], met[:, 0])
        self.sums[:, 2] = np.maximum(self.sums[:, 2], met[:, 1])
        self.sums[:, 3:] += met[:, 2:]
        done = self.index == N
        self.frames[:, :-1] = self.frames[:, 1:]
        self.frames[:, -1] = self._frame()
        obs = self._obs()
        # Tank/gain constraints have the shared shield; force/actuator violations
        # receive a disclosed training penalty, and remain evaluation outcomes.
        reward = -cost - 10*np.maximum(-met[:, 0], 0)
        info = [{} for _ in range(self.num_envs)]
        if done.any():
            self.last_metrics = self.sums.copy()
            for i in np.flatnonzero(done):
                info[i] = dict(terminal_observation=obs[i].copy(), episode=dict(r=-self.sums[i, 0], l=N), metrics=self.sums[i].tolist())
                if self.fixed is None:
                    p, pu, no = cases(int(self.rng.integers(2**31)), 1)
                    self.ps[i], self.pulses[i], self.noise[i] = p[0], pu[0], no[0]
                self.states[i] = BASE.y0
                self.index[i] = 0
                self.sums[i] = [0., 1e99, 3., 0., 0., 0.]
            fresh = self._frame()
            self.frames[done] = fresh[done, None, :]
            obs = self._obs()
        return obs, reward.astype(np.float32), done, info

    def close(self):
        pass

    def get_attr(self, attr_name, indices=None):
        return [getattr(self, attr_name)]*len(self._get_indices(indices))

    def set_attr(self, attr_name, value, indices=None):
        setattr(self, attr_name, value)

    def env_method(self, method_name, *args, indices=None, **kwargs):
        return [getattr(self, method_name)(*args, **kwargs)]

    def env_is_wrapped(self, wrapper_class, indices=None):
        return [False]*len(self._get_indices(indices))


def baseline(obs, parameters, kind, history):
    frame = obs[:, 6*(history-1):6*history]
    e, v, fe = frame[:, 0]*.003, frame[:, 1]*.03, frame[:, 2]
    if kind == 'fixed':
        return np.tile(action_from_target(np.asarray(parameters)), (len(obs), 1))
    if kind == 'feedback':
        ak, ad, le, lv, bias = parameters
        signal = (e/np.exp(le))**2+(v/np.exp(lv))**2+(fe/.5)**2
        gate = signal/(1+signal) * (.5+.5*np.tanh(bias*e*v/np.exp(le+lv)))
        target = np.column_stack([50+ak*gate, 3+ad*gate])
    elif kind in ('preview', 'prior_preview'):
        zeta, lead, mult = parameters
        # Ridge regression from the same causal observation window, not theta truth.
        t = obs[:, -7]*4
        programs = obs[:, -6:].reshape(-1, 2, 3)*np.array([4., 1., .003])
        num = np.zeros(len(obs)); den = np.zeros(len(obs))
        for j in range(history):
            tj = np.maximum(0., t-(history-1-j)*H)
            gam = surface(tj, programs)
            x = gam-obs[:, 6*j]*.003
            f = obs[:, 6*j+2]
            num += x*f
            den += x*x
        ridge = (.005/150)**2
        theta = np.clip((num+300*ridge)/(den+ridge), 100, 1500)
        if kind == 'prior_preview':
            theta[:] = 300.  # Known-program schedule, no measured material update.
        # Quasi-static physical cost minimizer; preview uses known surface program.
        targetk = np.full(len(obs), 20.)
        for delta in np.linspace(0., lead, 9):
            gam = surface(t+delta, programs)
            low = np.full(len(obs), 20.); high = np.full(len(obs), 1200.)
            for _ in range(24):
                k = (low+high)/2
                derivative = (theta*gam)**2*(-1e6+1.2*theta*k)/(theta+k)**3+.3*(k-50)/400**2
                low = np.where(derivative<0, k, low)
                high = np.where(derivative>=0, k, high)
            targetk = np.maximum(targetk, mult*(low+high)/2-400*delta)
        target = np.column_stack([targetk, 2*zeta*np.sqrt(.2*(theta+targetk))-1])
    else:
        raise ValueError(kind)
    return action_from_target(target)


def surface(t, programs):
    phase = (t[:, None]-programs[:, :, 0])/programs[:, :, 1]
    return np.sum(np.where((phase>0)&(phase<1), programs[:, :, 2]*np.sin(np.pi*phase)**2, 0), axis=1)


def evaluate(controller, fixed, history=4, seed=97531, sub=SUB):
    env = ContactVec(len(fixed[0]), history, seed, fixed=fixed, sub=sub)
    obs = env.reset()
    for _ in range(N):
        obs, _, _, _ = env.step(controller(obs))
    return env.last_metrics
