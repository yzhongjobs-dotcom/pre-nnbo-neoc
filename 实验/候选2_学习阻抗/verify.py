"""Independent DOP853 equations and a zero-excitation observability check."""
import json
import time
import numpy as np
from scipy.integrate import solve_ivp
import torch
from stable_baselines3 import PPO
from contact_env import ROOT, LO, HI, H, N, ContactVec, cases, evaluate, baseline, action_from_target


def independent(fn, history, fixed):
    env = ContactVec(1, history, 97531, fixed=fixed)
    obs = env.reset()
    y = env.states[0].copy()
    p = fixed[0][0]; program = fixed[1][0]
    total = 0.

    def external(t):
        theta = p[4]+p[5]*(1-np.cos(np.pi*np.clip(t/4, 0, 1)))/2
        phase = (t-program[:, 0])/program[:, 1]
        gamma = np.sum(np.where((phase>0)&(phase<1), program[:, 2]*np.sin(np.pi*phase)**2, 0))
        return theta, gamma

    for j in range(N):
        a = np.clip(fn(obs)[0], -1, 1)
        target = np.exp(np.log(LO)+(a+1)/2*np.log(HI/LO))
        desired = np.clip((target-y[2:4])/H, -p[6:8], p[6:8])

        def rhs(t, x):
            e, v, k, d, E = x[:5]
            lower = np.maximum(-p[6:8], -20*(x[2:4]-LO))
            upper = np.minimum(p[6:8], 20*(HI-x[2:4]))
            if e*e>1e-25:
                upper[0] = min(upper[0], (p[10]*(p[1]+d)*v*v+p[11]*(E-p[8]))/(.5*e*e))
            rate = np.maximum(lower, np.minimum(desired, upper))
            theta, gamma = external(t)
            fe = theta*(gamma-e); fr = -k*e-d*v
            cost = p[18]*(e/p[14])**2+p[19]*(v/p[15])**2+p[20]*(fe/p[16])**2
            cost += p[21]*(fr/p[17])**2+p[22]*(((k-p[24])/p[26])**2+((d-p[25])/p[27])**2)
            cost += p[23]*np.sum((rate/p[6:8])**2)
            return [v, (fr+fe-p[1]*v)/p[0], *rate,
                    p[10]*(p[1]+d)*v*v-.5*e*e*rate[0], cost]

        result = solve_ivp(rhs, [j*H, (j+1)*H], np.r_[y, 0.], method='DOP853',
                           rtol=2e-10, atol=2e-12, max_step=.0025)
        if not result.success:
            raise RuntimeError(result.message)
        y = result.y[:5, -1]; total += result.y[5, -1]
        env.states[0] = y; env.index[0] = j+1
        env.frames[:, :-1] = env.frames[:, 1:]
        env.frames[:, -1] = env._frame()
        obs = env._obs()
    th, ga = external(4.)
    total += p[29]*((y[0]/p[14])**2+p[20]*(th*(ga-y[0])/p[16])**2)
    return total


def main():
    torch.set_num_threads(1)
    out = ROOT/'results'
    bs = json.loads((out/'baselines.json').read_text())
    policies = {k:(lambda obs, k=k: baseline(obs, bs[k]['parameters'], k, 4), 4) for k in bs}
    for h in [4, 1]:
        m = PPO.load(out/f'ppo_h{h}_s11')
        policies[f'ppo_h{h}_s11'] = (lambda obs, m=m: m.predict(obs, deterministic=True)[0], h)
    fixed = cases(5505, 3, 'material_ood')
    rows = {}
    for name, (fn, h) in policies.items():
        a = evaluate(fn, fixed, h)
        refs = []
        for i in range(3):
            # Match each case's sensor stream to the batched evaluator by replaying
            # a one-case evaluator on both sides; batch RNG order differs.
            one = tuple(x[i:i+1].copy() for x in fixed)
            computed = evaluate(fn, one, h)[0, 0]
            direct = independent(fn, h, one)
            refs.append(dict(rk4_cost=float(computed), dop853_cost=float(direct),
                             relative_error=float(computed/direct-1)))
        rows[name] = refs
        print('independent', name, refs, flush=True)
    # Two materials are indistinguishable before any excitation in this carrier.
    one = cases(9909, 1)
    other = tuple(x.copy() for x in one)
    other[0][0, 4:6] = [600, 650]
    a = ContactVec(1, 4, 77, fixed=one).reset()
    b = ContactVec(1, 4, 77, fixed=other).reset()
    rows['resting_observations_max_difference'] = float(np.max(np.abs(a-b)))
    assert np.array_equal(a, b)
    rows['maximum_independent_relative_error'] = max(abs(v['relative_error']) for k, rs in rows.items() if isinstance(rs, list) for v in rs)
    assert rows['maximum_independent_relative_error']<.001
    (out/'independent_checks.json').write_text(json.dumps(rows, indent=2), encoding='utf8')


if __name__ == '__main__':
    main()
