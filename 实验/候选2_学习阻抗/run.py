"""Predeclared quick screen. Run from any directory; no B0 optimizer is claimed."""
import argparse
import json
import platform
import time
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import differential_evolution, minimize
import torch
from stable_baselines3 import PPO
from stable_baselines3.common.callbacks import BaseCallback
from contact_env import ROOT, BASE, cases, evaluate, baseline, ContactVec, action_from_target, advance

OUT = ROOT / 'results'


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False), encoding='utf8')


def summary(a):
    return dict(mean_cost=float(np.mean(a[:, 0])), median_cost=float(np.median(a[:, 0])),
                numerical_pass_fraction=float(np.mean(a[:, 1]>=-1e-7)),
                minimum_margin=float(a[:, 1].min()), mean_peak_force_N=float(a[:, 2].mean()),
                force_rms_N=float(np.sqrt(a[:, 3].mean()/4)),
                position_rms_mm=float(1000*np.sqrt(a[:, 4].mean()/4)),
                shield_intervention_fraction=float(a[:, 5].mean()/32))


class Select(BaseCallback):
    def __init__(self, history, seed, fixed):
        super().__init__()
        self.history, self.seed, self.fixed = history, seed, fixed
        self.best = np.inf
        self.rows = []
        self.started = time.perf_counter()

    def _on_step(self):
        if self.num_timesteps % 131072 == 0:
            result = evaluate(lambda obs: self.model.predict(obs, deterministic=True)[0], self.fixed, self.history)
            score = np.mean(result[:, 0]+100*np.maximum(-result[:, 1], 0))
            row = dict(steps=self.num_timesteps, seconds=time.perf_counter()-self.started,
                       selection_score=float(score), **summary(result))
            self.rows.append(row)
            if score<self.best:
                self.best = score
                self.model.save(OUT/f'ppo_h{self.history}_s{self.seed}')
            dump(OUT/f'learning_h{self.history}_s{self.seed}.json', self.rows)
            print('validation', self.history, self.seed, row, flush=True)
        return True


def calibrate(only=None):
    fixed = cases(1101, 64)
    valid = cases(2202, 64)
    # Same deployment-side observations and same deterministic sensor streams.
    bounds = dict(fixed=[(20, 1000), (1, 25)],
                  feedback=[(0, 1150), (0, 27), (-9.21, -4.61), (-6.91, -2.30), (-4, 4)],
                  preview=[(.15, 1.2), (0., .8), (.5, 1.5)],
                  prior_preview=[(.15, 1.2), (0., .8), (.5, 1.5)])
    out = json.loads((OUT/'baselines.json').read_text()) if only else {}
    if only:
        bounds = {only:bounds[only]}
    for kind, bound in bounds.items():
        tic = time.perf_counter()
        def obj(x):
            a = evaluate(lambda obs: baseline(obs, x, kind, 4), fixed)
            return np.mean(a[:, 0]+100*np.maximum(-a[:, 1], 0))
        r = differential_evolution(obj, bound, seed=3303, popsize=6, maxiter=18, tol=1e-4, polish=False)
        polished = minimize(obj, r.x, bounds=bound, method='Powell', options=dict(maxiter=35, ftol=1e-5, xtol=1e-4))
        options = [r.x, polished.x]
        vals = [evaluate(lambda obs: baseline(obs, x, kind, 4), valid) for x in options]
        scores = [np.mean(a[:, 0]+100*np.maximum(-a[:, 1], 0)) for a in vals]
        pick = int(np.argmin(scores))
        out[kind] = dict(parameters=options[pick].tolist(), validation=summary(vals[pick]),
                         calibration_seconds=time.perf_counter()-tic,
                         objective_evaluations=int(r.nfev+polished.nfev))
        dump(OUT/'baselines.json', out)
        print('baseline', kind, out[kind], flush=True)


def train():
    valid = cases(2202, 64)
    for history in [4, 1]:
        for seed in [11, 22, 33]:
            env = ContactVec(64, history, seed)
            model = PPO('MlpPolicy', env, seed=seed, n_steps=64, batch_size=512, n_epochs=6,
                        learning_rate=3e-4, gamma=1., gae_lambda=.95, ent_coef=.005,
                        policy_kwargs=dict(net_arch=dict(pi=[64, 64], vf=[64, 64]),
                                           log_std_init=-.7), verbose=0, device='cpu')
            with torch.no_grad():
                model.policy.action_net.bias.copy_(torch.tensor(action_from_target(np.array([50., 3.])), dtype=torch.float32))
            cb = Select(history, seed, valid)
            model.learn(total_timesteps=524288, callback=cb)
            print('completed', history, seed, 'seconds', time.perf_counter()-cb.started, flush=True)


def test():
    base = json.loads((OUT/'baselines.json').read_text())
    groups = dict(id=cases(4404, 96), material_ood=cases(5505, 96, 'material_ood'),
                  noise_mass_ood=cases(6606, 96, 'noise_mass_ood'))
    controllers = {kind:(lambda obs, k=kind: baseline(obs, base[k]['parameters'], k, 4), 4) for kind in base}
    for h in [4, 1]:
        for seed in [11, 22, 33]:
            m = PPO.load(OUT/f'ppo_h{h}_s{seed}', device='cpu')
            controllers[f'ppo_h{h}_s{seed}'] = (lambda obs, model=m: model.predict(obs, deterministic=True)[0], h)
    allmetrics = {}; results = {}
    for group, fixed in groups.items():
        results[group] = {}
        for name, (fn, h) in controllers.items():
            a = evaluate(fn, fixed, h)
            allmetrics[f'{group}__{name}'] = a
            results[group][name] = summary(a)
        for name in controllers:
            a = allmetrics[f'{group}__{name}']
            for rival in base:
                ratio = a[:, 0]/allmetrics[f'{group}__{rival}'][:, 0]-1
                results[group][name][f'excess_vs_{rival}'] = dict(mean=float(ratio.mean()),
                       median=float(np.median(ratio)), p90=float(np.quantile(ratio, .9)))
        print('test', group, results[group], flush=True)
    np.savez_compressed(OUT/'paired_metrics.npz', **allmetrics)
    dump(OUT/'test.json', results)
    # Fine-step integration is evaluated with the SAME sampled policy frequency.
    convergence = {}
    tiny = tuple(a[:12].copy() for a in groups['material_ood'])
    for name, (fn, h) in controllers.items():
        a = evaluate(fn, tiny, h)
        b = evaluate(fn, tiny, h, sub=100)
        convergence[name] = dict(max_relative_J_error=float(np.max(np.abs(a[:, 0]/b[:, 0]-1))),
                                 fine_minimum_margin=float(b[:, 1].min()))
    dump(OUT/'convergence.json', convergence)
    timing = {}
    for name, (fn, h) in controllers.items():
        env = ContactVec(1, h, 99)
        obs = env.reset()
        fn(obs)
        ticks = []
        for _ in range(500):
            tic = time.perf_counter_ns(); fn(obs); ticks.append((time.perf_counter_ns()-tic)*1e-6)
        timing[name] = dict(batch1_median_ms=float(np.median(ticks)), batch1_p95_ms=float(np.quantile(ticks, .95)))
    dump(OUT/'inference_timing.json', timing)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('stage', choices=['calibrate', 'train', 'test'])
    p.add_argument('--only', choices=['fixed', 'feedback', 'preview', 'prior_preview'])
    args = p.parse_args()
    OUT.mkdir(exist_ok=True)
    torch.set_num_threads(1)
    if not (OUT/'protocol.json').exists():
        dump(OUT/'protocol.json', dict(parent='Q2 continuous projected actual-gain carrier',
             training_steps_per_run=524288, seeds=[11, 22, 33], histories=[4, 1],
             checkpoints=[131072, 262144, 393216, 524288], train_calibration_cases=64,
             validation_seed=2202, validation_cases=64, tests=dict(id=4404, material_ood=5505, noise_mass_ood=6606),
             test_cases_per_group=96, python=platform.python_version(), numpy=np.__version__,
             scipy=scipy.__version__, torch=torch.__version__,
             objective='Q2 physical J; training adds 10 * negative step-margin penalty; validation selection adds 100 * negative episode-margin penalty',
             selection='minimum penalized validation cost, no test-driven retraining',
             known_task_input='two surface pulse programs; no material truth in actor/critic or causal baselines',
             reduced_scope='normal e,v,K,D,E, not full B0/MHE/4-gain/robot/cylinder; no sim2real claim'))
    if args.stage == 'calibrate':
        calibrate(args.only)
    else:
        globals()[args.stage]()


if __name__ == '__main__':
    main()
