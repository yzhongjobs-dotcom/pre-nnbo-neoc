"""Same-trajectory sensor diagnosis; never feeds truth to the controller."""
import json
from pathlib import Path
import numpy as np
import method3 as M
import physics_plant3 as P3
import scenes3 as S3
from sensing_repro import O, CH


def amp(t, x):
    w = 2*np.pi*20.8
    A = np.c_[np.cos(w*t), np.sin(w*t), np.ones(len(t)), t, t*t, t*t*t]
    c = np.linalg.lstsq(A, x, rcond=None)[0]
    return complex(c[0], -c[1])


def job(name):
    sc = S3.cell_scene3('shift_fast', 'relaxed', 'nominal', 1, base=470900, overrides=dict(O, **CH[name]))
    pol = M.build('L', sc)
    original = pol.step.step
    seen = []
    def capture(o):
        seen[:] = [o.log.copy()]
        return original(o)
    pol.step.step = capture
    gates = []
    est = pol.estimator.estimate
    def estimate(o, n, p):
        k, sig, chi = est(o, n, p)
        ref = max(p.k_good, .5*k) if np.isfinite(p.k_good) else k
        gates.append((o.log[n-1, 0], k, sig, chi, sig <= p.params['rel_max']*ref, chi <= p.params['chi2_max'], k > 200.))
        return k, sig, chi
    pol.estimator.estimate = estimate
    log = P3.episode(sc, pol, 801)
    sensed = seen[0]; n = len(sensed); raw = log[:n].copy()
    # Offline copies on identical commands: isolate mass, lag and alignment.
    def corrected(mass, delay, tau):
        p = M.build('L', sc); f = M.ForceSensing(p, mass, delay)
        f.step = lambda o: 0.
        f.reset(dict(P3.PLATFORM, tau_a=tau))
        o = P3.Obs(raw, n); f(o)
        return o.log
    exact = corrected(sc['tool_mass'], sc['sensor_delay_ms'], sc['tau_a_true'])
    mass = corrected(sc['tool_mass'], sc['sensor_delay_hat'], sc['tau_a_hat'])
    # Compare in the regression domain to respect sampling and the force filter.
    j0 = pol.ndel+1; count = n-j0-1
    def build(a):
        return M.common._build(a, j0, count, 0, pol.ns, pol.eb, pol.afc, pol.ndel, 0., 0., 1., True)[1:3]
    Y, z = build(sensed); _, ze = build(exact); _, zm = build(mass)
    t = raw[j0:j0+count, 0]; keep = (t > 1.) & (t < 3.)
    ya = amp(t[keep], Y[keep]); phase = np.conj(ya)/abs(ya)
    out = dict(name=name, metric=M.metric(log, pol, sc), params=pol.params, id_hist=pol.id_hist, true_tau=sc['tau_a_true'], calibrated_tau=sc['tau_a_hat'])
    out['phasors_N'] = {key: [v.real, v.imag] for key, v in
        ((key, amp(t[keep], val[keep])*phase) for key, val in
         [('force', z), ('exact_comp_force', ze), ('mass_residual', z-zm), ('lag_delay_residual', zm-ze)])}
    g = np.array(gates); g = g[g[:, 0] >= 1.]
    out['gates'] = dict(windows=len(g), sigma_fail=int((g[:, 4] == 0).sum()), chi_fail=int((g[:, 5] == 0).sum()),
                        floor_fail=int((g[:, 6] == 0).sum()),
                        sigma_median=float(np.median(g[:, 2])), chi_median=float(np.median(g[:, 3])), chi_limit=pol.params['chi2_max'])
    np.savez_compressed(Path('results/sensing_fix')/('trajectory_'+name.replace(' ', '_')+'_'+LABEL+'.npz'), raw=raw, sensed=sensed, gates=np.array(gates))
    return out


if __name__ == '__main__':
    import sys
    LABEL = sys.argv[1] if len(sys.argv)>1 else 'after'
    if '--baseline' in sys.argv:
        import subprocess
        source = subprocess.check_output(['git', 'show', '2d6e4b8eef169d777e02ae5cd399ea003006e24b:experiments/rossi_failure_study/method3.py'])
        baseline = {}; exec(compile(source, 'baseline_method3.py', 'exec'), baseline)
        M.ForceSensing = baseline['ForceSensing']
    out = [job(n) for n in ('as drawn', 'offset 0', 'offset exact', 'mass exact', 'all exact')]
    p = Path('results/sensing_fix')/('diagnostic_'+LABEL+'.json')
    p.write_text(json.dumps(out, indent=2), encoding='utf-8')
    for r in out:
        print(r['name'], r['gates'], r['phasors_N'], flush=True)
