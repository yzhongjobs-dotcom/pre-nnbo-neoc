"""SLS convolution against independent ODE quadrature, and synthetic recovery."""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
import memory_contact as M


def main():
    results = []
    t = np.arange(72)*.008
    y = .002+.0001*t+.00004*np.sin(2*np.pi*20.8*t)+.00002*np.sin(2*np.pi*7.3*t)
    n = 64; te = t[n-1]; ye = y[n-1]; tw = .256
    measured = np.c_[t, .02*t, y, np.full(len(t), 4.), np.zeros((len(t), 3))]
    jj = np.arange(n-33, n+4)
    for model, tm in M.TIMES.items():
        beta = np.array([4., .03, -.02, .01, 0., .03])
        theta = np.array([12000., 1000.] if model.startswith('drift_') else [12000., 1000., 6000.])
        B, U, _ = M.features(t[jj], y[jj], y[jj+1], te, ye, tw, .019, tm, t[jj[0]])
        if model.startswith('drift_'):
            U = U[:, :2]
        z = B@beta+U@theta; decay = np.exp(-2*np.pi*30*.008)
        for j, v in zip(jj, z):
            measured[j+1, 3] = decay*measured[j, 3]+v
        r = M.fit(measured[:n], model=model)
        np.testing.assert_allclose(r['coefficients'], theta, rtol=1e-6, atol=2e-3)
        for h in (1, 4):
            pred, _ = M.predict(r, measured, n, h, .019, model)
            np.testing.assert_allclose(pred, measured[n-1+h, 3], atol=1e-8)
        # Independently integrate qdot=ydot-q/tm across one command interval.
        delta, q0, ts = .00007, .00003, .019
        speed = delta/(ts*(-np.expm1(-.008/ts)))
        sol = solve_ivp(lambda tt, q: [speed*np.exp(-tt/ts)-q[0]/tm], [0., .008], [q0],
                        t_eval=np.r_[np.arange(8)*.001, .008], rtol=1e-11, atol=1e-14)
        _, testU, qend = M.features(np.array([0.]), np.array([0.]), np.array([delta]), .008, delta, .256, ts, tm, 0., q0)
        af = -np.expm1(-2*np.pi*30*.001); weights = af*(1-af)**(7-np.arange(8))
        np.testing.assert_allclose(testU[0, 2], weights@sol.y[0, :8], rtol=1e-8, atol=1e-12)
        np.testing.assert_allclose(qend, sol.y[0, -1], rtol=1e-8, atol=1e-12)
        results.append(dict(model=model, recovery=True, conditional_prediction=True, independent_ode=True))
    # Physical closure: moving cubic surface + initially nonequilibrated SLS.
    # Generate raw contact force independently, then sensor-filter at 1 ms.
    tm, ts = .2, .019
    samples = np.full_like(measured, 0.); samples[:, :3] = measured[:, :3]
    samples[0, 3] = 4.
    q = .000041  # unobserved, not equal to an equilibrium indentation
    def surface(tt):
        return .00175+.00007*tt+.0001*tt**3
    for j in range(len(t)-1):
        speed = (y[j+1]-y[j])/(ts*(-np.expm1(-.008/ts)))
        sol = solve_ivp(lambda h, state: [speed*np.exp(-h/ts)-(.00007+.0003*(t[j]+h)**2)-state[0]/tm],
                        [0., .008], [q], t_eval=np.r_[np.arange(8)*.001, .008], rtol=1e-11, atol=1e-14)
        force_filtered = samples[j, 3]
        for i, h in enumerate(np.arange(8)*.001):
            yy = y[j]+speed*ts*(-np.expm1(-h/ts))
            ke = 12000.+1000.*(t[j]+h-te)/tw
            force = ke*(yy-surface(t[j]+h))+6000.*sol.y[0, i]
            force_filtered += af*(force-force_filtered)
        samples[j+1, 3] = force_filtered
        q = sol.y[0, -1]
    physical_fit = M.fit(samples[:n], model='memory_200')
    np.testing.assert_allclose(physical_fit['coefficients'], [12000., 1000., 6000.], rtol=2e-6, atol=.01)
    results.append(dict(moving_cubic_surface=True, nonequilibrium_initial_memory=True,
                        recovered=physical_fit['coefficients'].tolist(), independent_physical_force_generation=True))
    out = Path(__file__).resolve().parent/'results/contact_memory_verification.json'
    out.write_text(json.dumps(results, indent=2), encoding='utf8')
    print('memory recovery, independent ODE, held-out state continuation passed')


if __name__ == '__main__':
    main()
