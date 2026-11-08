"""Monte Carlo check of the derived properties of memory filters. Purely synthetic: ln k = truth + window noise, no plant.

Measurement every hop h=48 ms. Window noise: variance s^2 per window of 0.192 s, modelled as a moving average over 4 hops (overlapping windows), R = 4 s^2 in the filter.
For each filter: stationary sd and its implied memory tau = N/(2 sd^2), lag on ramps, overshoot-free delay on steps.
Theory: Kalman random walk tau = sqrt(N/q), N = s^2*T0; sd^2 = N/(2 tau); ramp lag = rho*tau.
ARD (flat hyperprior): alpha* = max(nu^2-S, 0), gain 1-r/nu^2 above threshold; thresholded ARD: alpha* = max(nu^2-kappa^2 S, 0).
Usage: python -X utf8 memory_theory_check.py
"""
import numpy as np

import memory_filters as MF

S = .2; T0 = .192; H = .048; OV = 4; N = S*S*T0; R = OV*S*S


def noise(rng, n):
    w = rng.standard_normal(n+OV)*S*np.sqrt(OV)
    return np.convolve(w, np.ones(OV)/OV, 'valid')[:n]


def run(make, truth, seed=0, reps=40):
    out = []
    for r in range(reps):
        rng = np.random.default_rng(seed+r)
        f = make(); e = noise(rng, len(truth)); m = np.zeros(len(truth))
        for i in range(len(truth)):
            m[i], _ = f.update(truth[i]+e[i], R, 0., True)
        out.append(m)
    return np.array(out)


def main():
    filters = dict(
        window=lambda: MF.Window(),
        kal_q002=lambda: MF.Kalman(q_t=.002, q_s=0., dt=H, huber=0.),
        kal_q01=lambda: MF.Kalman(q_t=.01, q_s=0., dt=H, huber=0.),
        kal_q002_huber=lambda: MF.Kalman(q_t=.002, q_s=0., dt=H),
        forgetting=lambda: MF.Forgetting(),
        ard_flat=lambda: MF.HuangARD(learn_rho=False),
        ard_learn=lambda: MF.HuangARD(),
        lepski_c2=lambda: MF.Lepski(c=2.),
        tard_2_3=lambda: MF.ThresholdARD(q_t=.002, q_s=0., kappa_up=2., kappa_down=3.),
        tard_25_25=lambda: MF.ThresholdARD(q_t=.002, q_s=0., kappa_up=2.5, kappa_down=2.5),
        tard_15_3=lambda: MF.ThresholdARD(q_t=.002, q_s=0., kappa_up=1.5, kappa_down=3.),
        tard_3_3=lambda: MF.ThresholdARD(q_t=.002, q_s=0., kappa_up=3., kappa_down=3.),
    )
    n = 2000
    print('name | stationary sd | tau=N/(2sd^2) s | ramp .29/s lag | ramp .73/s lag | step .3: delay to half s | step 1.1: delay to half s')
    for name, mk in filters.items():
        st = run(mk, np.zeros(n))[:, 300:]
        sd = float(st.std())
        lags = []
        for rho in (.29, .73):
            tr = np.concatenate([np.zeros(300), rho*H*np.arange(n-300)])[:800]
            m = run(mk, tr)
            lags.append(float(np.mean(tr[600:]-m[:, 600:])))
            if rho == .73:
                rms73 = float(np.sqrt(np.mean((tr[600:]-m[:, 600:])**2)))
        delays = []
        for d in (.3, 1.1):
            tr = np.concatenate([np.zeros(300), d*np.ones(300)])
            m = run(mk, tr).mean(0)
            hit = np.where(m[300:] >= d/2)[0]
            delays.append(float(hit[0]*H) if len(hit) else np.inf)
        print(f"{name:16s} {sd:.3f} {N/(2*sd*sd):6.2f} {lags[0]:.3f} {lags[1]:.3f} rms {rms73:.3f} {delays[0]:.2f} {delays[1]:.2f}", flush=True)
    print('theory: q=.002 tau', np.sqrt(N/.002), 'sd', np.sqrt(N/(2*np.sqrt(N/.002))), 'lag .29', .29*np.sqrt(N/.002), '| q=.01 tau', np.sqrt(N/.01), 'sd', np.sqrt(N/(2*np.sqrt(N/.01))))


if __name__ == '__main__':
    main()
