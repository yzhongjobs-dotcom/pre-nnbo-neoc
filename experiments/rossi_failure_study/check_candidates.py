"""Acceptance of contact_candidates.py before any comparison.

Regressions: L equals the earlier host on the zero-physics plant; zero phase reduces W, S, H to L; hiding truth columns changes nothing.
Mechanism: online phase against the plant's measured probe phase; identified friction against the set values; coupling c against the
geometric value tan(theta) + xi k/K_t in reversal scenes.
"""
import json
import math
from pathlib import Path

import numpy as np

import sweeps as S
import physics_plant as PP1
import physics_plant2 as PP2
import scenes2 as SC2
import contact_candidates as C
from contact_compare import ModelPolicy, RIDGE, CTL
from failure_replay import Measured, cell_scene
from ridge_estimator import RidgeProjection

HERE = Path(__file__).resolve().parent


def zero_physics(sc):
    sc = dict(sc)
    sc.update(mu_k=0., mu_s=0., k_d=0., k1=0., k2=0., theta0=0., c_slope=0., xi_c=0., K_t=1e12, v_s=.005, v_e=.0005, nu=.5, m_t=1., zeta_t=.5)
    t = np.arange(int(sc['T']*1000)+2)*1e-3
    sc['xcmd'] = PP2.X0+sc['v']*np.maximum(t-.008*PP2.PLATFORM['ndel'], 0.)    # the old plant applies the command delay to x
    return sc


def plant_phase(sc, tick=20.8):
    """Phase and gain of the plant's response to a probe, measured from noise-free replay with the motion drives removed."""
    s = dict(sc)
    for k in ('u', 'fext', 'act', 'rip'):
        s[k] = np.zeros_like(sc[k])
    s0, _ = PP2._init(s, PP2.PLATFORM)
    n = int(s['T']/.008); t = np.arange(n)*.008
    log = PP2.replay(s, s0[2]+1e-3/max(1., 1.)*0+(.15/2000.)*np.sin(2*np.pi*tick*t))
    m = log[:, 0] > 2.
    a = np.c_[np.sin(2*np.pi*tick*log[m, 0]), np.cos(2*np.pi*tick*log[m, 0]), np.ones(m.sum()), log[m, 0]]
    cy = np.linalg.lstsq(a, log[m, 2], rcond=None)[0]; cf = np.linalg.lstsq(a, log[m, 18], rcond=None)[0]
    g = complex(cf[0], cf[1])/complex(cy[0], cy[1])
    return float(np.angle(g)), abs(g)


def main():
    out = {}
    # regressions
    sc = cell_scene('shift_fast', 'mid', 2, 420000); seed = 420002
    p = dict(S.FINAL_C); p.update(CTL)
    a = ModelPolicy(RidgeProjection(**RIDGE), single_model=False, **p)
    la = PP1.episode(sc, Measured(a), seed)
    z = zero_physics(sc)
    lb, _, _ = C.run_one('L', z, seed)
    n = min(len(la), len(lb)); print('L against earlier host, command max gap', float(np.max(np.abs(la[:n, 11]-lb[:n, 11]))))
    for name in ('W', 'S', 'H'):
        pol = C.build(name); pol.update_phase = lambda o, n: None
        lc = PP2.episode(z, C.Measured2(pol), seed)
        print(name, 'with phase frozen at 0, gap to L', float(np.max(np.abs(lc[:n, 11]-lb[:n, 11]))), '(H keeps its F^.206, so a gap there is expected)')
    # mechanism
    rows = []
    for regime in ('nominal', 'rate', 'memory'):
        for tissue in ('soft', 'mid', 'hard'):
            sc = SC2.cell_scene2('rest', tissue, regime, 0, base=421000)
            true, gain = plant_phase(sc)
            log, pol, _ = C.run_one('W', sc, 1)
            rows.append(dict(regime=regime, tissue=tissue, phi_true=true, phi_hat=pol.phi, gain=gain, updates=pol.phi_count))
            print(f"{regime:8s} {tissue:5s} phase plant {math.degrees(true):6.1f} deg, online {math.degrees(pol.phi):6.1f} deg, updates {pol.phi_count}, gain {gain:.0f}, k_good {pol.k_good:.0f}")
    out['phase'] = rows
    rows = []
    for i in range(6):
        sc = SC2.cell_scene2('rest', 'mid', 'reversal', i, base=421000)
        for name in ('AS', 'TR', 'LF'):
            log, pol, _ = C.run_one(name, sc, 1)
            k = pol.k_good; target = math.tan(float(np.median(log[200:, 17])))+sc['xi_c']*k/sc['K_t']
            rows.append(dict(i=i, name=name, mu=sc['mu_k'], mu_hat=pol.mu_hat, kd=sc['k_d'], kd_hat=pol.kd_hat, c_target=target, c_hat=pol.c_hat))
            print(f"reversal {i} {name:3s} mu {sc['mu_k']:.3f} hat {pol.mu_hat:.3f} | k_d {sc['k_d']:5.2f} hat {pol.kd_hat:5.2f} | c target {target:+.3f} hat {pol.c_hat:+.3f} | RMSE Fn {np.sqrt(np.mean((log[125:, 3]-4)**2)):.3f}")
    out['friction'] = rows
    (HERE/'results'/'candidate_checks.json').write_text(json.dumps(out, indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
