"""Isolate the calibrated mass error on a fixed sinusoidal trajectory."""
import math
import json
from pathlib import Path
import numpy as np
import method3 as M
import physics_plant3 as P3
import scenes3 as S3
from sensing_repro import O
from diagnose_sensing import amp


def main():
    rows = []
    for delay in (0., 1.88):
        sc = S3.cell_scene3('shift_fast', 'relaxed', 'nominal', 1, T=4., base=470900,
                             overrides=dict(O, sensor_delay_ms=delay, delay_error=0., tau_error=0.))
        s, _ = P3._init(sc, P3.PLATFORM)
        t = np.arange(500)*.008
        raw = P3.replay(sc, s[2]+.0002*np.sin(2*np.pi*20.8*t))
        zs = []
        for error in (0., .047):
            p = M.build('L', sc)
            f = M.ForceSensing(p, sc['tool_mass']*(1+error), delay)
            f.step = lambda o: 0.
            f.reset(dict(P3.PLATFORM, tau_a=sc['tau_a_true']))
            o = P3.Obs(raw, len(raw)); f(o)
            _, Y, z, _ = M.common._build(o.log, 3, len(raw)-4, 0, 8, p.eb, p.afc, 1, 0., 0., 1., True)
            zs.append(z)
        tt = raw[3:-1, 0]; sel = tt > 1.
        k = amp(tt[sel], (zs[1]-zs[0])[sel])/amp(tt[sel], Y[sel])
        rows.append(dict(delay_ms=delay, mass_error=.047, mass_residual_equivalent_gain_N_m=[k.real, k.imag],
                         mass_residual_phase_deg=float(np.angle(k)*180/np.pi),
                         continuous_mass_gain_N_m=.047*sc['tool_mass']*(2*np.pi*20.8)**2))
    Path('results/sensing_fix/mass_phase.json').write_text(json.dumps(rows, indent=2), encoding='utf-8')
    print(json.dumps(rows, indent=2))


if __name__ == '__main__':
    main()
