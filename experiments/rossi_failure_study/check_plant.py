"""Independent checks on the physics plant before any method is compared. No method ranking here."""
import json, math, os, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
for s in ('motion_aware_scan', 'method_redesign', 'contact_information', 'cobot_scan_benchmark'):
    sys.path.insert(0, str(HERE.parent/s))
import physics_plant as PP
import scenes as SC
from scipy import signal


class Hold:
    direct = True
    def __init__(self, f=None): self.f = f; self.t_call = []
    def reset(self, p): self.n = 0
    def __call__(self, o):
        self.n += 1
        return float(o.log[o.n-1, 11]) if self.f is None else float(self.f(o))


def quiet_scene(**kw):
    rng = np.random.default_rng(1)
    sc = SC.make('rest', 'mid', rng, T=6., v=.0, overrides=dict(field_sigma=0.))
    sc['h_s'][:] = 0.; sc.update(kw)
    return sc


def main():
    out = {}
    # 1 equilibrium: a held command at the target force should stay near it, with slow relaxation only.
    sc = quiet_scene(); log = PP.episode(sc, Hold(), 3, sensor=False)
    out['equilibrium_force_start_end'] = [float(log[5, 3]), float(log[-1, 3])]
    # 2 dynamic gain: a small sine on the command; the measured transfer should equal the plant's scored gain.
    res = []
    for f, in [(1.,), (5.,), (20.8,), (40.,)]:
        sc = quiet_scene()
        A = 2e-6
        y0 = {}
        def cmd(o, y0=y0):
            if 'y' not in y0: y0['y'] = float(o.log[0, 11])
            return y0['y']+A*math.sin(2*np.pi*f*o.log[o.n-1, 0])
        log = PP.episode(sc, Hold(cmd), 3, sensor=False)
        m = log[:, 0] > 2.
        t = log[m, 0]; X = np.c_[np.sin(2*np.pi*f*t), np.cos(2*np.pi*f*t), np.ones(len(t)), t]
        cy = np.linalg.lstsq(X, log[m, 2], rcond=None)[0]; cf = np.linalg.lstsq(X, log[m, 3], rcond=None)[0]
        G = complex(cf[0], cf[1])/complex(cy[0], cy[1])
        res.append(dict(f=f, measured_inphase=G.real, measured_quad=G.imag, scored_probe=float(np.median(log[m, 8])), scored_1hz=float(np.median(log[m, 9])), inst=float(np.median(log[m, 10]))))
    out['dynamic_gain'] = res
    # 3 determinism: replaying logged commands reproduces the encoder exactly, with and without the sensor noise.
    rng = np.random.default_rng(2); sc = SC.make('combo', 'mixed', rng, T=8.)
    log = PP.episode(sc, Hold(), 5, sensor=True); rep = PP.replay(sc, log[:, 11], seed=5, sensor=True)
    out['replay_gap_force_encoder'] = [float(np.max(np.abs(rep[:, 5]-log[:, 5]))), float(np.max(np.abs(rep[:, 2]-log[:, 2])))] if rep is not None and len(rep) == len(log) else 'lost'
    # 4 reaction: a held probe on a stiff contact must push the body by about F/k.
    sc = quiet_scene(); log = PP.episode(sc, Hold(), 3, sensor=False)
    out['body_sag_vs_F_over_k'] = [float(log[-1, 4]), float(-4./sc['k'])]
    # 5 body response: free-response frequency to a force pulse compared with the design resonance.
    sc = quiet_scene(); sc['fext'] = np.zeros_like(sc['fext']); sc['fext'][1500:1520] = 200.; sc['c'] *= .1
    log = PP.episode(sc, Hold(), 3, sensor=False)
    z = log[:, 4]-log[-1, 4]; seg = z[(log[:, 0] > 1.6) & (log[:, 0] < 4.)]
    f, p = signal.periodogram(seg-seg.mean(), fs=125.)
    out['body_peak_hz_vs_design'] = [float(f[np.argmax(p)]), float(math.sqrt(sc['k']/sc['m'])/(2*np.pi))]
    # 6 not a cubic: fraction of event-window surface motion a cubic cannot describe; and in-band energy.
    rows = {}
    for cond in ('breath', 'shift_slow', 'shift_fast', 'flinch', 'cough', 'tremor'):
        vals = []; band = []; base = []
        for s in range(8):
            rng = np.random.default_rng(100+s); sc = SC.make(cond, 'mid', rng, T=10., overrides=dict(breath_amp=.003))
            log = PP.episode(sc, Hold(), s, sensor=False)
            if log is None:
                continue
            z = log[:, 4]; t = log[:, 0]
            ev = sc['events'][0][1] if sc['events'] else 3.
            m = (t > ev-.05) & (t < ev+.5)
            r = []
            for i in range(np.flatnonzero(m)[0], np.flatnonzero(m)[-1]-16, 4):
                w = z[i:i+16]; x = np.linspace(-1, 1, 16)
                r.append(np.std(w-np.polyval(np.polyfit(x, w, 3), x))/max(np.ptp(w), 1e-9))
            f, p = signal.welch(z[m]-z[m].mean(), fs=125., nperseg=min(64, m.sum()))
            vals.append(float(np.median(r))); band.append(float(np.sqrt(p[(f > 16) & (f < 26)].sum()*(f[1]-f[0]))))
            base.append(float(np.ptp(z[m])))
        rows[cond] = dict(cubic_residual_fraction=float(np.median(vals)), body_motion_range_mm=1e3*float(np.median(base)), rms_in_probe_band_um=1e6*float(np.median(band)))
    out['motion_character'] = rows
    print(json.dumps(out, indent=1))
    (HERE/'results').mkdir(exist_ok=True); (HERE/'results'/'plant_checks.json').write_text(json.dumps(out, indent=1), encoding='utf8')


if __name__ == '__main__':
    main()
