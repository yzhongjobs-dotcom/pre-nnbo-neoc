"""Capture, per tick, the exact regression window (Y, z) the estimators saw in the Huber-D4 closed loop.

Noisy windows are what every estimator in table 2 receives. A second pass replays the same commands through
the plant with sensor noise off, so z - z0 and Y - Y0 are the sensor noise alone and (Y0, z0) carry only the plant.
Prior-level analysis then needs no plant: huber_prior_check.py reads results/huber_windows/*.npz.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_key] = '1'
import sys
sys.dont_write_bytecode = True
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import corrected_comparison as CC
import method3 as M
import overlay_compare as O
import memory_segment_eval as SE
from d_robust import common

P3, B = SE.P3, O.B
HERE = Path(__file__).resolve().parent
OUT = HERE/'results/huber_windows'
NW = 32
COLS = [0, 1, 2, 3, 5, 11, 14, 18, 21, 26, 27]


class Replay:
    """Open-loop replay of recorded commands through the same sensing stage the Method's estimator reads from."""
    def __init__(self, commands, pol, sc):
        self.commands = commands; self.i = 0; self.seen = []
        self.Tc = pol.Tc; self.ndel = pol.ndel; self.afc = pol.afc
        self.sensing = M.ForceSensing(self, sc.get('tool_mass_hat', 0.), sc.get('sensor_delay_hat', 0.), interp=sc.get('cmd_interp', 0.) > 0., t2=sc.get('actuator_t2'))

    def reset(self, platform):
        self.i = 0; self.seen = []

    def step(self, o):
        self.seen.append(o.log[o.n-1].copy())
        u = self.commands[min(self.i, len(self.commands)-1)]; self.i += 1
        return u


class ReplayPolicy:
    def __init__(self, replay):
        self.replay = replay

    def reset(self, platform):
        self.replay.sensing.reset(platform)

    def __call__(self, o):
        return self.replay.sensing(o)


def hc_log(log, tau_h, mag_h, tc):
    """The H preprocessing of contact_candidates.Preprocess._log with a recorded tau_h, mag_h per row."""
    out = log.copy()
    y_last = float(log[0, 2]); y_out = y_last; prev = 0.
    for j in range(1, len(log)):
        y = float(log[j, 2]); dy = y-y_last; y_last = y
        y_out += (dy+tau_h[j]*(dy-prev)/tc)/mag_h[j]
        prev = dy
        out[j, 2] = y_out
    return out


def time_weighted(log, j0, m, ns, eb, afc, te, ye, enc):
    """h1 of the V2 note, eq. (1): the same sub-samples, weights and encoder interpolation as common._build, each times (t_jm - T)."""
    a = 1.-afc; ebn = eb**ns
    h1 = np.zeros(m)
    for jj in range(m):
        j = j0+jj
        for mm in range(1, ns+1):
            w = afc*a**(ns-mm)
            fr = (1.-eb**(mm-1))/(1.-ebn) if ebn < 1. else (mm-1)/ns
            if enc:
                p = log[j, 2]+fr*(log[j+1, 2]-log[j, 2])-ye
            else:
                raise ValueError('only the encoder form is used')
            h1[jj] += w*(log[j, 0]+(mm-1)*1e-3-te)*p
    return h1


def windows(log, policy, enc, eb):
    """eb[n-1]: the online servo estimate in force when the estimator was called at n. Also returns the time-weighted regressor h1."""
    n_all = len(log)
    Y = np.full((n_all, NW), np.nan); z = np.full((n_all, NW), np.nan); Y1 = np.full((n_all, NW), np.nan)
    for n in range(1, n_all+1):
        j0 = max(policy['ndel']+1, n-1-NW); m = n-2-j0+1
        if m < 8:
            continue
        _, y1, z1, keep = common._build(log, j0, m, 0, policy['ns'], eb[n-1], policy['afc'], policy['ndel'], log[n-1, 0], log[n-1, 2], NW*policy['ns']*1e-3, enc)
        t1 = time_weighted(log, j0, m, policy['ns'], eb[n-1], policy['afc'], log[n-1, 0], log[n-1, 2], enc)
        broken = np.flatnonzero(~keep)
        if len(broken):
            y1, z1, t1 = y1[broken[-1]+1:], z1[broken[-1]+1:], t1[broken[-1]+1:]
        if len(y1):
            Y[n-1, NW-len(y1):] = y1; z[n-1, NW-len(y1):] = z1; Y1[n-1, NW-len(y1):] = t1
    return Y, z, Y1


def job(case):
    O.setup()
    sc = SE.MC.scene_for(tuple(case['cell']), CC.BASE, case['variant'])
    pol = B.build('Huber-D4', O.chosen()['Huber-D4'], sc)
    est = pol.estimator; inner = est.estimate
    sent = []; seen = []; ebs = {}; tau = {}; mag = {}; raw = {}; policy = {}

    def record(o, n, p):
        policy.update(ndel=p.ndel, ns=p.ns, eb=p.eb, afc=p.afc, enc=bool(p.params['enc_y']), noise=float(p.params['noise']), Tc=p.Tc, a_force=float(p.params['a_force']))
        ebs[n-1] = float(p.eb)
        out = inner(o, n, p)
        raw[n-1] = out[0]
        return out

    est.estimate = record
    inner_log = est._log

    def record_log(o, n):
        # H coefficients in force when each row enters the preprocessed log, whoever asks first.
        own = est.owner
        for j in range(est.count, n):
            seen.append(o.log[j].copy())
            tau[j] = float(own.tau_h); mag[j] = float(own.mag_h)
        return inner_log(o, n)

    est._log = record_log

    class Rec:
        def reset(self, platform):
            pol.reset(platform)

        def __call__(self, o):
            u = pol(o); sent.append(float(u))
            return u

    log = P3.episode(sc, Rec(), case['seed'])
    n = len(log)
    tau_h = np.array([tau.get(j, np.nan) for j in range(n)]); mag_h = np.array([mag.get(j, np.nan) for j in range(n)])
    tau_h = np.nan_to_num(tau_h, nan=0.); mag_h = np.where(np.isfinite(mag_h), mag_h, 1.)
    # The rows as the estimator first read them: after ForceSensing, before later ticks finish them.
    rows = np.asarray(seen)
    noisy = hc_log(rows, tau_h[:len(rows)], mag_h[:len(rows)], policy['Tc'])
    used = est.private[:len(rows), 2]; last = max(tau)
    worst = float(np.max(np.abs(noisy[:last+1, 2]-used[:last+1])))
    assert worst < 1e-12, worst
    eb = np.array([ebs.get(j, np.nan) for j in range(len(rows))])
    eb = np.where(np.isfinite(eb), eb, policy['eb'])
    Y, z, H1 = windows(noisy, policy, policy['enc'], eb)
    # Same commands, sensor noise off, same preprocessing coefficients.
    replay = Replay(sent, pol, sc)
    log0 = P3.episode(dict(sc, noise_scale=0.), ReplayPolicy(replay), case['seed'])
    # The estimator is not called after a final contact loss, so its rows may stop before the episode ends.
    rows0 = np.asarray(replay.seen)[:len(rows)]
    assert len(rows0) == len(rows), (len(rows0), len(rows))
    Y0, z0, H10 = windows(hc_log(rows0, tau_h[:len(rows)], mag_h[:len(rows)], policy['Tc']), policy, policy['enc'], eb)
    # Same commands, noise off, sensing stage given the true calibration: what the calibration error alone contributes.
    ideal = dict(sc, noise_scale=0.); CC.perfect_calibration(ideal)
    replay_i = Replay(sent, pol, ideal)
    P3.episode(ideal, ReplayPolicy(replay_i), case['seed'])
    rows_i = np.asarray(replay_i.seen)[:len(rows)]
    Y0i, z0i, _ = windows(hc_log(rows_i, tau_h[:len(rows)], mag_h[:len(rows)], policy['Tc']), policy, policy['enc'], eb)
    meta = dict(index=case['index'], group=case['group'], cell=case['cell'], seed=case['seed'], variant=case['variant'],
                event_info=sc['event_info'], cols=COLS, policy=policy, preprocess_max_error=worst,
                K_t=float(sc['K_t']), kr=float(sc['kr']), cover=float(sc['cover']),
                tool_mass=float(sc['tool_mass']), tool_mass_hat=float(sc['tool_mass_hat']), delay=float(sc['sensor_delay_ms']), delay_hat=float(sc['sensor_delay_hat']), tau_a=float(sc['tau_a_true']), tau_a_hat=float(sc['tau_a_hat']))
    np.savez_compressed(OUT/f"case_{case['index']:02d}.npz", Y=Y.astype(np.float64), z=z.astype(np.float64), H1=H1, H10=H10, Y0=Y0.astype(np.float64), z0=z0.astype(np.float64), Y0i=Y0i.astype(np.float64), z0i=z0i.astype(np.float64),
                        log=log[:len(rows), COLS], log0=log0[:len(rows), COLS], raw=np.array([raw.get(j, np.nan) for j in range(len(rows))]), meta=np.asarray(json.dumps(meta, default=float)))
    return meta


if __name__ == '__main__':
    OUT.mkdir(exist_ok=True)
    with ProcessPoolExecutor(8) as ex:
        metas = list(ex.map(job, CC.cases(), chunksize=1))
    print(len(metas), 'cases; noise', metas[0]['policy']['noise'], 'a_force', metas[0]['policy']['a_force'], 'policy', metas[0]['policy'], flush=True)
