"""Does the dynamic truth hold, against the records it came from and against what it did not see, and is it a twin of any candidate?

A. the tables on the plant's grid against the four Neumann protocols, including peak decay and dissipation over cycles, which were not fitted as such
B. other records, not used in the fit: Greiner human brain, low frequency
C. the plant's own code: open-loop sine and step against the analytic response
D. each candidate form, matched to the truth's phase at the probe as it would be online, against the truth over .3-20 Hz
"""
import json
import math
from pathlib import Path

import numpy as np

import contact_candidates as C
import neumann_dynamics as D
import physics_plant3 as P3
import scenes3 as S3

HERE = Path(__file__).resolve().parent
TAU = np.array(S3.DYNAMIC['tau'])
QUIET = dict(contact_law='halfspace', probe='linear', slide='long', force=8., field_sigma=0., surface_amp=0., pose_error=0., mu=0., k_d=0., rho=0.,
             asymmetry=False, region='R3', subject='P04', accel=.1, speed=.0001, tool_mass=0., mass_error=0., tau_error=0., tau_a_true=.019, sensor_delay_ms=0., delay_error=0.)


def gain(alpha, tau, f):
    s = 2j*np.pi*np.asarray(f, float)[..., None]*tau
    return 1.+(s/(1.+s))@alpha


def on_grid(alpha):
    """Weights on the identification grid are not needed: rebuild the columns for the plant's own branches."""
    return np.r_[1., alpha]


def predict(r, alpha):
    t, d = r['t'], np.maximum(r['d'], 0.)
    u = d**D.EXP
    h = np.diff(t, prepend=t[0]); h[0] = h[1]
    du = np.diff(u, prepend=u[0])
    f = u.copy()
    for a, tau in zip(alpha, TAU):
        q = np.zeros_like(u); e = np.exp(-h/tau); g = tau/h*(1.-e)
        for k in range(1, len(u)):
            q[k] = e[k]*q[k-1]+du[k]*g[k]
        f = f+a*q
    return f


def against_neumann():
    records = [r for r in D.prepare() if D.usable(r)]
    out = {}
    for shape, tissue, phase in (('neural_white', 'white', 27), ('neural_grey', 'grey', 30)):
        alpha = np.array(S3.DYNAMIC['shapes'][shape][str(phase)]['alpha'])
        row = {}
        for kind in ('rate', 'relax', 'cyclic'):
            errors = []
            for r in [x for x in records if x['kind'] == kind and x['tissue'] == tissue]:
                c = D.in_contact(r); m = predict(r, alpha); s = float(r['f'][c]@m[c])/float(m[c]@m[c])
                errors.append(float(np.sqrt(np.mean((np.maximum(s*m[c], 0.)-r['f'][c])**2))/r['peak_force']))
            row[kind+'_error_pct'] = 100.*float(np.mean(errors))
        ph, ga = [], []
        for r in [x for x in records if x['kind'] == 'dma' and x['tissue'] == tissue]:
            f = np.array([x['f'] for x in r['rows']]); K = np.array([x['K'] for x in r['rows']]); m = gain(alpha, TAU, f)
            ph.append(np.degrees(np.angle(m))-np.degrees(np.angle(K))); e = np.log(np.abs(m))-np.log(np.abs(K)); ga.append(e-e.mean())
        row['dma_phase_rmse_deg'] = float(np.sqrt(np.mean(np.array(ph)**2))); row['dma_phase_error_at_20Hz_deg'] = float(np.mean(np.array(ph)[:, -1])); row['dma_gain_rmse_pct'] = 100.*float(np.sqrt(np.mean(np.array(ga)**2)))
        # Rate: the slope of force against speed, as measured
        rate = sorted([x for x in records if x['kind'] == 'rate' and x['tissue'] == tissue], key=lambda x: x['speed'])
        depth = .85*min(x['depth_peak'] for x in rate); v = [x['speed'] for x in rate]
        real = [float(np.interp(depth, x['d'][x['t'] <= x['t_peak']], x['f'][x['t'] <= x['t_peak']])) for x in rate]
        model = [float(np.interp(depth, x['d'][x['t'] <= x['t_peak']], predict(x, alpha)[x['t'] <= x['t_peak']])) for x in rate]
        row['rate_slope_measured'] = float(np.polyfit(np.log(v), np.log(real), 1)[0]); row['rate_slope_model'] = float(np.polyfit(np.log(v), np.log(model), 1)[0])
        # Relaxation against the median of the files
        pts = np.array([1., 10., 100.]); rm = []; rr = []
        for x in [x for x in records if x['kind'] == 'relax' and x['tissue'] == tissue]:
            t = x['t']-x['t_peak']; m = predict(x, alpha)
            rr.append([np.mean(x['f'][(t > p*.92) & (t < p*1.08)])/x['f'][t > .02][0] for p in pts]); rm.append([np.mean(m[(t > p*.92) & (t < p*1.08)])/m[t > .02][0] for p in pts])
        row['relaxation_measured_1_10_100s'] = np.round(np.median(rr, axis=0), 3).tolist(); row['relaxation_range_measured'] = [np.round(np.min(rr, axis=0), 3).tolist(), np.round(np.max(rr, axis=0), 3).tolist()]
        row['relaxation_model_1_10_100s'] = np.round(np.median(rm, axis=0), 3).tolist()
        # Cycles: peak decay and dissipation, real against model
        decay = []; loss = []; absent = []
        for x in [x for x in records if x['kind'] == 'cyclic' and x['tissue'] == tissue]:
            t, d = x['t'], x['d']; m = predict(x, alpha); c = D.in_contact(x); s = float(x['f'][c]@m[c])/float(m[c]@m[c]); m = np.maximum(s*m, 0.)
            T = 10.036
            for series, store in ((x['f'], 0), (m, 1)):
                peaks = [float(series[(t > x['t_peak']+i*T-1.5) & (t < x['t_peak']+i*T+1.5)].max()) for i in range(5)]
                up = (t >= x['t_peak']+T-T/2.) & (t <= x['t_peak']+T); dn = (t > x['t_peak']+T) & (t <= x['t_peak']+T+T/2.)
                eu = np.trapz(series[up], d[up]); ed = -np.trapz(series[dn], d[dn])
                decay.append((store, peaks[4]/peaks[0])); loss.append((store, 1.-ed/eu)); absent.append((store, float(np.mean(series[dn] < .02*peaks[0]))))
        for name, seq in (('cycle5_peak_ratio', decay), ('cycle2_dissipated', loss), ('unloading_no_force', absent)):
            row[name+'_measured'] = float(np.mean([v for k, v in seq if k == 0])); row[name+'_model'] = float(np.mean([v for k, v in seq if k == 1]))
        out[shape] = row
    return out


def against_greiner():
    """Human brain in uniaxial test, a different tissue and a different scale; only the form can be compared."""
    out = {}
    white = np.array(S3.DYNAMIC['shapes']['neural_white']['27']['alpha']); grey = np.array(S3.DYNAMIC['shapes']['neural_grey']['30']['alpha'])
    for name, a in (('neural_white', white), ('neural_grey', grey)):
        # the triangle of about 80 s period: its fundamental
        m = gain(a, TAU, np.array([1./80.]))
        step = lambda t: (1.+np.sum(a*np.exp(-t/TAU)))
        out[name] = dict(phase_at_period_80s_deg=float(np.degrees(np.angle(m))[0]), relaxation_1_10_100s=[round(step(t)/step(1.), 3) for t in (1., 10., 100.)])
    out['measured'] = dict(phase_loops_deg='5.5-10.1, area over 2 pi U', relaxation_compression_10_100s_of_1s=[[.899, .71], [.822, .531]], source='Greiner corona radiata and visual cortex, read by the scratch script in this analysis')
    return out


def plant_response():
    rows = []
    for shape, phase in (('neural_white', 21.), ('neural_grey', 28.), ('viscous', 14.)):
        sc = S3.cell_scene3('rest', 'relaxed', 'nominal', 0, T=30., base=493000, overrides=dict(QUIET, dynamics=shape, phase=phase))
        s0, _ = P3._init(sc, P3.PLATFORM)
        tick = np.arange(int(sc['T']/.008))*.008
        for f in (.5, 2., 8., 20.8):
            cmd = s0[2]+2.e-5*np.sin(2.*np.pi*f*tick)
            fine = dict(P3.PLATFORM, log_dt=.001, tau_a=sc['tau_a_true'])
            base = P3.replay(sc, np.full(len(tick), s0[2]), platform=fine)
            log = P3.replay(sc, cmd, platform=fine)
            q = log[:, 0] > 20.
            t = log[q, 0]
            a = np.c_[np.sin(2.*np.pi*f*t), np.cos(2.*np.pi*f*t)]
            ff = np.linalg.lstsq(a, log[q, 3]-base[q, 3], rcond=None)[0]; dd = np.linalg.lstsq(a, log[q, 6]-base[q, 6], rcond=None)[0]
            K = complex(*ff)/complex(*dd)
            analytic = gain(sc['alpha'], sc['tau'], np.array([f]))[0]
            rows.append(dict(shape=shape, f=f, phase_plant=float(np.degrees(np.angle(K))), phase_analytic=float(np.degrees(np.angle(analytic))), force_noise_free=True))
        # hold: force drift over 10 s at fixed position from the opening state, as the curve held by the scene start
        for wait in (2., 10., 60.):
            held = S3.cell_scene3('rest', 'relaxed', 'nominal', 0, T=30., base=493000, overrides=dict(QUIET, dynamics=shape, phase=phase, load_wait=wait))
            s1, _ = P3._init(held, P3.PLATFORM)
            hold = P3.replay(held, np.full(len(tick), s1[2]))
            rows.append(dict(shape=shape, f=0., wait_s=wait, force_at_start=float(hold[0, 3]), force_after_10s=float(np.interp(10., hold[:, 0], hold[:, 3])), force_after_30s=float(hold[-1, 3]), equilibrium=float(held['equilibrium']), k_static_over_curve=float(P3.contact(held, s1[13], float(held['xcmd'][0]))[1]/held['equilibrium']/P3.contact(dict(held, equilibrium=1.), P3.solve_depth(dict(held, equilibrium=1.), 8., float(held['xcmd'][0])), float(held['xcmd'][0]))[1])))
    return rows


def candidate_gap():
    """Each candidate at the phase the truth shows at the probe: how well does its assumed response follow the truth's over .3-20 Hz? Relative to the probe frequency, as the Method anchors there."""
    f = np.array([.3, 1., 3., 8., 14., 20.8])
    rows = []
    for shape in S3.SHAPES:
        for phase in S3.PHASES:
            a = np.array(S3.DYNAMIC['shapes'][shape][str(int(phase))]['alpha']); true = gain(a, TAU, f)
            for name in ('W', 'S', 'H'):
                phi = math.radians(phase)
                if name == 'H':
                    # rate-type phase lead through a constant multiplier
                    model = 1.+1j*2.*np.pi*f/(2.*np.pi*20.8)*math.tan(phi)
                else:
                    weights, tau = C.memory_from_phase(name, phi)
                    model = gain(np.asarray(weights), np.asarray(tau), f)
                g = np.abs(model)/abs(model[-1])/(np.abs(true)/abs(true[-1]))
                rows.append(dict(shape=shape, phase=phase, candidate=name, gain_error_at_0p3Hz_pct=100.*float(g[0]-1.), gain_rmse_pct=100.*float(np.sqrt(np.mean((g-1.)**2))), phase_error_at_3Hz_deg=float(np.degrees(np.angle(model[2]))-np.degrees(np.angle(true[2]))), phase_error_at_0p3Hz_deg=float(np.degrees(np.angle(model[0]))-np.degrees(np.angle(true[0])))))
    return rows


def main():
    out = dict(neumann=against_neumann(), greiner=against_greiner(), plant=plant_response(), candidates=candidate_gap())
    print('A. tables against the four protocols')
    for k, v in out['neumann'].items():
        print(' ', k, json.dumps({a: (round(b, 3) if isinstance(b, float) else b) for a, b in v.items()}))
    print('B. Greiner', json.dumps(out['greiner']))
    print('C. plant against analytic')
    for r in out['plant']:
        print(' ', {a: (round(b, 3) if isinstance(b, float) else b) for a, b in r.items()})
    print('D. candidate against truth, phase matched at the probe: gain error at .3 Hz %, gain rmse %, phase error at 3 Hz, at .3 Hz')
    for r in out['candidates']:
        print(f"  {r['shape']:13s} {r['phase']:4.0f} {r['candidate']} {r['gain_error_at_0p3Hz_pct']:+7.1f} {r['gain_rmse_pct']:6.1f} {r['phase_error_at_3Hz_deg']:+6.1f} {r['phase_error_at_0p3Hz_deg']:+6.1f}")
    (HERE/'results'/'dynamics3_checks.json').write_text(json.dumps(out, indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
