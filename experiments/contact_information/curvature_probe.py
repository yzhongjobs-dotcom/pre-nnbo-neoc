"""Curvature-information design with the original probe inside the candidate set.

Nominal servo/sensor only, no tissue parameters or simulation queries. Actual
nonlinear force cost still requires paired closed-loop measurement.
"""
import json
import math
from pathlib import Path
import numpy as np
from local_contact import periodic_encoder, predicted_force_rms, sine_reference_rms, residualize


def information(y, t, tw=.256):
    te = t[-1]+.008
    # Unit-response model; curvature parameter is H*a_force/k^2.
    m = np.arange(8)*.001
    af = -math.expm1(-2*math.pi*30*.001)
    w = af*(1-af)**(7-np.arange(8))
    frac = -np.expm1(-m/.019)/(-math.expm1(-.008/.019))
    v = (t[:, None]+m-te)/tw
    center = np.polynomial.polynomial.polyfit((np.r_[t, te]-te)/tw, y-y[-1], 3)
    u = y[:-1, None]+(y[1:]-y[:-1])[:, None]*frac-y[-1]-np.polynomial.polynomial.polyval(v, center)
    B = np.column_stack([(v**p)@w for p in range(5)])
    U = np.column_stack([u@w, (u*v)@w, (.5*u*u)@w])
    A = residualize(B, U)
    scale = np.linalg.norm(A, axis=0)
    _, s, vt = np.linalg.svd(A/scale, full_matrices=False)
    inverse = ((vt.T/s**2)@vt)/scale[:, None]/scale[None, :]
    functional = np.array([1., 0., -center[0]])
    return float(1/(functional@inverse@functional)), float(1/inverse[2, 2])


def summarize(wave=None):
    n = 32; t = np.arange(n)*.008
    if wave is None:
        w = 2*np.pi*20.8*.008; rho = math.exp(-.008/.019)
        G = (1-rho)*np.exp(-1j*w)/(np.exp(1j*w)-rho)
        values = [information(np.imag(G*np.exp(1j*(w*np.arange(n+1)+phase))), t)
                  for phase in np.linspace(0, 2*np.pi, 32, endpoint=False)]
    else:
        y = periodic_encoder(wave)
        values = [information(y[(np.arange(n+1)+phase)%n], t) for phase in range(n)]
    return dict(gain_information=float(min(v[0] for v in values)),
                curvature_information=float(min(v[1] for v in values)))


def design(count=192):
    rng = np.random.default_rng(73622)
    phase = np.arange(32)*2*np.pi/32
    hs = np.arange(1, 9)
    original = summarize()
    reference = sine_reference_rms()
    best = dict(wave=None, selected_kind='original_sine', **original)
    frontier = []
    for i in range(count):
        coef = rng.normal(size=(2, len(hs)))
        if i < len(hs):
            coef[:] = 0.; coef[0, i] = 1.
        wave = coef[0]@np.sin(hs[:, None]*phase)+coef[1]@np.cos(hs[:, None]*phase)
        scale = max(np.max(abs(wave)), np.sqrt(np.mean(wave**2))/.7071067811865476,
                    np.max(abs(wave-np.roll(wave, 1)))/(2*math.sin(np.pi*20.8*.008)),
                    predicted_force_rms(wave)/reference)
        wave /= scale
        info = summarize(wave)
        row = dict(index=i, **info)
        frontier.append(row)
        if info['gain_information'] >= original['gain_information'] and info['curvature_information'] > best['curvature_information']:
            best = dict(wave=wave.tolist(), selected_kind='periodic_fourier', **info,
                        peak=float(max(abs(wave))), rms=float(np.sqrt(np.mean(wave**2))),
                        predicted_force_rms=predicted_force_rms(wave))
    return dict(original=original, selected=best, all_candidates=frontier,
                protocol=dict(count=count, seed=73622, phases=32, gain_information_not_lower_than_original=True,
                              budget='original command peak/RMS/slew and nominal physical-force RMS',
                              plant_or_material_access=False, curvature_parameter='H*a_force/k^2',
                              phase_minimum_is_sampled_not_continuous_bound=True))


if __name__ == '__main__':
    target = Path(__file__).resolve().parent/'results/curvature_probe_design_v1.json'
    if target.exists():
        raise ValueError('preserve primary result')
    result = design()
    target.write_text(json.dumps(result, indent=2), encoding='utf8')
    print(json.dumps({k:v for k,v in result.items() if k != 'all_candidates'}, indent=2))
