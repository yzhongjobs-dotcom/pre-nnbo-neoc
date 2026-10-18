"""B0 footprint signal/timing and Hertz-scale audit, not a control benchmark.

Equations: B0 technical appendix B.2, D.1 and synthetic example G.2.
No measured sensor noise, validated tissue model or universal clinical loads
are assumed. Published load examples only probe the synthetic geometry's scale.
"""
from pathlib import Path
import json
import hashlib
import numpy as np
from scipy.integrate import quad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
CONFIG = dict(R_m=.03, L_m=.02, reference_distance_m=.02, a_m=.003,
              E0_Pa=30000., Ef_Pa=45000., speed_m_s=.008,
              provenance='B0 appendix G.2 geometry/materials; v=.008 m/s from its traction example',
              assumption='known a, old material and separated normal moment; optimistic information case',
              conditioning_noise_model='independent equal channel noise after force/moment normalization; no calibrated sensor SNR is implied',
              loads_N=[.1413716694115407, 1., 2.5, 3., 5., 10., 15.])


def signal(rho):
    # Second channel is the signed normal contribution to M_y, i.e. -m_1.
    w = np.sqrt(np.maximum(0., 1-rho*rho))
    return np.array([(np.arccos(rho)-rho*w)/np.pi, -2*w**3/(3*np.pi)])


def jacobian(rho):
    w = np.sqrt(1-rho*rho)
    # Columns: normalized material contrast A and normalized boundary position.
    return np.column_stack([signal(rho), [-2*w/np.pi, 2*rho*w/np.pi]])


def main():
    c = CONFIG
    a, R, L = c['a_m'], c['R_m'], c['L_m']
    full_force_contrast = np.pi*L*(c['Ef_Pa']-c['E0_Pa'])*a*a/(4*R)
    rows = []
    for rho in [.99, .95, .9, .75, .5, .25, 0.]:
        s = signal(rho)
        J = jacobian(rho)
        inverse = np.linalg.inv(J)
        rows.append(dict(rho=rho, entered_distance_fraction=1-rho,
             force_fraction=float(s[0]), normal_moment_fraction=float(s[1]),
             synthetic_delta_force_N=float(full_force_contrast*s[0]),
             synthetic_normal_moment_Nm=float(a*full_force_contrast*s[1]),
             time_to_boundary_center_s=float(rho*a/c['speed_m_s']),
             known_boundary_contrast_noise_amplification=float(1/np.linalg.norm(s)),
             unknown_boundary_contrast_noise_amplification=float(np.linalg.norm(inverse[0])),
             normalized_jacobian_condition=float(np.linalg.cond(J))))

    geometries = []
    for force in c['loads_N']:
        ar = np.sqrt(4*R*force/(np.pi*L*c['E0_Pa']))
        geometries.append(dict(load_N=force, required_uniform_halfwidth_m=float(ar),
                 a_over_R=float(ar/R), fullwidth_over_L=float(2*ar/L),
                 a_less_than_reference_distance=bool(ar<c['reference_distance_m'])))

    # Independent integral evaluation across the entire partial-coverage interval.
    quadrature_error = 0.
    for rho in np.linspace(-.99, .99, 81):
        q0 = quad(lambda x: np.sqrt(1-x*x), rho, 1, epsabs=1e-12)[0]*2/np.pi
        q1 = -quad(lambda x: x*np.sqrt(1-x*x), rho, 1, epsabs=1e-12)[0]*2/np.pi
        quadrature_error = max(quadrature_error, float(np.max(np.abs(signal(rho)-[q0,q1]))))
    eps = np.geomspace(1e-5, 1e-3, 60)
    slopes = [float(np.polyfit(np.log(eps), np.log(np.abs(signal(1-eps)[j])), 1)[0]) for j in range(2)]
    source = ROOT.parents[1]/'研究'/'技术附录.txt'
    result = dict(config=c, source_appendix_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  partial_coverage=rows, uniform_geometry_load_screen=geometries,
                  full_contrast_force_N=float(full_force_contrast),
                  ideal_first_entry_to_center_s=a/c['speed_m_s'],
                  synthetic_uniform_domain_max_force_N=float(np.pi*L*c['E0_Pa']*c['reference_distance_m']**2/(4*R)),
                  max_closed_form_quadrature_error=quadrature_error,
                  entry_loglog_signal_slopes=slopes,
                  causal_bound='No local observation depends on material outside the union of already-contacted footprints, absent a spatial prior.',
                  timing_bound='For adjustment before center arrival: available = rho*a/v - estimation/compute/communication delay; a required gain change needs at least max_j(abs(delta_gain_j)/rate_j), possibly more under tank constraints.',
                  interpretation='Rechecks/quantifies conditions already acknowledged by B0; not new identifiability theory, an MHE run, a clinical force limit, or proof that preview is necessary.')
    assert quadrature_error<1e-9
    assert all(abs(s-1.5)<.002 for s in slopes)
    (ROOT/'results.json').write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf8')

    rho = np.linspace(0, .995, 500)
    sig = signal(rho)
    fig, ax = plt.subplots(1,2,figsize=(10,3.8),constrained_layout=True)
    ax[0].plot(1-rho,sig[0],label='Force contrast')
    ax[0].plot(1-rho,-sig[1],label='Normal moment / halfwidth')
    ax[0].set(xlabel='Distance since first entry / halfwidth', ylabel='Signal / full-footprint force contrast')
    ax[0].legend()
    known = [1/np.linalg.norm(signal(r)) for r in rho]
    unknown = [np.linalg.norm(np.linalg.inv(jacobian(r))[0]) for r in rho]
    ax[1].semilogy(1-rho,known,label='Known boundary')
    ax[1].semilogy(1-rho,unknown,label='Jointly unknown boundary')
    ax[1].set(xlabel='Distance since first entry / halfwidth', ylabel='Local contrast-noise amplification')
    ax[1].legend()
    fig.savefig(ROOT/'footprint_information.png',dpi=180)
    print(json.dumps(dict(quadrature_error=quadrature_error,entry_slopes=slopes,
                         ideal_lead_s=a/c['speed_m_s'],coverage_at_rho_09=rows[2],geometry=geometries[3:5]),ensure_ascii=False))


if __name__=='__main__':
    main()
