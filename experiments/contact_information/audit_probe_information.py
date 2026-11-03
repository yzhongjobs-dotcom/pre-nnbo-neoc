"""Does apparent gain information come from the imposed probe frequency?

Reuses frozen measurements. No plant truth, no new waveforms, no adoption.
Projection on deterministic probe instruments is a diagnostic of frequency
support, not a causal decomposition of every encoder displacement.
"""
import hashlib
import json
import math
import argparse
from pathlib import Path
import numpy as np
from local_contact import design, features, residualize, fit

HERE = Path(__file__).resolve().parent
OUT = HERE/'results'


def instrument_fit(measured, nw=32, tau=.019, fc=30., noise=.04, recorded_probe=None):
    data = design(measured, nw=nw, tau=tau, fc=fc)
    if data is None:
        return None
    B, U, z = data
    n = len(measured); j = np.arange(max(1, n-1-nw), n-1)
    j = j[(measured[j, 3] > .5) & (measured[j+1, 3] > .5)]
    te, ye, tw = measured[-1, 0], measured[-1, 2], nw*.008
    phases = 2*np.pi*20.8*measured[:, 0]
    columns = []
    for c, fun in enumerate((np.sin, np.cos)):
        signal = fun(phases) if recorded_probe is None else recorded_probe[:n, c]
        _, q = features(measured[j, 0], signal[j], signal[j+1], te, 0., tw, tau, fc)
        columns.extend([q[:, 0], q[:, 1]])
    Q = residualize(B, np.column_stack(columns))
    q, s, _ = np.linalg.svd(Q, full_matrices=False)
    q = q[:, s > s[0]*1e-10]
    X = residualize(B, U[:, :2]); target = residualize(B, z)
    xp = q@(q.T@X)
    xp_norm = np.linalg.norm(xp, axis=0)
    if min(xp_norm) < 1e-13:
        return None
    _, ss, vt = np.linalg.svd(xp/xp_norm, full_matrices=False)
    if ss[-1] < ss[0]*1e-8:
        return None
    theta = np.linalg.lstsq(xp/xp_norm, q@(q.T@target), rcond=1e-8)[0]/xp_norm
    inverse = ((vt.T/ss**2)@vt)/xp_norm[:, None]/xp_norm[None, :]
    covariance = noise**2*inverse
    full_inverse = np.linalg.inv(X.T@X)
    information_full = 1/(noise**2*full_inverse[0, 0])
    information_probe = 1/covariance[0, 0]
    residual = target-X@theta
    j_stat = float(np.linalg.norm(q.T@residual)**2/noise**2)
    return dict(k=float(theta[0]), sigma=float(np.sqrt(covariance[0, 0])),
                coefficients=theta, beta=np.linalg.lstsq(B, z-U[:, :2]@theta, rcond=1e-10)[0],
                information_full=information_full, information_probe=information_probe,
                probe_fraction=information_probe/information_full,
                overidentification_reference=j_stat, overidentification_df=q.shape[1]-2,
                te=float(te), ye=float(ye), tw=tw,
                # The same-model IID reference is not a closed-loop guarantee.
                chi2=float(residual@residual/(max(len(z)-7, 1)*noise**2)))


def recover_probe(measured):
    """Replay the SAME causal policy; reconstruct its own known probe command.

    The final logged row holds the newly issued command; at observation time
    that cell held the previously sent command. Restore that provisional cell
    while calling the policy. No physical state or score labels are provided.
    """
    from run_joint import JointPolicy, common, CE
    pol = JointPolicy(rate=False, single_model=False, **CE.FINAL)
    pol.reset(dict(common.P.PLATFORMS['PC125']))
    log = np.full((len(measured), 16), np.nan)
    log[:, [0, 1, 2, 5, 11, 12, 13]] = measured
    command = np.zeros((len(measured), 2)); encoder = np.zeros_like(command)
    tau = np.empty(len(measured)); delay = np.empty(len(measured), int)
    errors = []; last_cos = 0.
    for j in range(len(measured)):
        log[j, 11] = measured[j-1, 4] if j else measured[0, 2]
        previous = pol.x_prev
        issued = pol(common.P.Obs(log, j+1))
        log[j, 11] = measured[j, 4]
        errors.append(float(issued-measured[j, 4]))
        active = pol.params['excitation'] and pol.mode == 0 and measured[j, 3] >= .3
        k = pol.k_good if np.isfinite(pol.k_good) else pol.k_seek
        angle = 2*np.pi*pol.params['f_exc']*(measured[j, 0]+pol.Tc)
        cosine = pol.params['a_force']*pol.exc_gain*math.cos(angle) if active else 0.
        increment = np.array([(pol.x_prev-previous)/k if active else 0., (cosine-last_cos)/k if active else 0.])
        command[j] = (command[j-1] if j else np.zeros(2))+increment
        last_cos = cosine; tau[j] = pol.tau; delay[j] = pol.ndel
        if j:
            rho = math.exp(-pol.Tc/tau[j-1]); source = j-1-delay[j-1]
            forcing = command[source] if source >= 0 else np.zeros(2)
            encoder[j] = rho*encoder[j-1]+(1-rho)*forcing
    maximum = float(np.max(np.abs(errors)))
    if maximum > 1e-10:
        raise AssertionError(f'causal replay mismatch {maximum}')
    return encoder, dict(max_command_error=maximum, controller_parameters=pol.params,
                         trace='own causal probe command, nominal/identified servo only; not actual plant counterfactual')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--label', default='probe_information_audit_v2')
    parser.add_argument('--recorded-probe', action='store_true')
    args = parser.parse_args()
    target = OUT/(args.label+'.json')
    if target.exists():
        raise ValueError('preserve primary output')
    rows, summary, inputs, traces = [], [], {}, {}
    config = json.loads((HERE.parent/'motion_aware_scan/tuned_v2/Proposed-IMMb.json').read_text())['best']['params']
    for path in sorted(OUT.glob('lean_input_*.npz')):
        inputs[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
        data = np.load(path); m = data['measured']; taus = data['causal_tau']
        family = json.loads(str(data['metadata']))['family']
        recorded = None
        if args.recorded_probe:
            recorded, traces[family] = recover_probe(m)
        for end in range(192, len(m)-1, 8):
            tau = float(taus[end-1])
            r = instrument_fit(m[:end], tau=tau, recorded_probe=recorded)
            old = fit(m[:end], nw=32, tau=tau, rate=False, return_aux=True)
            if r is None or old is None:
                continue
            accepted = bool(200 <= old['k'] <= 3e5 and old['start_gain'] > 0 and
                            old['risk'] <= math.log(1.2) and old['chi2'] <= config['chi2_max'])
            B, U = features(m[end-1:end, 0], m[end-1:end, 2], m[end:end+1, 2],
                            r['te'], r['ye'], r['tw'], tau, 30.)
            decay = math.exp(-2*np.pi*30*.008)
            predictions = {}
            for name, result in [('least_squares', old), ('probe_moments', r)]:
                predicted = decay*m[end-1, 3]+float((B@result['beta']+U[:, :2]@result['coefficients'])[0])
                predictions[name] = predicted-m[end, 3]
            rows.append(dict(family=family, end=end, time=float(m[end-1, 0]),
                accepted_by_lean=accepted, k_lean=old['k'], sigma_lean=old['sigma'],
                k_probe=r['k'], sigma_probe=r['sigma'], probe_fraction=r['probe_fraction'],
                overidentification_reference=r['overidentification_reference'],
                overidentification_df=r['overidentification_df'], chi2_lean=old['chi2'],
                conditional_force_error=predictions))
        rr = [r for r in rows if r['family'] == family]
        accepted = [r for r in rr if r['accepted_by_lean']]
        summary.append(dict(family=family, windows=len(rr), accepted=len(accepted),
            probe_information_fraction_P10_P50_P90=np.percentile([r['probe_fraction'] for r in rr], [10, 50, 90]).tolist(),
            accepted_probe_fraction_P10_P50_P90=np.percentile([r['probe_fraction'] for r in accepted], [10, 50, 90]).tolist(),
            moment_reference_P50_P90=np.percentile([r['overidentification_reference'] for r in rr], [50, 90]).tolist(),
            gain_disagreement_relative_P50_P90=np.percentile([abs(r['k_probe']-r['k_lean'])/max(abs(r['k_lean']), 1.) for r in rr], [50, 90]).tolist(),
            conditional_force_rmse={name:float(np.sqrt(np.mean([r['conditional_force_error'][name]**2 for r in rr])))
                                    for name in ('least_squares', 'probe_moments')}))
        print(json.dumps(summary[-1]), flush=True)
    target.write_text(json.dumps(dict(protocol=dict(inputs_sha256=inputs,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        observations='seven measured channels and causal servo ID only; no scoring labels',
        replayed_traces=traces, chi2_gate=config['chi2_max'],
        baseline='preserved lean candidate', stride=8, overlapping_windows=True,
        frequency=20.8, probe_basis=('actual causal probe input through known servo plus quadrature and linear envelopes' if args.recorded_probe else
                                    'sine/cosine and linear envelopes; same encoder/filter domain'),
        moments_are_conditional_reference_not_calibrated_tests=True,
        held_out_prediction='future encoder given, next force hidden', plant_changed=False), summary=summary, rows=rows), indent=2), encoding='utf8')


if __name__ == '__main__':
    main()
