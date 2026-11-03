"""Does an added frequency provide the missing slow contact-response information?

Pure design calculation. No plant, task trajectory, hidden state or tuned data.
Fixed relaxation times are favorable conditional hypotheses, not fitted truth.
"""
import json
import math
from pathlib import Path
import numpy as np
from memory_contact import features
from local_contact import residualize, sine_reference_rms

HERE = Path(__file__).resolve().parent
TAU = .019
LOW = 1/(2*math.pi*TAU)


def bands(kind):
    if kind == 'single':
        return [(20.8, 1.)]
    ratio = sine_reference_rms()/sine_reference_rms(LOW)
    if kind == 'lower':
        return [(LOW, ratio)]
    tones = [(LOW, ratio/math.sqrt(2)), (20.8, 1/math.sqrt(2))]
    # Analytic bounds on command peak, first and second sample differences,
    # and infinite-time servo-output RMS. All no larger than original probe.
    scale = max(sum(a for f,a in tones),
                sum(a*abs(math.sin(math.pi*f*.008)) for f,a in tones)/math.sin(math.pi*20.8*.008),
                sum(a*math.sin(math.pi*f*.008)**2 for f,a in tones)/math.sin(math.pi*20.8*.008)**2,
                math.sqrt(sum((a*sine_reference_rms(f))**2 for f,a in tones))/sine_reference_rms())
    return [(f, a/scale) for f,a in tones]


def information(kind, tau_m, nw):
    tones = bands(kind)
    phases = [(p,) for p in np.arange(8)*2*np.pi/8] if len(tones) == 1 else [
        (p,q) for p in np.arange(4)*np.pi/2 for q in np.arange(4)*np.pi/2]
    rows = []
    t = np.arange(nw+1)*.008
    rho = math.exp(-.008/TAU)
    for phase in phases:
        y = np.zeros(nw+1)
        for (frequency, amplitude), angle in zip(tones, phase):
            w = 2*np.pi*frequency*.008
            G = (1-rho)*np.exp(-1j*w)/(np.exp(1j*w)-rho)
            y += .16740639598337909*amplitude*np.imag(G*np.exp(1j*(w*np.arange(nw+1)+angle)))
        B,U,_ = features(t[:-1], y[:-1], y[1:], t[-1], y[-1], t[-1], TAU, tau_m, 0.)
        X = residualize(B,U)
        norm = np.linalg.norm(X,axis=0)
        _,s,vt = np.linalg.svd(X/norm,full_matrices=False)
        cov = .04**2*((vt.T/s**2)@vt)/norm[:,None]/norm[None,:]
        std = {}
        for f in (0.,3.,20.8):
            q = (2*np.pi*f*tau_m)**2
            c = np.array([1.,0.,q/(1+q)])
            std[str(f)] = math.sqrt(float(c@cov@c))
        rows.append(std)
    return dict(kind=kind, tau_m=tau_m, window=nw*.008, tones=tones,
                conditional_gain_sd_P50={f:float(np.median([r[f] for r in rows])) for f in rows[0]},
                conditional_gain_sd_max={f:max(r[f] for r in rows) for f in rows[0]})


def main():
    target = HERE/'results/response_band_information_v1.json'
    if target.exists():
        raise ValueError('preserve result')
    result = dict(protocol=dict(question='can one extra frequency make control-band response identifiable under the existing moving-surface nuisance?',
        physically_selected_low_frequency=LOW, no_truth=True, no_plant=True,
        target_force_proxy='nominal unit contact gain; original amplitude .1674 N and noise .04 N',
        sd_units='relative to the gain used for command conversion; NOT relative error of an unknown ke',
        covariance='fixed design, IID .04 N force innovation, relaxation time known; more favorable than unknown tau',
        nuisance='B4 plus unknown initial memory exponential; endpoint relaxed gain+slope, constant memory modulus',
        extra_window='information/time tradeoff only; assumes local model persists for the whole window',
        no_closed_loop_benefit_claim=True), rows=[])
    for nw in (32,64):
        for tau_m in (.04,.2,.6):
            for kind in ('single','lower','dual'):
                row = information(kind,tau_m,nw)
                result['rows'].append(row)
                print(nw,tau_m,kind,row['conditional_gain_sd_P50'])
    target.write_text(json.dumps(result,indent=2),encoding='utf8')


if __name__ == '__main__':
    main()
