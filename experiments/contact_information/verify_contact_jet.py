"""Checks of changed contact algebra, endpoint functional and filter ordering."""
import json
from pathlib import Path
import numpy as np
import local_contact
from contact_jet import features, fit, predict, MODELS


def main():
    t = np.arange(80)*.008
    y = .002+.0004*t+.00003*np.sin(2*np.pi*20.8*t)+.00001*np.sin(2*np.pi*9.1*t)
    n = 64
    measured = np.c_[t, t*.02, y, np.full(len(t), 4.), np.zeros((len(t), 3))]
    records = []
    for model in MODELS:
        te, ye, tw = t[n-1], y[n-1], .256
        j = np.arange(n-33, n-1)
        center = None
        if model == 'centered_curvature':
            ix = np.arange(n-33, n)
            center = np.polynomial.polynomial.polyfit((t[ix]-te)/tw, y[ix]-ye, 3)
        B, U = features(t[j], y[j], y[j+1], te, ye, tw, model=model, center=center)
        theta = np.array([12000., 1200.] + ([18000.] if U.shape[1] == 3 else []))
        beta = np.zeros(B.shape[1]); beta[0] = 4.; beta[1] = .04
        dec = np.exp(-2*np.pi*30*.008)
        # Generate both fitting and unseen intervals from the stated force law.
        jj = np.arange(n-33, n+4)
        BB, UU = features(t[jj], y[jj], y[jj+1], te, ye, tw, model=model, center=center)
        zz = BB@beta+UU@theta
        for q, z in zip(jj, zz):
            measured[q+1, 3] = dec*measured[q, 3]+z
        r = fit(measured[:n], model=model)
        np.testing.assert_allclose(r['coefficients'], theta, rtol=1e-7, atol=1e-4)
        for h in (1, 4):
            fp, _ = predict(r, measured, n, h, .019, model)
            np.testing.assert_allclose(fp, measured[n-1+h, 3], atol=1e-9)
        # Independent fixed-design influence-row covariance of the endpoint gain.
        X = np.column_stack([B, U]); norm = np.linalg.norm(X, axis=0)
        inverse = np.linalg.pinv(X/norm, rcond=1e-12)/norm[:, None]
        functional = np.r_[np.zeros(B.shape[1]), r['functional']]
        influence = functional@inverse
        np.testing.assert_allclose(r['sigma']**2, .04**2*(influence@influence), rtol=1e-7)
        if model == 'lean':
            old = local_contact.fit(measured[:n], nw=32, rate=False)
            np.testing.assert_allclose(r['coefficients'], old['coefficients'], rtol=1e-10)
        records.append(dict(model=model, recovered=r['coefficients'].tolist(),
                            k=r['k'], endpoint_covariance_checked=True, held_out_prediction_checked=True))
    # A nonlinear feature is filter(u^2), never filter(u)^2.
    B, U = features(np.array([0.]), np.array([0.]), np.array([.0001]), .008, 0., .256, model='normal_curvature')
    mass = B[0, 0]
    assert abs(U[0, 2]-.5*(U[0, 0]/mass)**2/.001*mass) > 1e-8
    target = Path(__file__).resolve().parent/'results/contact_jet_verification.json'
    target.write_text(json.dumps(dict(checks=records, nonlinear_filter_order=True), indent=2), encoding='utf8')
    print('contact jet recovery, predictive filtering, endpoint covariance and lean equivalence passed')


if __name__ == '__main__':
    main()
