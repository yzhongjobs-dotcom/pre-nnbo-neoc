"""Static response against measured curves, independent of every candidate: Remus (human abdomen) and Beber Fig. 2 (phantom, extracted from the figure's vector data)."""
import json
import math
from pathlib import Path

import numpy as np

import physics_plant3 as P
import scenes3 as S

HERE = Path(__file__).resolve().parent
# Beber Fig. 2, loading branch: penetration mm, force N, sphere r = 20 mm, chest phantom
BEBER_D = np.arange(1., 15.)
BEBER_F = np.array([.52, 1.08, 1.78, 2.54, 3.37, 4.27, 5.25, 6.25, 7.27, 8.31, 9.42, 10.47, 11.54, 12.61])


def main():
    out = {}
    # round trip: the transfer on the hemisphere itself returns the measured curve
    hemi = []
    d = np.geomspace(1.e-4, .025, 200)
    for name, poly in S.REMUS.items():
        f = np.polyval(poly, d[d <= .025]*1000.)
        k = np.polyval(np.polyder(poly), d*1000.)*1000.
        fi = np.r_[0., np.cumsum(.5*(k[1:]+k[:-1])*np.diff(d))]+f[0]
        hemi.append(float(np.max(np.abs(fi[d > .002]/f[d > .002]-1.))))
    out['hemisphere_integral_max_error'] = max(hemi)
    print('polynomial tangent integrates back to the curve, max error', round(max(hemi), 4))
    rows = []
    print('\nprobe   gamma region | depth mm at 4 8 12 N | tangent N/mm at 4 8 12 N | exponent over 4-12 N')
    for probe, geo in S.PROBES.items():
        for gamma in (.5, .75, 1.):
            for region in ('R3', 'R4', 'R1', 'A5'):
                depths = np.geomspace(1.e-5, .06, 800)
                F, k = P.probe_curve(np.array(S.REMUS[region]), geo, gamma, depths)
                dd = np.interp([4., 8., 12.], F, depths); kk = np.interp([4., 8., 12.], F, k)
                n = math.log(3.)/math.log(dd[2]/dd[0])
                rows.append(dict(probe=probe, gamma=gamma, region=region, depth_mm=(dd*1e3).tolist(), k_N_mm=(kk/1e3).tolist(), exponent=n))
                print(f"{probe:6s} {gamma:.2f} {region} | {dd[0]*1e3:5.1f} {dd[1]*1e3:5.1f} {dd[2]*1e3:5.1f} | {kk[0]/1e3:.2f} {kk[1]/1e3:.2f} {kk[2]/1e3:.2f} | {n:.2f}")
    out['probe_curves'] = rows
    # plant reproduces the table when uniform and aligned
    errs = []
    for law in ('foundation', 'halfspace'):
        for probe in S.PROBES:
            sc = S.cell_scene3('rest', 'relaxed', 'nominal', 0, base=450000, overrides=dict(contact_law=law, probe=probe, slide='long', field_sigma=0., pose_error=0., surface_amp=0., region='R1', force=8.))
            depths = np.geomspace(1.e-3, .03, 40)
            F, k = P.probe_curve(sc['remus_poly'], sc['geometry'], sc['transfer_gamma'], np.geomspace(1.e-5, .06, 2000))
            Fi = np.interp(depths, np.geomspace(1.e-5, .06, 2000), F)
            Fp = np.array([P.contact(sc, x)[0] for x in depths])
            errs.append(float(np.max(np.abs(Fp/Fi-1.))))
    out['plant_vs_curve_max_error'] = max(errs)
    print('\nplant force against the curve when uniform and aligned, max error', round(max(errs), 4))
    # size scaling check: the same transfer from the r = 10 mm curves to Beber's r = 20 mm sphere. Different material, so a plausibility range only.
    print('\nsphere r = 20 mm, 15 mm penetration, from the abdomen curves: force N; Beber phantom measured 13.2 N')
    for region in ('R1', 'R2', 'R3', 'R4', 'R5'):
        d = np.geomspace(1.e-5, .015, 500); A = np.pi*(2.*.02*d-d*d)
        s = d/np.sqrt(A); dh = np.geomspace(1.e-5, .05, 4000); sh = dh/np.sqrt(P._hemisphere_area(dh)); delta = np.maximum(np.interp(s, sh, dh), .005)
        k = np.polyval(np.polyder(S.REMUS[region]), delta*1000.)*1000.*(A/P._hemisphere_area(delta))**.5
        F = np.r_[0., np.cumsum(.5*(k[1:]+k[:-1])*np.diff(d))]
        print(region, 'F', round(F[-1], 1), 'exponent', round(np.polyfit(np.log(d[d > .002]), np.log(F[d > .002]), 1)[0], 2))
    print('Beber exponent', round(np.polyfit(np.log(BEBER_D), np.log(BEBER_F), 1)[0], 2))
    (HERE/'results'/'curve3_checks.json').write_text(json.dumps(out, indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
