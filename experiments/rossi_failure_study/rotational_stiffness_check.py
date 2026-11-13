import os, sys, math
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'): os.environ[k]='1'
sys.dont_write_bytecode=True; sys.path.insert(0, os.getcwd())
import numpy as np
import corrected_comparison as CC, overlay_compare as O, memory_segment_eval as SE
P3=SE.P3; O.setup()
print('case group law slide L(mm) R(mm) | x(mm) theta0(deg) depth(mm) F(N) k_n(N/m) | K_rot/(k_n L^2)  [Winkler flat: 1/12=.0833] | dF for +-0.2deg (N) | M0 (mN m) | moment per 0.2deg (mN m)')
for idx in (0,3,12,36,39,24):
    case=CC.cases()[idx]; sc=SE.MC.scene_for(tuple(case['cell']),CC.BASE,case['variant'])
    g=sc['geometry']; L=g['length']
    xs=np.linspace(sc['xcmd'][0], sc['xcmd'][-1], 7)[1:-1]
    for x in xs[::2]:
        try:
            th0=float(P3._theta(x, P3._constants(sc, P3.PLATFORM)) ) if hasattr(P3,'PLATFORM') else 0.
        except Exception:
            th0=0.
        for th in (th0, th0+math.radians(1.5)):
            d=P3.solve_depth(sc, sc['FD'], x, th)
            f0,k0,c0,_,_=P3.contact(sc, d, x, th)
            e=math.radians(.1)
            fp,kp,cp,_,_=P3.contact(sc, d, x, th+e); fm,km,cm,_,_=P3.contact(sc, d, x, th-e)
            krot=-(fp*cp-fm*cm)/(2*e)
            e2=math.radians(.2)
            f2p=P3.contact(sc, d, x, th+e2)[0]; f2m=P3.contact(sc, d, x, th-e2)[0]
            print(f"{idx:2d} {case['group']:8s} {sc['contact_law']:10s} {sc['slide']:5s} {1e3*L:5.1f} {1e3*g['radius'] if math.isfinite(g['radius']) else float('inf'):6.1f} | {1e3*x:6.1f} {math.degrees(th):+5.2f} {1e3*d:5.2f} {f0:5.2f} {k0:7.0f} | {krot/(k0*L*L):+.4f} | {f2p-f0:+.3f} {f2m-f0:+.3f} | {1e3*f0*c0:+7.2f} | {1e3*krot*e2:6.2f}")
