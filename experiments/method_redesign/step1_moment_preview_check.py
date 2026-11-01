import sys, numpy as np
sys.path.insert(0, r'E:/github_project/scan-impedance-update/experiments/cobot_scan_benchmark')
import plant as P
from scipy.optimize import brentq
def pvec(E0, ratio=20.):
    p = np.zeros(P.NP); p[0]=E0; p[1]=0.2; p[6]=1; p[7]=0.05; p[8]=0.008; p[9]=0.0; p[10]=E0*ratio; p[11]=0.001; p[13]=1.0; p[14]=1.0
    return p
for E0 in (1e4, 4e4):
    p = pvec(E0)
    print('E0=%.0e  (constant 4 N force; depth solved per position)' % E0)
    print(' dist_to_edge[mm]  depth[mm]  halfwidth[mm]  m1[N.m]  m1/Fn[mm]  k/k_base  Ft_moment_term[N]')
    k0=None
    for dist in (-20,-15,-10,-7,-4,-2,0,3,6):
        x = P.X0-0.008+dist*1e-3
        d = brentq(lambda dd: P.footprint(x, dd, p)[0]-4.0, 1e-7, 0.03)
        fn, m1, mu, a = P.footprint(x, d, p)
        k = P.k_eff(x, d, p); k0 = k0 or k
        print('  %6.0f %10.2f %10.2f %11.2e %9.3f %8.2f %10.3f' % (dist, d*1e3, a*1e3, m1, m1/fn*1e3, k/k0, m1/P.R))
