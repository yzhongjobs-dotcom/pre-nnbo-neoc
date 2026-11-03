"""Print the bank internals on one scene, against truth, for diagnosis only."""
import json
import math
import sys
import numpy as np
from bank_eval import cohort, FINAL, BankContact, MeasuredOnly, common


class Tr(BankContact):
    def _stiffness_est(self, o, n):
        r = super()._stiffness_est(o, n)
        if self.ready:
            self.rec.append((float(o.log[n-1, 0]), self.slow.k, self.slow.sig, self.fast.k, self.fast.sig,
                             float(self.k_good), float(self.kbar), float(self.exceed), self.fast.nis))
        return r

    def reset(self, p):
        super().reset(p); self.rec = []


def main():
    fam, idx = sys.argv[1], int(sys.argv[2]); kw = json.loads(sys.argv[3]) if len(sys.argv) > 3 else {}
    scene = [s for s in cohort('diag', 50) if s[0] == fam][idx]
    _, sid, p, T = scene
    pol = Tr(single_model=False, feedback='projection', **kw, **FINAL)
    log, _ = common.run(MeasuredOnly(pol), p, T, sid, platform=dict(common.P.PLATFORMS['PC125']))
    t = log[:, 0]; kt = log[:, 8]; f = log[:, 3]
    a = np.asarray(pol.rec)
    print(sid, 'rmse', round(common.task_metrics(log, 20.8)['rmse'], 3), 'peak', round(f.max(), 1), 'reseeds', pol.reseeds)
    step = max(int(len(a)/40), 1)
    print('   t   k_true  slow    sig   fast    sig   k_good  k_bar  gap   unc  exc  f_true')
    for row in a[::step]:
        i = int(np.argmin(abs(t-row[0])))
        gap = abs(math.log(row[3]/row[1])); unc = math.hypot(row[4]/row[3], row[2]/row[1])
        print(f'{row[0]:5.2f} {kt[i]:7.0f} {row[1]:7.0f} {row[2]:6.0f} {row[3]:7.0f} {row[4]:6.0f} {row[5]:7.0f} {row[6]:7.0f} {gap:5.2f} {unc:5.2f} {row[7]:3.0f} {f[i]:6.2f}')


if __name__ == '__main__':
    main()
