"""Re-select estimator parameters after two corrections, then repeat the paired closed-loop evaluation on fresh seeds.

1. C_neutral estimator tuning ran with true gains, so estimate quality could not influence the selection.
2. The band-pass estimators lacked the contact test window_fit applies; approach motion entered their memory.
Public control parameters are unchanged. R2 and Projection under C_projection keep their earlier parameters.
"""
import json
from pathlib import Path

import sweeps as S

HERE = Path(__file__).resolve().parent


def main():
    path = HERE/'results'/'tuning.json'
    t = json.loads(path.read_text(encoding='utf-8'))
    pub = t['public']
    for name in ('Projection', 'BandPass', 'R2'):
        t['estimators']['C_neutral'][name] = S.tune_estimator(name, pub['C_neutral']['best']['params'], False)
        print('C_neutral', name, t['estimators']['C_neutral'][name]['best'], flush=True)
    t['estimators']['C_projection']['BandPass'] = S.tune_estimator('BandPass', pub['C_projection']['best']['params'], False)
    print('C_projection BandPass', t['estimators']['C_projection']['BandPass']['best'], flush=True)
    t.pop('expanded_boundary', None)
    path.write_text(json.dumps(t, indent=1, allow_nan=True), encoding='utf-8')
    S.run_main(t, 15, base=386000)


if __name__ == '__main__':
    main()
