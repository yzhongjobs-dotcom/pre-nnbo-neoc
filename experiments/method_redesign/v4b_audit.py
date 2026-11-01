"""Post-selection mechanism audit. Truth is read only AFTER closed-loop runs.

Stores aggregate JSON only. Does not tune parameters or read final seeds.
"""
import json
import math
import numpy as np
import v4b_evaluate as E
from v4b_architecture import V4b, P, CE


def stats(values):
    x = np.asarray(values, float)
    x = x[np.isfinite(x)]
    if not len(x):
        return dict(n=0, median=None, mean=None, p95=None)
    return dict(n=len(x), median=float(np.median(x)), mean=float(x.mean()),
                p95=float(np.percentile(x, 95)))


class Audit(V4b):
    def reset(self, platform):
        super().reset(platform)
        self.audit = []
        self.loss_entries = []

    def __call__(self, o):
        oldmode = self.c1.mode
        command = super().__call__(o)
        if oldmode == 0 and self.c1.mode == 1:
            self.loss_entries.append(dict(j=o.n-2, residual_velocity_used=self.c1.motion_before_loss,
                last_IMM_surface_velocity=-float(getattr(self.observer,'v_cmd_motion',0.)),
                recovery_command_velocity=float(np.clip(-self.c1.motion_before_loss+self.c1.approach,-.04,.06))))
        c = self.c1
        self.audit.append(dict(j=o.n-2, k=c.kuse, kg=self.observer.k_good,
                               velocity=c.th[2], mode=c.mode, level=c.th[1],
                               reference=c.g[5], offset=self.ff_offset))
        return command


def main():
    rows = json.loads((E.OUT/'v4b_episode_summary.json').read_text(encoding='utf-8'))
    selected = []
    scenes = {(f,sid):(p,T) for f,sid,p,T in CE.cohort('diag',50)}
    for fam in CE.FAMS:
        rr = [r for r in rows if r['method']=='V4b' and r['split']=='diag'
              and r['family']==fam and r['rmse'] is not None]
        selected.extend([(fam,min(r['seed'] for r in rr)),
                         (fam,max(rr,key=lambda r:r['rmse'])['seed'])])
    out = []
    for fam,sid in dict.fromkeys(selected):
        p,T = scenes[fam,sid]
        pol = Audit()
        log = P.episode(p,pol,T=T,seed=sid)
        if log is None:
            out.append(dict(family=fam,seed=sid,diverged=True))
            continue
        scores = {key:[] for key in ('kg_over_true','model_k_over_true',
                      'residual32_error_N','residual128_error_N',
                      'residual32_prediction_N','residual128_prediction_N',
                      'reference_abs_m','ff_increment_abs_m')}
        bad_gain = n_gain = 0
        for a in pol.audit:
            j = a['j']
            if j < 0 or log[j,0] < 1. or a['mode'] or log[j,3] < .5:
                continue
            kt = log[j,8]
            if kt > 0:
                scores['model_k_over_true'].append(a['k']/kt)
                if np.isfinite(a['kg']) and a['kg'] > 0:
                    ratio = a['kg']/kt
                    scores['kg_over_true'].append(ratio)
                    bad_gain += int(ratio < 1/1.2 or ratio > 1/.75)
                    n_gain += 1
            scores['reference_abs_m'].append(abs(a['reference']))
            for steps,label in ((4,'32'),(16,'128')):
                h = steps*.008
                q = j+steps
                if q >= len(log):
                    continue
                # Freeze k, subtract the realized residual-coordinate encoder
                # displacement. Includes stiffness/nonlinearity and future FF
                # policy changes; this is not pure anatomical-motion error.
                dy = (log[q,2]-pol.responses[q])-(log[j,2]-pol.responses[j])
                realized = log[q,3]-log[j,3]-a['k']*dy
                predicted = a['k']*a['velocity']*min(h,pol.c1.tv)
                scores['residual'+label+'_error_N'].append(abs(predicted-realized))
                scores['residual'+label+'_prediction_N'].append(abs(predicted))
        scores['ff_increment_abs_m'] = np.abs(np.diff(pol.offsets)).tolist()
        losses=[]
        for event in pol.loss_entries:
            j=event['j']
            if j>=0:
                _,ht,_=P.surface(log[j,1],log[j,0],np.asarray(p,float))
                losses.append(dict(event,time_s=float(log[j,0]),physical_surface_velocity_posthoc=float(ht)))
        out.append(dict(family=fam,seed=sid,loss_entries=losses,metrics=E.common.task_metrics(log,20.8),
                        counts=pol.counts,first_motion_s=pol.first_motion_time,
                        accepted_gain_outside_historical_band_fraction=bad_gain/max(n_gain,1),
                        scores={key:stats(v) for key,v in scores.items()}))
        print('AUDIT',fam,sid,flush=True)
    # An execution/startup sensitivity check, fixed first diagnostic scene only.
    # This is not included in the requested four-method ranking or selection.
    interface = []
    from contact_update import ContactUpdate
    for fam,sid,p,T in CE.cohort('diag',1):
        for label,pol in [('direct-noFF-safePrior',V4b(motion=False)),
                           ('gain-noFF-originalPrior',E.build('C1+V1'))]:
            log=P.episode(p,pol,T=T,seed=sid)
            m=E.common.task_metrics(log,20.8)
            interface.append(dict(family=fam,seed=sid,configuration=label,
                                  rmse=float(m['rmse']) if np.isfinite(m['rmse']) else None))
    E.save('v4b_audit.json',dict(episodes=out,startup_interface_sensitivity=interface,
         protocol='First and worst V4b finite diagnostic scene per family, selected after evaluation. Truth is post hoc only. No retuning. Forecast errors include model and future-policy changes, not pure motion errors. First-scene sensitivity combines startup prior and interface timing; it cannot independently attribute their effects.'))


if __name__ == '__main__':
    main()
