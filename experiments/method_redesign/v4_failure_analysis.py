"""Post-selection diagnostics. No tuning; truth is used ONLY after each episode.

Replays first and worst finite diagnostic scene per family; retains only
aggregate errors, gain ratios, saturation/recovery fractions, not tick files.
"""
import json, math
import numpy as np
import v4_evaluate as E
from v4_architecture import V4, P, CE


class AuditV4(V4):
    def reset(self, platform):
        super().reset(platform)
        self.audit = []

    def _observe(self, o):
        self.audit_n = o.n
        super()._observe(o)

    def _motion_state(self, st):
        out = super()._motion_state(st)
        ob = self.observer
        # Contains causal internal estimates only. No truth seen by controller.
        self.audit.append(dict(j=self.audit_n-1, k=self.kuse, g0=float(st[7]),
                    y=float(st[4]), kg=float(ob.k_good), ready=bool(self.motion_ready),
                    v=float(ob.tracker.x[1]), a=float(ob.tracker.x[2]),
                    df32=float(out[16+31]), df128=float(out[16+127])))
        return out


def stats(x):
    a=np.asarray(x,float); a=a[np.isfinite(a)]
    return dict(n=len(a),median=float(np.median(a)) if len(a) else None,
                p95=float(np.percentile(a,95)) if len(a) else None,
                mean=float(a.mean()) if len(a) else None)


def main():
    params=json.loads((E.OUT/'v4_tuning.json').read_text())['best']['params']
    rows=json.loads((E.OUT/'v4_episode_summary.json').read_text())
    scenes={(f,sid):(f,sid,p,T) for f,sid,p,T in CE.cohort('diag',50)}
    selected=[]
    for fam in CE.FAMS:
        rr=[r for r in rows if r['split']=='diag' and r['family']==fam and r['method']=='V4' and r['rmse'] is not None]
        if rr:
            selected.extend([(fam,min(r['seed'] for r in rr)),(fam,max(rr,key=lambda r:r['rmse'])['seed'])])
    out=[]
    for fam,sid in dict.fromkeys(selected):
        _,_,p,T=scenes[fam,sid]
        pol=AuditV4(**params)
        log=P.episode(p,pol,T=T,seed=sid)
        if log is None:
            out.append(dict(family=fam,seed=sid,diverged=True)); continue
        scores={k:[] for k in ('kg_over_true','model_k_over_true','g0_error_abs_N',
                               'df32_error_abs_N','df128_error_abs_N','df32_abs_N','df128_abs_N')}
        ready=0; n=0
        for a in pol.audit:
            j=a['j']
            if log[j,0] < 1. or log[j,3] < .5:
                continue
            n+=1; ready+=a['ready']
            kt=log[j,8]
            if kt<=0:
                continue
            scores['model_k_over_true'].append(a['k']/kt)
            scores['g0_error_abs_N'].append(abs(a['g0']-log[j,3]))
            if not a['ready']:
                continue
            scores['kg_over_true'].append(a['kg']/kt)
            for steps,h in ((4,32),(16,128)):
                if j+steps>=len(log) or log[j+steps,3]<.5:
                    continue
                # Realised-trajectory drift diagnostic using frozen TRUE local
                # tangent stiffness. Includes contact nonlinearity and changing
                # stiffness; not a pure anatomical-motion error.
                actual=log[j+steps,3]-log[j,3]-kt*(log[j+steps,2]-log[j,2])
                scores[f'df{h}_error_abs_N'].append(abs(a[f'df{h}']-actual))
                scores[f'df{h}_abs_N'].append(abs(a[f'df{h}']))
        result=dict(family=fam,seed=sid,metrics=E.common.task_metrics(log,20.8),
                    motion_ready_fraction=ready/max(n,1), counters=pol.counts,
                    first_recovery_s=next((float(d[0]) for d in pol.details if d[5]),None),
                    first_motion_s=next((float(log[a['j'],0]) for a in pol.audit if a['ready']),None),
                    gain_outside_screened_band_fraction=(float(np.mean((np.asarray(scores['kg_over_true'])<1/1.2) | (np.asarray(scores['kg_over_true'])>1/.75))) if scores['kg_over_true'] else None),
                    errors={k:stats(v) for k,v in scores.items()})
        out.append(result)
        print(fam,sid,'ready',result['motion_ready_fraction'],'kg',result['errors']['kg_over_true'],flush=True)
    E.save('v4_failure_analysis.json',dict(selection='first and worst finite RMSE diagnostic scene per family after frozen evaluation',
        scope='Post hoc descriptive diagnostics, not additional selection/tuning; no causal module-isolation claim.',
        truth_boundary='True force/stiffness accessed only in scoring after the policy episode completes.',
        forecast_metric='delta_f(h) compared to Ftrue(t+h)-Ftrue(t)-ktrue(t)*(y(t+h)-y(t)); frozen-local-tangent approximation, includes nonlinearity/stiffness changes.',
        rows=out))


if __name__=='__main__':
    main()
