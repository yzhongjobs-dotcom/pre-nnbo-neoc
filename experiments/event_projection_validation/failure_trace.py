"""Trace decisive paired episodes, scored with truth only after the run."""
import json, numpy as np
from validate_closed_loop import cohort, common, FINAL, EventV1, EventLean, ServoID, JointPolicy, MeasuredOnly
from pathlib import Path

class TraceMixin:
    def reset(self,platform):
        super().reset(platform)
        self.monitor=[]
    def _trace(self,t,mode):
        super()._trace(t,mode)
        self.monitor.append([float(t),float(self.kbar),float(self.k_hat),
             float(self.k_good),float(self.kg_rel),str(self.fit.get('event_type','base')) if self.fit else 'none',
             int(mode)])
class TraceBaseline(TraceMixin,ServoID):pass
class TraceCandidate(TraceMixin,EventV1):pass

def eval_case(family,sid):
    scene=next(s for s in cohort('diag',5) if s[0]==family and s[1]==sid)
    _,_,p,T=scene;ret={}
    for name,cls in [('base',TraceBaseline),('candidate',TraceCandidate)]:
        pol=cls(single_model=False,**FINAL)
        log,_=common.run(MeasuredOnly(pol),p,T,sid)
        m=common.task_metrics(log,20.8)
        data=[]
        for rec in pol.monitor:
            i=int(np.searchsorted(log[:,0],rec[0]))
            i=min(max(0,i),len(log)-1)
            data.append(dict(t=rec[0],kbar=rec[1],khat=rec[2],kgood=rec[3],
                rel=rec[4],event=rec[5],mode=rec[6],
                truek=float(log[i,8]),force=float(log[i,3]),
                y=float(log[i,2]),cmd=float(log[i,11])))
        forces=np.asarray([x['force'] for x in data])
        i=int(np.nanargmax(forces))
        around=data[max(0,i-15):min(len(data),i+16)]
        first_gt20=next((x['t'] for x in data if x['force']>20),None)
        ret[name]=dict(rmse=float(m['rmse']),peak=float(m['peak']),
                       first_gt20=first_gt20,first_peak_time=data[i]['t'],
                       max_peak_context=around,
                       events=[x for x in data if x['event'] not in ('base','none')][:60],
                       trace=data)
    return ret

def main():
    cases=[('T2f-stiff',114001),('T2m-stiff',113004)]
    result={}
    for family,sid in cases:
        r=eval_case(family,sid)
        result[family+'|'+str(sid)]=r
        for m,v in r.items():
            print(family,sid,m,'RMSE',v['rmse'],'peak',v['peak'],
                  'first>20',v['first_gt20'],'peak_t',v['first_peak_time'],flush=True)
            for x in v['max_peak_context'][::5]:
                print('   ',x,flush=True)
    path=Path(__file__).resolve().parent/'failure_trace.json'
    path.write_text(json.dumps(result,ensure_ascii=False,indent=1),encoding='utf8')

if __name__=='__main__':main()
