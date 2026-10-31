import sys,json,time
from pathlib import Path
sys.dont_write_bytecode=True
sys.path.insert(0,'D:/桌面/科研Agent Harness设计/ARIS/tools')
import semantic_scholar_fetch as s
ROOT=Path(__file__).parent
qs=['"ultrasound" +"motion-aware"','"ultrasound" +Jiang +motion','"stiffness" +"non-stationary"','"stiffness" +"moving environment"','"force" +"motion estimation"','"ultrasound" +"force control" +2026','"ultrasound" +"model predictive"','"ultrasound" +"disturbance observer"','"contact" +"stiffness estimation"','"beating heart" +"force"']
rows=[]
for q in qs:
 try:r={'query':q,**s.search_bulk(q,max_results=100,year='2022-2026',sort='citationCount:desc')}
 except Exception as e:r={'query':q,'error':str(e),'data':[]}
 rows.append(r);(ROOT/'recent_targeted.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2),encoding='utf-8')
 print(q,'N',len(r.get('data',[])),r.get('error',''),flush=True)
 for w in r.get('data',[])[:30]:print(w.get('year'),w['title'],w.get('externalIds',{}).get('DOI'),flush=True)
 time.sleep(1.5)
