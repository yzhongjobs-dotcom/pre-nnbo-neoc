import sys,json,time,concurrent.futures
from pathlib import Path
sys.dont_write_bytecode=True
sys.path.insert(0,'D:/桌面/科研Agent Harness设计/ARIS/tools')
from openalex_fetch import OpenAlexClient
ROOT=Path(__file__).parent
QUERIES=[
'robotic ultrasound motion aware',
'robotic ultrasound respiratory motion compensation',
'robotic ultrasound patient motion force control',
'robotic ultrasound force control Jiang',
'robotic ultrasound force control reinforcement learning',
'robotic ultrasound adaptive stiffness estimation',
'robotic ultrasound disturbance observer',
'robotic ultrasound model predictive force control',
'contact stiffness estimation non stationary environment',
'online impedance estimation moving environment force control',
'viscoelastic estimation motion compensation robotic surgery',
'force based physiological motion estimation prediction',
'contact force motion estimation feedforward compensation',
'robotic polishing stiffness estimation moving surface',
'robotic palpation online tissue stiffness motion',
'motion aware robotic ultrasound Jiang Navab',
]
def run(q):
 t=time.time()
 try:
  r=OpenAlexClient().search_works(q,max_results=20,publication_year='2022-2026')
  return {'source':'OpenAlex','query':q,'year':'2022-2026','max':20,'elapsed':round(time.time()-t,1),'results':r}
 except Exception as e:return {'source':'OpenAlex','query':q,'error':str(e),'results':[]}
if __name__=='__main__':
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
  rows=list(ex.map(run,QUERIES))
 (ROOT/'recent_openalex.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2),encoding='utf-8')
 for i,r in enumerate(rows,1):
  print(f'Q{i:02d} {r["query"]} n={len(r["results"])} {r.get("error","")}')
  for w in r['results']: print(w.get('publication_year'),w.get('title'),w.get('doi'))
