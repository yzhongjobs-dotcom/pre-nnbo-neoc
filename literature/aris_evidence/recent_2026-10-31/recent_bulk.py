import sys,json,time,requests
from pathlib import Path
sys.dont_write_bytecode=True
sys.path.insert(0,'D:/桌面/科研Agent Harness设计/ARIS/tools')
import semantic_scholar_fetch as s
ROOT=Path(__file__).parent
rows=[]
for q in ['robotic ultrasound','stiffness estimation moving environment','force motion compensation beating heart']:
 try:
  r=s.search_bulk(q,max_results=100,year='2022-2026',sort='citationCount:desc')
  rows.append({'query':q,**r})
  print('QUERY',q,'TOTAL',r.get('total'),'N',len(r.get('data',[])),flush=True)
  for w in r.get('data',[]):print(w['year'],w['title'],w.get('externalIds',{}).get('DOI'),flush=True)
 except Exception as e:rows.append({'query':q,'error':str(e),'data':[]});print(type(e).__name__,str(e),flush=True)
 (ROOT/'recent_s2_bulk.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2),encoding='utf-8')
 time.sleep(2)
r=requests.get('https://api.openalex.org/works/https://doi.org/10.1109/ICRA.2016.7487243',timeout=20)
print('OA single',r.status_code,r.text[:300] if r.status_code!=200 else r.json().get('id'),flush=True)
