import sys,json,requests,time
from pathlib import Path
sys.dont_write_bytecode=True
sys.path.insert(0,'D:/桌面/科研Agent Harness设计/ARIS/tools')
from openalex_fetch import OpenAlexClient
ROOT=Path(__file__).parent;c=OpenAlexClient();rows=[]
filters=['title.search:robotic ultrasound','title.search:stiffness estimation','title.search:motion compensation','cites:W2120315429']
try:
 d=requests.get('https://api.openalex.org/works/https://doi.org/10.1109/TIM.2022.3200360',timeout=20).json()
 aid=d['authorships'][0]['author']['id'].split('/')[-1];print('JIANG',aid,flush=True);filters.append('authorships.author.id:'+aid)
except Exception as e:print('author error',str(e),flush=True)
for f in filters:
 try:
  r=requests.get('https://api.openalex.org/works',params={'filter':f+',publication_year:2022-2026','per_page':200},timeout=25);r.raise_for_status();d=r.json()
  row={'filter':f,'url':r.url,'count':d.get('meta',{}).get('count'),'results':[c._parse_work(w) for w in d.get('results',[])]}
 except Exception as e:row={'filter':f,'error':str(e),'results':[]}
 rows.append(row);(ROOT/'recent_oa_filters.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2),encoding='utf-8')
 print(f,'TOTAL',row.get('count'),'N',len(row['results']),row.get('error',''),flush=True)
 if f.startswith(('cites','authorships')):
  for w in row['results']:print(w['publication_year'],w['title'],w['doi'],flush=True)
 time.sleep(1.5)
