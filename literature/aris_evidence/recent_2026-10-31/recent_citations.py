import sys,json,requests,time
from pathlib import Path
sys.dont_write_bytecode=True
sys.path.insert(0,'D:/桌面/科研Agent Harness设计/ARIS/tools')
from openalex_fetch import OpenAlexClient
ROOT=Path(__file__).parent
rows=[]
c=OpenAlexClient()
seeds=[('Rossi','W2416622902'),('Moreira','W2019057824'),('DOSMAC','W7123713864'),('Beber','W4401415597')]
try:
 r=requests.get('https://api.crossref.org/works',params={'query.title':'Force tracking with feed-forward motion estimation for beating heart surgery','rows':3},timeout=20).json()['message']['items']
 (ROOT/'yuen_crossref.json').write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8')
 for w in r:print('YUEN LOOKUP',w.get('title'),w.get('DOI'),flush=True)
 doi=r[0]['DOI'];w=c.get_work(doi);seeds.append(('Yuen',w['openalex_id']))
except Exception as e:print('Yuen identity error',str(e),flush=True)
for name,wid in seeds:
 params={'filter':f'cites:{wid},publication_year:2022-2026','per_page':200}
 try:
  rr=requests.get('https://api.openalex.org/works',params=params,timeout=25);rr.raise_for_status();d=rr.json()
  row={'seed':name,'id':wid,'query':rr.url,'meta':d.get('meta'),'results':[c._parse_work(w) for w in d.get('results',[])]}
 except Exception as e:row={'seed':name,'id':wid,'error':str(e),'results':[]}
 rows.append(row);(ROOT/'recent_citations.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2),encoding='utf-8')
 print('CITES',name,'N',len(row['results']),row.get('meta'),row.get('error',''),flush=True)
 for w in row['results']:print(w['publication_year'],w['title'],w['doi'],flush=True)
 time.sleep(1)
