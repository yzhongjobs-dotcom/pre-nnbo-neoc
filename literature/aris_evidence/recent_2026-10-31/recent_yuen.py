import requests,json,sys
from pathlib import Path
ROOT=Path(__file__).parent
url='https://api.semanticscholar.org/graph/v1/paper/DOI:10.1109/TRO.2010.2053734/citations'
r=requests.get(url,params={'fields':'title,year,venue,authors,abstract,externalIds,url','limit':1000},timeout=25)
d={'source':'S2','url':r.url,'status':r.status_code}
if r.status_code==200:
 d['response']=r.json()
 for x in r.json().get('data',[]):
  w=x.get('citingPaper',{})
  if 2022<=(w.get('year') or 0)<=2026:print(w['year'],w['title'],w.get('externalIds',{}).get('DOI'));print(w.get('abstract'))
else:d['error']=r.text[:500];print(r.status_code,r.text[:500])
(ROOT/'recent_yuen_s2.json').write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
