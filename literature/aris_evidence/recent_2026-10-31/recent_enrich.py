import sys,json,time,requests,concurrent.futures
from pathlib import Path
sys.dont_write_bytecode=True
sys.path.insert(0,'D:/桌面/科研Agent Harness设计/ARIS/tools')
from openalex_fetch import OpenAlexClient
import semantic_scholar_fetch as s
ROOT=Path(__file__).parent
dois=['10.1109/TIM.2022.3200360','10.1109/LRA.2022.3180440','10.1109/TMECH.2024.3443310','10.1109/ICRA57147.2024.10611406','10.1016/j.conengprac.2026.107163','10.1109/TIM.2025.3581663','10.1109/IROS60139.2025.11246299','10.1109/TII.2024.3477549','10.1109/TIE.2024.3451138','10.1109/LRA.2024.3349959','10.1109/LRA.2025.3604730','10.1109/IROS58592.2024.10802530','10.1109/TRO.2026.3734315']
def one(doi):
 row={'doi':doi}
 try:row['openalex']=OpenAlexClient().get_work(doi)
 except Exception as e:row['oa_error']=str(e)
 try:
  r=requests.get('https://api.crossref.org/works/'+doi,timeout=20);r.raise_for_status();w=r.json()['message'];row['crossref']={k:w.get(k) for k in ['DOI','title','container-title','author','published','published-online','published-print','resource','abstract','type']}
 except Exception as e:row['cr_error']=str(e)
 try:row['s2']=s.get_paper('DOI:'+doi)
 except Exception as e:row['s2_error']=str(e)
 return row
rows=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
 for row in ex.map(one,dois):
  rows.append(row);(ROOT/'recent_enrich.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2),encoding='utf-8')
  print(row['doi'],row.get('openalex',{}).get('title'), 'OA',row.get('oa_error','OK'),'CR',row.get('cr_error','OK'),'S2',row.get('s2_error','OK'),flush=True)
  if row['doi'].endswith(('107163','10802530','3734315')):
   print('ABSTRACT',row.get('openalex',{}).get('abstract') or row.get('s2',{}).get('abstract') or row.get('crossref',{}).get('abstract'),flush=True)
