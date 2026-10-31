import sys,json,time
from pathlib import Path
sys.dont_write_bytecode=True
sys.path.insert(0,'D:/桌面/科研Agent Harness设计/ARIS/tools')
import semantic_scholar_fetch as s2
from recent_search import QUERIES,ROOT
rows=[]
for i,q in enumerate(['robotic ultrasound motion']+QUERIES):
 try:
  r=s2.search(q,max_results=20,year='2022-2026',fields_of_study='Computer Science,Engineering,Medicine',publication_types='JournalArticle,Conference')
  row={'source':'S2','query':q,**r}
 except Exception as e: row={'source':'S2','query':q,'error':str(e),'data':[]}
 rows.append(row)
 (ROOT/'recent_s2.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2),encoding='utf-8')
 print(f'S{i:02d} {q} n={len(row["data"])} {row.get("error","")}',flush=True)
 for w in row['data']: print(w.get('year'),w.get('title'),w.get('externalIds',{}).get('DOI'),flush=True)
 time.sleep(1.2)
