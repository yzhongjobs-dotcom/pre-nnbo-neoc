"""Fetch only the pinned, screened source files; large PDFs remain outside the repository."""
import argparse,concurrent.futures,hashlib,json,pathlib,re,subprocess,urllib.request,urllib.parse
ROOT=pathlib.Path(__file__).resolve().parent

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cache',type=pathlib.Path,required=True);ap.add_argument('--download',action='store_true');ap.add_argument('--ids',nargs='*');ap.add_argument('--idea-screen',action='store_true');args=ap.parse_args()
    m=json.loads((ROOT/'阻抗文献重建_来源与筛选.json').read_text())
    records=m['files']+m['idea_stage']['files']+[x for x in m['extra_primary_files'] if x.get('ok')]
    if args.ids:records=[x for x in records if x.get('id') in args.ids]
    def get(r):
        pdf=r['path'].endswith('.pdf');dest=args.cache/('pdf' if pdf else 'idea')/pathlib.Path(r['path']).name;dest.parent.mkdir(parents=True,exist_ok=True)
        url=r.get('raw_url') or r.get('url') or 'https://raw.githubusercontent.com/'+m['repository']+'/'+m['commit']+'/'+urllib.parse.quote(r['path'])
        if not dest.exists():
            if not args.download:raise FileNotFoundError(dest)
            dest.write_bytes(urllib.request.urlopen(url,timeout=120).read())
        b=dest.read_bytes();sha256=hashlib.sha256(b).hexdigest()
        if sha256!=r['sha256']:raise ValueError('SHA256 mismatch: '+r['path'])
        if r.get('sha') and hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()!=r['sha']:raise ValueError('Git blob mismatch: '+r['path'])
        if pdf:
            out=args.cache/'text'/(dest.stem+'.txt');out.parent.mkdir(exist_ok=True);subprocess.run(['pdftotext',str(dest),str(out)],check=True,stderr=subprocess.PIPE)
        return r.get('id',dest.name)
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:done=list(pool.map(get,records))
    print(json.dumps({'verified_files':len(done),'ids':done},ensure_ascii=False))
    if args.idea_screen:
        rows=[json.loads(s) for s in (args.cache/'idea/LITERATURE_CORPUS.jsonl').read_text().splitlines() if s.strip()];titles={}
        for r in rows:
            title=r.get('title','')
            if re.search(r'(impedance|admittance|hunt.crossley|real.time iteration|neighboring.extremal)',title,re.I):
                key=re.sub('[^a-z0-9]','',title.lower())
                if key not in titles or len(str(r))>len(str(titles[key])):titles[key]=r
        (args.cache/'idea-title-screen.json').write_text(json.dumps(list(titles.values()),ensure_ascii=False,indent=2))
        print(json.dumps({'discovery_records':len(rows),'broad_title_candidates':len(titles),'note':'Candidates only, not admitted priors; contains aliases and unrelated electrical impedance.'},ensure_ascii=False))
if __name__=='__main__':main()
