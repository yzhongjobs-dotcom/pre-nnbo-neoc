"""Reproduce the pinned corpus and two limited mathematical boundary checks.

PDFs remain in an external cache. They are not copied into the project/Git.
These checks do not purport to refute a paper's full robot model or theorem.
"""
from pathlib import Path
import argparse,hashlib,json,math,subprocess,urllib.request

HERE=Path(__file__).resolve().parent


def math_checks():
    # Inner approximation of a disk by a diamond: outside the approximation is
    # not necessarily outside the true feasible region.
    point=(.6,.6)
    disk_norm=math.hypot(*point);diamond_norm=sum(abs(x) for x in point)
    assert disk_norm<1 and diamond_norm>1
    # Sum of squared distances to vertices is minimized at their centroid,
    # which need not maximize distance to the boundary (Chebyshev radius).
    vertices=[(0.,0.),(4.,0.),(0.,1.)]
    centroid=(4/3,1/3);radius=(5-math.sqrt(17))/2;incenter=(radius,radius)
    def margin(p):return min(p[0],p[1],(4-p[0]-4*p[1])/math.sqrt(17))
    def objective(p):return sum(sum((a-b)**2 for a,b in zip(p,v)) for v in vertices)
    assert objective(centroid)<objective(incenter)
    assert margin(centroid)<margin(incenter)
    return dict(
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        source_paper='R21, PDF page 10, equations (34)-(37), Proposition 2',
        source_page_visually_checked=True,
        inner_approximation=dict(point=point,true_disk_norm=disk_norm,inner_diamond_norm=diamond_norm),
        centroid_objective=dict(vertices=vertices,centroid=centroid,incenter=incenter,
            sum_squared_distances_centroid=objective(centroid),sum_squared_distances_incenter=objective(incenter),
            boundary_margin_centroid=margin(centroid),boundary_margin_incenter=margin(incenter)),
        conclusion='General construction alone does not imply outside-inner-set means infeasible, or centroid means maximum boundary margin. Additional conditions/definitions are needed.',
        scope='Independent algebraic counterexamples to generic inference only; no audit/refutation of the full underwater robot dynamics, experiments or all proofs; no new impedance-controller performance evidence.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache',type=Path,default=Path('/tmp/r15-paper-reverse'))
    parser.add_argument('--download',action='store_true')
    parser.add_argument('--ids',nargs='*')
    parser.add_argument('--write-math-checks',action='store_true')
    args=parser.parse_args()
    manifest=json.loads((HERE/'顶刊样本来源与阅读范围.json').read_text())
    selected=[d for d in manifest['files'] if not args.ids or d['id'] in args.ids]
    if args.ids:assert set(args.ids)<={d['id'] for d in selected},'Unknown paper id'
    for item in selected:
        path=args.cache/'pdf'/item['filename']
        if args.download and not path.exists():
            path.parent.mkdir(parents=True,exist_ok=True)
            req=urllib.request.Request(item['raw_url'],headers={'User-Agent':'research-source-audit'})
            with urllib.request.urlopen(req,timeout=120) as response:path.write_bytes(response.read())
        if not path.exists():raise FileNotFoundError(f'{path}: use --download if needed')
        data=path.read_bytes()
        assert hashlib.sha256(data).hexdigest()==item['sha256'],item['id']
        assert hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==item['git_blob_sha'],item['id']
        plain=args.cache/'plain'/(path.stem+'.txt');plain.parent.mkdir(exist_ok=True)
        if not plain.exists():subprocess.run(['pdftotext',str(path),str(plain)],check=True)
    checks=math_checks()
    if args.write_math_checks:
        (HERE/'顶刊样本数学核对.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(verified_pdfs=len(selected),source_commit=manifest['commit'],math_checks=checks),ensure_ascii=False,indent=2))


if __name__=='__main__':main()
