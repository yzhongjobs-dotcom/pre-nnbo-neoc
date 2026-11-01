"""Independent summary reconstruction and serial representative replay/timing.

Run after row3_gls.py completes. Adds only audit and serial_timing summaries
to its JSON and regenerates MD. No hyperparameter selection or final access.
"""
import hashlib
import json
from pathlib import Path
import numpy as np
import row3_gls as R
from row3_gls_report import markdown


def scientific(obj):
    if isinstance(obj, dict):
        return {k:scientific(v) for k,v in obj.items() if k not in ('p50_ms','p99_ms','wall_s','ms')}
    if isinstance(obj, list):
        return [scientific(v) for v in obj]
    return obj


def main():
    p=R.OUTPUT
    r=json.loads(p.read_text(encoding='utf-8'))
    assert r['status']=='complete'
    assert len(r['tuning'])==12
    assert r['selected']==min(r['tuning'], key=lambda x:x['objective'])
    assert len(r['open_loop']['episodes'])==300
    assert len(r['diagnostic']['episodes'])==500
    assert len(r['sweep']['episodes'])==480
    for name in ('diagnostic','sweep'):
        assert R.summarize_closed(r[name]['episodes'])==r[name]['summary']
        assert len({(a['method'],a['group'],a['dose'],a['sid']) for a in r[name]['episodes']})==len(r[name]['episodes'])
    # Compare all Proposed sweep episodes with the pre-existing primary result.
    old=json.loads((R.S.ROOT/'experiments/motion_aware_scan/results_v2/sweep.json').read_text(encoding='utf-8'))
    reference={(a['kind'],a['dose'],a['sid']):a for a in old if a['cfg']=='Proposed'}
    for a in r['sweep']['episodes']:
        if a['method']!='Proposed': continue
        b=reference[a['group'],a['dose'],a['sid']]
        assert a['ok']==b['ok']
        for key in ('rmse','peak','loss','loss_max'):
            assert a[key]==b[key], (a['sid'], key,a[key],b[key])
    # All phases ran from the recorded scientific code/data snapshot.
    for rel,digest in r['metadata']['sources'].items():
        assert hashlib.sha256((R.S.ROOT/rel).read_bytes()).hexdigest()==digest,rel
    timing={}
    for scene in R.S.C.cohort('diag',1):
        fam,sid,params,T=scene
        # One common trajectory, all methods, run serially and warm-started.
        batch=R.open_job((scene,r['selected']['params'],r['baseline_E2']['params']))
        for row in batch:
            prev=next(a for a in r['open_loop']['episodes'] if a['sid']==sid and a['method']==row['method'])
            summary={k:R.stats(v) for k,v in row['segments'].items()}
            assert scientific(summary)==scientific(prev['segments']), (sid,row['method'])
        timing[fam]={}
        for row in batch:
            a=row['segments']['all']['ms']
            timing[fam][row['method']+'_estimator']=dict(p50_ms=float(np.median(a)),p99_ms=float(np.quantile(a,.99)))
        for method in ('Proposed','E4'):
            a=R.closed_job((method,fam,None,sid,params,T,r['selected']['params']))
            prev=next(a0 for a0 in r['diagnostic']['episodes'] if a0['sid']==sid and a0['method']==method)
            for key in ('rmse','peak','loss','loss_max','ok'):
                assert a[key]==prev[key],(fam,method,key)
            timing[fam][method+'_controller']={key:a[key] for key in ('p50_ms','p99_ms')}
        print('audited',fam,flush=True)
    r['serial_timing']=dict(protocol='one prespecified first diagnostic scene per family; this task uses no worker pool; external machine load uncontrolled; exclude first second',summary=timing)
    dummy=np.zeros((40,16)); dummy[:,0]=np.arange(40)*R.S.DT; dummy[:,5]=R.S.common.FD
    X,Y,Z,idx=R.S.design(dummy,40)
    from gls_estimator import GLS
    est=GLS(r['selected']['params'])
    V,eigenvalues,_=est.geometry(X,idx-39)
    r['prior_example']=dict(force_N=R.S.common.FD, retained_intervals=len(idx),
        complement_dimension=len(eigenvalues),unit_kernel_eigenvalues=eigenvalues.tolist(),
        projected_residual_rms_N=float(np.sqrt((est.alpha*R.S.common.FD)**2*eigenvalues.mean())))
    r['metadata']['audit_source_sha256']={str(p.relative_to(R.S.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
        for p in [Path(__file__),R.HERE/'row3_gls_report.py',
                  R.S.ROOT/'experiments/motion_aware_scan/cohort_eval.py',
                  R.S.ROOT/'experiments/motion_aware_scan/scenes_v2.py',
                  R.S.ROOT/'experiments/motion_aware_scan/sweep.py',
                  R.S.ROOT/'experiments/motion_aware_scan/variants.py',
                  R.S.ROOT/'experiments/motion_aware_scan/proposed.py',
                  R.S.ROOT/'experiments/motion_aware_scan/imm.py',
                  R.S.ROOT/'experiments/cobot_scan_benchmark/scenarios.py']}
    r['interpretation']=dict(replace_projection=False,
        reason='Tail improvement in stiff push estimation does not by itself justify replacement; posterior gate loses availability and is conditionally miscalibrated. See paired closed-loop tables.',
        final_confirmation_used=False)
    r['audit']=dict(counts_and_unique_pairing=True,closed_summary_and_bootstrap_exact=True,
        all_240_Proposed_sweep_results_exact=True,recorded_source_hashes=True,
        representative_15_open_and_10_closed_scientific_replays_exact=True)
    R.save(r)
    p.with_suffix('.md').write_text(markdown(r),encoding='utf-8')
    print('AUDIT PASS',flush=True)


if __name__=='__main__':
    main()
