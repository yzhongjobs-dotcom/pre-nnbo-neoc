"""Pinned public phantom exports: observational task audit, not controller replay.

Downloads go to an external cache. No raw B-mode images are redistributed.
Original high-rate verdicts and independently calculated ~1 Hz proxies stay
separate. The exported five-pose revision is not the later seven-pose paper.
"""
from pathlib import Path
from urllib.request import urlopen
import argparse
import csv
import hashlib
import io
import json
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
COMMIT = 'c765e23be3b2ccebdd971b521d9e1f8c2662d9f4'
BASE = f'https://raw.githubusercontent.com/suung23/FR5-for-RUS/{COMMIT}/'
ROOT = 'Paper/icra_bladder/'
EPISODES = ['ep0001_pose01_rule', 'ep0002_pose01_search', 'ep0003_pose01_hold',
    'ep0004_pose02_search', 'ep0005_pose02_hold', 'ep0006_pose02_rule',
    'ep0007_pose03_rule', 'ep0008_pose03_hold', 'ep0009_pose03_search',
    'ep0010_pose04_rule', 'ep0011_pose04_search', 'ep0012_pose04_hold',
    'ep0013_pose05_rule', 'ep0014_pose05_search', 'ep0015_pose05_hold']


def clean(obj):
    if isinstance(obj, dict): return {k: clean(v) for k,v in obj.items()}
    if isinstance(obj, list): return [clean(v) for v in obj]
    if isinstance(obj, float) and not np.isfinite(obj): return None
    return obj


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cache', default='/tmp/r15-real-scan/cache')
    parser.add_argument('--inspect', action='store_true', help='inspect fixed pose-1 frames in cache; no raw images stored in project')
    args = parser.parse_args()
    cache = Path(args.cache)
    manifest = []
    def get(path):
        dst = cache/path
        if not dst.exists():
            dst.parent.mkdir(parents=True, exist_ok=True)
            dst.write_bytes(urlopen(BASE+path, timeout=45).read())
        content = dst.read_bytes()
        return content, dict(path=path, bytes=len(content),
                             sha256=hashlib.sha256(content).hexdigest())
    paths = [ROOT+'data/verdicts.json', ROOT+'README.md', ROOT+'export_results.py',
             'policy_learning/results/experiment_freeze.json',
             'policy_learning/rus_policy/view_quality.py']
    paths += [ROOT+f'data/frames/{ep}/index.csv' for ep in EPISODES]
    downloaded = list(ThreadPoolExecutor(6).map(get, paths))
    sources = {}
    for path, (content, meta) in zip(paths, downloaded):
        sources[path] = content
        manifest.append(meta)
    verdicts = json.loads(sources[ROOT+'data/verdicts.json'])
    records, traces = [], {}
    keys = ['found','find_time_s','hold_window_s','hold_good_fraction',
        'hold_worst_loss_s','force_mean_n','force_max_n',
        'time_above_limit_s','hold_force_max_n','area_max']
    for ep in EPISODES:
        rows = list(csv.DictReader(io.StringIO(sources[ROOT+f'data/frames/{ep}/index.csv'].decode())))
        val = lambda key: np.array([float(r[key]) if r[key] else np.nan for r in rows])
        t, area, dx, fn = (val(k) for k in ['t','area','centroid_dx','force_n'])
        # Necessary image geometry conditions only: component connectivity
        # is absent from index.csv. Not a diagnostic or original success test.
        proxy = (area >= .04) & (np.abs(dx) <= .30)
        band = np.isfinite(fn) & (fn >= 2) & (fn <= 4)
        # Duration-weight each exported sample forward; final sample has zero
        # weight. No interpolation over missing forces or high-rate peaks.
        dt = np.r_[np.diff(t), 0.]
        duration = lambda mask: float(dt[mask].sum())
        v = verdicts[ep]
        record = dict(episode=ep, pose=int(v['pose']), condition=v['condition'],
            exported_frames=len(rows), exported_duration_s=float(t[-1]-t[0]),
            force_available_frames=int(np.isfinite(fn).sum()),
            force_available_duration_s=duration(np.isfinite(fn)),
            force_band_2_4N_duration_s=duration(band),
            force_band_without_area_centroid_proxy_s=duration(band & ~proxy),
            force_band_with_zero_segmented_area_s=duration(band & (area==0)),
            proxy_duration_s=duration(proxy),
            exported_peak_force_n=float(np.nanmax(fn)) if np.isfinite(fn).any() else None,
            author_verdict={k:v.get(k) for k in keys})
        records.append(record)
        traces[ep] = (t, area, dx, fn)
    groups = {}
    for condition in ['rule','search','hold']:
        rr = [r for r in records if r['condition']==condition]
        vv = [r['author_verdict'] for r in rr]
        total_band = sum(r['force_band_2_4N_duration_s'] for r in rr)
        found = [v for v in vv if v['found']]
        groups[condition] = dict(episodes=len(rr),
            exported_frames=sum(r['exported_frames'] for r in rr),
            force_available_frames=sum(r['force_available_frames'] for r in rr),
            author_found_count=sum(bool(v['found']) for v in vv),
            author_median_find_s_among_found=float(np.median([v['find_time_s'] for v in found])) if found else None,
            author_median_hold_good_fraction_among_found=float(np.median([v['hold_good_fraction'] for v in found])) if found else None,
            author_max_force_n=max(v['force_max_n'] for v in vv),
            author_episodes_above_5N=sum(v['force_max_n']>5 for v in vv),
            export_force_band_duration_s=total_band,
            export_band_without_proxy_fraction=sum(r['force_band_without_area_centroid_proxy_s'] for r in rr)/total_band,
            export_band_zero_area_fraction=sum(r['force_band_with_zero_segmented_area_s'] for r in rr)/total_band)
    protocol = dict(repository='suung23/FR5-for-RUS', commit=COMMIT,
        branch_at_discovery='eval/rule-v1',episodes=EPISODES,
        sample_rate='approximately 1 Hz, exported; missing force rows are not interpolated',
        unit_of_replication='5 paired starting poses, not individual frames',
        force_band_N=[2,4], force_band_status='descriptive band, not a safety or clinical threshold',
        independent_proxy='segmented area >=0.04 and abs(centroid_dx)<=0.30; necessary conditions only',
        excluded='Other episodes lack matching exported frames; no pooling of later phantom position B',
        limitations=['No action/pose stream, raw force series, or segmentation truth labels',
            'Force matching uses 0.1 s rounded timestamps in the source exporter; missingness is large and unequal, including a whole episode',
            'Band-conditioned proportions describe available exported observations only, not unbiased episode prevalence',
            'Published high-rate verdicts are source summaries, not independently reproduced',
            'No causal controller comparison, clinical validation or dynamic simulation identification',
            'Static hold is a weak reference; rule/search are existing alternatives, not our methods'],
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (HERE/'protocol.json').write_text(json.dumps(protocol,indent=2,ensure_ascii=False)+'\n')
    (HERE/'source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (HERE/'results.json').write_text(json.dumps(clean(dict(groups=groups,episodes=records)),
        indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    # All five paired poses; no hand-picked successful example.
    fig, axes = plt.subplots(5,2,figsize=(10,11),sharex=True)
    colors = dict(rule='#0072B2',search='#D55E00',hold='#009E73')
    for pose in range(1,6):
        for ep,(t,area,dx,fn) in traces.items():
            if f'pose{pose:02d}' not in ep: continue
            condition=ep.split('_')[-1];color=colors[condition]
            axes[pose-1,0].plot(t,area,label=condition,color=color,lw=1)
            axes[pose-1,1].plot(t,fn,color=color,lw=1)
        axes[pose-1,0].axhline(.04,color='gray',ls='--',lw=.6)
        axes[pose-1,1].axhline(5,color='gray',ls='--',lw=.6)
        axes[pose-1,0].set_ylabel(f'Pose {pose}\nMask area / ROI')
        axes[pose-1,1].set_ylabel('Exported force (N)')
        for ax in axes[pose-1]: ax.grid(alpha=.2)
    axes[0,0].legend(loc='upper left',fontsize=8)
    for ax in axes[-1]: ax.set_xlabel('Time (s)')
    fig.suptitle('Public phantom exports: image geometry and contact force\n~1 Hz traces; source high-rate peaks require the original raw stream',fontsize=11)
    fig.tight_layout(rect=(0,0,1,.95));fig.savefig(HERE/'paired_traces.png',dpi=160);plt.close(fig)
    if args.inspect:
        fig,axes=plt.subplots(3,3,figsize=(9,9));inspection=[]
        for row,ep in enumerate(EPISODES[:3]):
            rows=list(csv.DictReader(io.StringIO(sources[ROOT+f'data/frames/{ep}/index.csv'].decode())))
            for col,target in enumerate([0,60,100]):
                rec=min(rows,key=lambda r:abs(float(r['t'])-target));n=int(rec['n'])
                images=[]
                for prefix in ['f','m']:
                    path=ROOT+f'data/frames/{ep}/{prefix}{n:04d}.png'
                    content,meta=get(path);images.append(plt.imread(io.BytesIO(content)))
                    inspection.append(dict(**meta,episode=ep,n=n,t=rec['t']))
                ax=axes[row,col];ax.imshow(images[0],cmap='gray')
                if images[1].max()>.5:ax.contour(images[1],levels=[.5],colors=['red'],linewidths=.7)
                ax.set_title(f'{ep.split("_")[-1]} t={rec["t"]}s; area={rec["area"]}',fontsize=9);ax.axis('off')
        fig.tight_layout();fig.savefig(cache/'inspection.png',dpi=150);plt.close(fig)
        (HERE/'inspection_manifest.json').write_text(json.dumps(inspection,indent=2)+'\n')
    print(json.dumps(groups,indent=2))


if __name__=='__main__': main()
