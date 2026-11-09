"""Three additive scenes, frozen baseline r2 methods; no development/tuning.

Reuses baseline_compare's plant, execution, calibration and causal replay.
Run --stage check/pilot/test/capture/replay/report, in that order.
"""
import argparse
import hashlib
import json
import math
import time
from pathlib import Path

import baseline_compare as B  # sets BLAS=1 before importing numpy
import baseline_report as F
import numpy as np

SCENES = ('breath_rib', 'breath_fast', 'breath_rib_fast')
LABELS = dict(zip(SCENES, ('呼吸＋跨肋骨', '呼吸＋突然推拉', '呼吸＋跨肋骨＋突然推拉')))
BASE, SENSOR_BASE = 674500, 957000
ORIGINAL_LOAD = B.load


def cells():
    return [(c, 'relaxed' if c == 'breath_fast' else 'rib', dy, ph, i)
            for c in SCENES for dy in ('neural_white', 'viscous') for ph in (14., 28.)
            for i in range(3 if c == 'breath_fast' else 9)]


def sensor_seed(cell, dev=False):
    assert not dev
    c, _, dy, ph, i = cell
    return SENSOR_BASE + 1000*SCENES.index(c) + 100*(dy == 'viscous') + 10*(ph == 28.) + i


def load(name):
    return ORIGINAL_LOAD('overlay_calibration.json' if name == 'baseline_calibration.json' else name)


def setup():
    B.TEST_BASE = BASE
    B.sensor_seed = sensor_seed
    B.load = load


def chosen():
    configs = dict(ORIGINAL_LOAD('baseline_test_r2.json')['chosen'])
    # The adopted default, not the tuned configuration that was rejected.
    configs['Huber-D4'] = configs.pop('Huber-current')
    return {m: configs[m] for m in B.METHODS}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def equal(a,b):
    if isinstance(a,np.ndarray):
        return np.array_equal(a,b)
    if isinstance(a,(list,tuple)):
        return len(a)==len(b) and all(equal(x,y) for x,y in zip(a,b))
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(equal(a[k],b[k]) for k in a)
    return a==b


def check():
    """Check all actual test drives, not merely scenario labels."""
    rows = []
    for cell in cells():
        sc = B.scene(cell, BASE, 0.)
        c = cell[0]
        assert sc['cond'] == c and sc['cell'] == cell[:2]
        assert sc['FD'] == 8. and sc['T'] == 10.
        fast = c != 'breath_rib'
        assert [e['kind'] for e in sc['event_info']] == (['fast'] if fast else [])
        assert not np.any(sc['act']) and not np.any(sc['fext'])
        total = sc['motion_breath'] + sc['motion_cardiac']
        for wave in sc['motion_components']:
            total = total + wave
        assert np.array_equal(total, sc['u'])
        assert len(sc['motion_components']) == int(fast)
        assert np.ptp(sc['motion_breath']) > .002
        crossed = [float(x) for x, *_ in sc['ribs'] if sc['xcmd'].min() <= x <= sc['xcmd'].max()]
        assert bool(crossed) == (c != 'breath_fast')
        if fast:
            e = sc['event_info'][0]
            assert .002 <= e['amplitude'] <= .012 and e['end'] < sc['T']
        if c == 'breath_fast':
            old = B.scene(('shift_fast', *cell[1:]), BASE, 0.)
            for k in old:
                if k in ('cond', 'cell'):
                    continue
                assert equal(old[k],sc[k]), k
        rows.append(dict(cell=cell,source=sc['source_cond'],ribs_crossed=len(crossed),
                         events=sc['event_info'],breath_pp=sc['breath_pp'],breath_period=sc['breath_period'],
                         sensor_seed=sensor_seed(cell)))
    assert len(set(r['sensor_seed'] for r in rows)) == len(rows)
    # Guard all old scenes bitwise against the pre-change generator, loaded
    # from Git (without touching the frozen reference or live source).
    import subprocess
    from types import ModuleType
    original = subprocess.check_output(['git', 'show', '858f3c1:experiments/rossi_failure_study/scenes3.py'], text=True, encoding='utf-8')
    module = ModuleType('original_scenes3'); module.__file__ = str(B.HERE/'scenes3.py')
    exec(compile(original, module.__file__, 'exec'), module.__dict__)
    old_checked = 0
    for c in B.S3.CONDITIONS:
        for tissue in B.S3.TISSUES:
            old = module.cell_scene3(c,tissue,base=BASE)
            new = B.S3.cell_scene3(c,tissue,base=BASE)
            assert old.keys() == new.keys()
            for k, v in old.items():
                assert equal(v,new[k]), (c,tissue,k)
            old_checked += 1
    obs = B.P3.Obs(np.ones((2,26)),2)
    assert np.isnan(np.delete(obs.log, B.OBS_COLS, axis=1)).all()
    protected = ['physics_plant3.py','method3.py','id_baselines.py','ctrl_baselines.py',
                 'd_robust.py','baseline_round2.py']
    plan = dict(passed=True,base=BASE,sensor_base=SENSOR_BASE,cells=cells(),rows=rows,
                chosen=chosen(),old_scenes_bitwise_equal=old_checked,truth_mask=True,
                protected_hashes={p:sha(B.HERE/p) for p in protected},
                config_source_sha256=sha(B.OUT/'baseline_test_r2.json'),
                note='No tuning; 9 seeds for rib-bearing scenes, 3 for fast shift; original 2 shapes x 2 phases.')
    B.dump('overlay_plan.json',plan)
    print('Checked 84 drives, 36 historical scenes, observation masking and frozen configurations.',flush=True)


def execute(job):
    setup()
    return B.execute(job)


def calibration_job(job):
    setup()
    return B.calibration_job(job)


def capture_job(row):
    setup()
    return B.capture_job(row)


def replay_job(job):
    setup()
    return B.replay_job(job)


def pilot():
    setup()
    cell = cells()[-1]
    cal = calibration_job((cell,False))
    B.dump('overlay_calibration.json',dict(positions={cal['key']:cal['position']},rows=[cal]))
    rows = B.run([(m,cfg,cell,BASE,0.,False,False) for m,cfg in chosen().items()],worker=execute)
    B.dump('overlay_pilot.json',dict(rows=rows))
    for r in rows:
        print(r['method'],round(r['seconds'],2),'s',round(math.log1p(r['J']),3),r['failed'],flush=True)


def test():
    plan = ORIGINAL_LOAD('overlay_plan.json')
    assert plan['chosen'] == chosen()
    setup()
    start = time.perf_counter()
    rr = B.run([(cell,False) for cell in cells()],worker=calibration_job)
    B.dump('overlay_calibration.json',dict(positions={r['key']:r['position'] for r in rr},rows=rr,
                                          threshold_N=.5,approach_speed=.005,extra_simulated_seconds=2.))
    path = B.OUT/'overlay_test.json'
    saved = ORIGINAL_LOAD(path.name) if path.exists() else dict(chosen=chosen(),base=BASE,sensor_base=SENSOR_BASE,rows=[])
    assert saved['chosen'] == chosen()
    for m,cfg in chosen().items():
        complete = {B.key(r) for r in saved['rows'] if r['method'] == m}
        jobs = [(m,cfg,c,BASE,0.,True,False) for c in cells() if c not in complete]
        if jobs:
            saved['rows'] += B.run(jobs,worker=execute)
            saved['seconds'] = saved.get('seconds',0.) + time.perf_counter()-start
            start = time.perf_counter()
            B.dump(path.name,saved)
        print(m,B.L([r for r in saved['rows'] if r['method'] == m]),flush=True)


def capture():
    sources = [r for r in ORIGINAL_LOAD('overlay_test.json')['rows'] if r['method'] in ('Huber-D4','D4')]
    rows = B.run(sources,worker=capture_job)
    B.dump('overlay_capture.json',dict(rows=rows))


def replay():
    test_data = ORIGINAL_LOAD('overlay_test.json')
    sources = [dict(r,retain_series=True) for r in test_data['rows'] if r['method'] in ('Huber-D4','D4')]
    path = B.OUT/'overlay_module.json'
    saved = ORIGINAL_LOAD(path.name) if path.exists() else dict(rows=[])
    for m in B.ID_METHODS:
        done = {(r['source'],B.key(r)) for r in saved['rows'] if r['method'] == m}
        jobs = [(s,m,test_data['chosen'][m],0.) for s in sources if (s['method'],B.key(s)) not in done]
        if not jobs:
            continue
        rr = B.run(jobs,worker=replay_job)
        folder = B.OUT/'baseline_replay';folder.mkdir(exist_ok=True)
        for r in rr:
            arrays = r.pop('series')
            for kind,g in r['groups'].items():
                for field in ('error_samples','accepted_errors','times'):
                    arrays[kind+'_'+field] = np.asarray(g.pop(field))
            name = f"overlay_{m}_{r['source']}_{r['cond']}_{r['dynamics']}_{int(r['phase'])}_{r['i']}.npz"
            np.savez_compressed(folder/name,**arrays)
            r['samples_file'] = 'baseline_replay/'+name
        saved['rows'] += rr
        B.dump(path.name,saved)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--stage',choices=('check','pilot','test','capture','replay'),required=True)
    args = p.parse_args()
    setup()
    globals()[args.stage]()
