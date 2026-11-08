"""Four fixed episodes, instance-only channel interventions; JSON to stdout.

Run: python -X utf8 experiments/rossi_failure_study/tail_channel_diagnosis.py
No result files are written. k means diagnostic truth_k = |K(20.8 Hz)|.
"""
import ast
import inspect
import json
import textwrap
from pathlib import Path
from types import MethodType

import overlay_compare as O
import numpy as np
import method3 as M
import physics_plant3 as P3
import memory_closed_loop as MC
import proposed

B = O.B
HERE = Path(__file__).resolve().parent


def cases():
    yield dict(label='A', cell=('breath_rib', 'rib', 'neural_white', 28., 6),
               base=674500, variant='nominal', seed=957016, stored_peak=35.5,
               rounded=True)
    rows = json.loads((HERE/'results/memory_closed_loop_true_base.json').read_text())
    for label, name, field, value in (
        ('slow_thin', 'Huber-true', 'variant', 'slow_thin'),
        ('combo', 'Huber-true six', 'scene', 'combo'),
        ('cough', 'Huber-true six', 'scene', 'cough'),
    ):
        r = max((r for r in rows if r['name'] == name and r[field] == value),
                key=lambda r: r['peak'])
        grid = ([c for c in O.cells() if c[0] == 'breath_rib']
                if label == 'slow_thin' else B.cells())
        q, cell = next((q, c) for q, c in enumerate(grid)
                       if (c[0], c[2], c[3], c[4]) ==
                       (r['scene'], r['dynamics'], r['phase'], r['i']))
        yield dict(label=label, cell=cell, base=774500, variant=r['variant'],
                   seed=(880000 if label == 'slow_thin' else 890000)+q,
                   seed_index=q, stored_peak=r['peak'], rounded=False)


def instrument_step(pol, run, capture):
    # Compile a private copy, bound ONLY to this policy instance. The inherited
    # ServoID/IMM/sensing paths still call it through their original super chain
    # using a private forwarding instance hook below (no class is patched).
    source = textwrap.dedent(inspect.getsource(proposed.MotionAware.step))
    tree = ast.parse(source)
    fn = tree.body[0]
    final = fn.body[-1]
    assert isinstance(final, ast.Return)
    assert ast.unparse(final.value) == 'sent + du / self.kbar + (dx - dq) / k_un'
    if run == 'ff_off':
        fn.body.insert(-1, ast.parse('dq = 0.0').body[0])
    fn.body[-1] = ast.Return(ast.Call(
        func=ast.Attribute(ast.Name('self', ast.Load()), '_tail_capture', ast.Load()),
        args=[ast.Call(ast.Name('locals', ast.Load()), [], []), final.value], keywords=[]))
    ast.fix_missing_locations(tree)
    ns = {}
    exec(compile(tree, '<instance tail diagnosis>', 'exec'), proposed.__dict__, ns)
    private_step = MethodType(ns['step'], pol)
    pol._tail_capture = capture

    # Reproduce the two thin inherited step wrappers explicitly, retaining the
    # already-attached ForceSensing wrapper, whose saved step is rebound here.
    def inherited_step(o):
        n = o.n
        if n > 60 and n % 25 == 0:
            pol._identify(o.log, n)
        if n % 125 == 0:
            pol.id_hist.append((float(o.log[n-1, 0]), pol.tau_hat, pol.ndel_hat))
        if pol.mode == 1 or o.log[o.n-1, 5] < .3:
            pol._clear_motion()
        return private_step(o)

    # Fail if the host changes: this forwarding is exactly ServoID + IMM step.
    from v1_servo_id import ServoID
    from variants import MotionAwareIMM
    assert type(pol).step is ServoID.step
    assert 'return super().step(o)' in inspect.getsource(MotionAwareIMM.step)
    sensing = pol.step
    if isinstance(sensing, M.ForceSensing):
        sensing.step = inherited_step
    else:
        pol.step = inherited_step


def episode(case, run):
    sc = MC.scene_for(case['cell'], case['base'], case['variant'])
    pol = B.build('Huber-D4', O.chosen()['Huber-D4'], sc)
    records = []
    checks = []
    original_stiffness = pol.stiffness

    def stiffness(o, n):
        val, ok = original_stiffness(o, n)
        truth = getattr(o, 'truth_k', np.nan)
        if np.isfinite(truth):
            if run == 'kbar':
                val = np.float64(1.2*truth)
            elif run == 'kgood':
                pol.k_good = np.float64(1.2*truth)
        return val, ok

    if run in ('kbar', 'kgood'):
        pol.stiffness = stiffness

    def capture(v, command):
        o = v['o']
        truth = getattr(o, 'truth_k', np.nan)
        expected = v['sent']+v['du']/pol.kbar+(v['dx']-v['dq'])/v['k_un']
        assert command == expected or (np.isnan(command) and np.isnan(expected))
        if np.isfinite(truth):
            if run == 'kbar':
                assert pol.kbar == 1.2*truth
            elif run == 'kgood':
                assert pol.k_good == v['k_un'] == 1.2*truth
            elif run == 'ff_off':
                assert v['dq'] == 0.
            checks.append(True)
        records.append(dict(t=float(v['t']), k=float(truth),
                            kbar=float(pol.kbar), kgood=float(pol.k_good),
                            fb=float(v['du']/pol.kbar), ff=float(-v['dq']/v['k_un']),
                            exc=float(v['dx']/v['k_un']), command=float(command),
                            mode=int(pol.mode)))
        return command

    instrument_step(pol, run, capture)
    original_step = pol.step

    def step(o):
        # Enforce the intervention even during approach ticks (which do not
        # call stiffness); ordinary estimator updates still occur in contact.
        truth = getattr(o, 'truth_k', np.nan)
        if np.isfinite(truth):
            if run == 'kbar':
                pol.kbar = np.float64(1.2*truth)
            elif run == 'kgood':
                pol.k_good = np.float64(1.2*truth)
        return original_step(o)

    pol.step = step
    log = P3.episode(sc, pol, case['seed'], diagnostic=True)
    assert log is not None and checks
    m = M.metric(log, pol, sc)
    peak_i = int(np.argmax(np.where(log[:, 0] >= 1., log[:, 18], -np.inf)))
    tp = float(log[peak_i, 0])
    out = dict(case=case['label'], run=run, peak=m['peak'], L=float(np.log1p(m['J'])),
               peak_t=tp, command_checks=len(checks),
               zero_truth_command_ticks=sum(r['k'] == 0. for r in records),
               nonfinite_command_ticks=sum(not np.isfinite(r['command']) for r in records))
    if run == 'own':
        tol = .05 if case['rounded'] else 1e-7
        if abs(m['peak']-case['stored_peak']) > tol:
            raise RuntimeError(f"STOP own mismatch: {case}; actual={m['peak']}")
        # Raw local evidence, not episode/population statistics. Include enough
        # prehistory to see growth, contact loss and stiffness-channel jumps.
        near = [r for r in records if tp-1.0 <= r['t'] <= tp+.05]
        for r in near:
            j = int(np.argmin(abs(log[:, 0]-r['t'])))
            r.update(force=float(log[j, 18]), normal=float(log[j, 3]),
                     measured=float(log[j, 5]),
                     k_over_kbar=r['k']/r['kbar'],
                     k_over_kgood=r['k']/r['kgood'])
        out['around_peak'] = near
        w = (log[:, 0] >= tp-1.) & (log[:, 0] <= tp)
        out['force_window'] = log[w][:, [0, 18, 3, 6]].tolist()
    return out


def main():
    O.setup()
    selected = list(cases())
    assert selected[0]['seed'] == O.sensor_seed(selected[0]['cell'])
    print(json.dumps(dict(cases=selected)), flush=True)
    # All own episodes must reproduce before any intervention is interpreted.
    for run in ('own', 'kbar', 'kgood', 'ff_off'):
        for case in selected:
            # Exact zero truth is retained (no invented stiffness floor). IEEE
            # nonfinite commands invoke the plant's existing hold-command guard.
            with np.errstate(divide='ignore', invalid='ignore'):
                result = episode(case, run)
            print(json.dumps(result, allow_nan=True), flush=True)


if __name__ == '__main__':
    main()
