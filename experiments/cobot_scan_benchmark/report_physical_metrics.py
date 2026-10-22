"""Replay frozen controllers and report physical outcomes without retuning.

Run from any directory with ``python -X utf8 report_physical_metrics.py``.
Original seed-3000 and supplementary seed-4000 cohorts are kept separate.
All existing deployable baseline variants, three truth-stiffness rules, and
both candidate versions are evaluated. Per-scene fixed-grid oracle B2 is a
historical gate diagnostic, not an implementable controller, and is not rerun.
No original code, weights, tuning, or result files are overwritten.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

for _key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[_key] = "1"

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
OUT = ROOT / "physical_results"
CELLS = {"T3": {"ramp": 0.5, "A": 0.001}, "T1": {"ratio": 30.0, "A": 0.001}}


def read(name):
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def specifications(task):
    gate = next(c for c in read("gate_PC125.json") if c["task"] == task and c["dose"] == CELLS[task])
    ad = next(c for c in read("adaptive_results.json")["cells"] if c["task"] == task)
    mp = next(c for c in read("mpc_results.json")["cells"] if c["task"] == task)
    specs = [
        ("Fixed_low", "fixed", [0.3, 20.0, 0.0, 10.0, 0.0], -1),
        ("Fixed_high", "fixed", [0.3, 2500.0, 0.0, 10.0, 0.0], -1),
        ("Best_fixed", "fixed", gate["B1"]["gains"], -1),
    ]
    for name in ("EtR-A", "EtR-A-Haninger", "EtR-R", "Duan2018", "Song2025"):
        specs.append((name, "adaptive", ad["methods"][name]["tuning"], -1))
    for name in ("Jin_OS", "Jin_MS", "Jin_MS_literal", "Xue_A1", "Xue_nominal_ys0"):
        specs.append((name, "mpc", mp[name]["spec"], -1))
    for seed in range(3):
        specs.append(("Learned_VIC", "learned", {}, seed))
    specs.append(("Candidate1_v1", "c1", read("candidate1_results.json")["tuning"][task]["C1"]["config"], -1))
    # This is the exact constructor call used by v2cfg(), including the tuned v1
    # configuration. The short 'configs' record in the fresh JSON is overrides only.
    c1cfg = dict(read("candidate1_results.json")["tuning"][task]["C1"]["config"])
    c1cfg.update(k0=3e4)
    c1cfg.update(read("candidate1_v2_results.json")["sup_final"])
    g = gate["B1"]["gains"]
    c1cfg["fb"] = [g[0], g[1], g[3], g[4]]
    specs.append(("Candidate1_v2", "c1v2", c1cfg, -1))
    for seed in range(3):
        specs.append(("Candidate2", "c2", {}, seed))
    for latency in (0.0, 0.1, 0.2):
        specs.append((f"TruthK_{latency:.1f}s", "truth", [latency] + gate["B3b"][str(latency)]["rule"], -1))
    return specs


class TimedPolicy:
    def __init__(self, policy):
        self.policy = policy

    def reset(self, pl):
        self.policy.reset(pl)
        self.times = []

    def __call__(self, obs):
        started = time.perf_counter()
        gains = self.policy(obs)
        self.times.append(time.perf_counter() - started)
        return gains


def build(spec):
    name, kind, cfg, seed = spec
    if kind == "fixed":
        return tuple(cfg), None
    if kind == "adaptive":
        import baselines_adaptive as A
        if name == "Duan2018":
            return A.duan_gains(cfg["m"], cfg["b"], cfg["sigma"]), None
        if name == "Song2025":
            return A.Song(*[cfg["mapped"][k] for k in ("M", "D", "ki", "KPF")]), None
        q = cfg["rule_params"]
        if name == "EtR-A-Haninger":
            return A.EtR(zeta=q[0], M=q[1], cP=q[2], ki=q[3], Dmin=20.0, rule="haninger", **cfg["est_params"]), None
        return A.EtR(*q, rule="truthk", **cfg["est_params"]), None
    if kind == "mpc":
        import baselines_mpc as M
        return M.CLASSES[cfg[0]](**cfg[1]), M.UPDATE[cfg[0]]
    if kind == "learned":
        from baselines_learned import LearnedVIC, UPDATE
        return LearnedVIC(np.load(ROOT / f"learned_vic_seed{seed}.npz")["theta"]), UPDATE
    if kind == "c1":
        from candidate1_model import Candidate1
        p = Candidate1(**cfg)
        return p, p.update
    if kind == "c1v2":
        from candidate1_v2 import Candidate1v2
        p = Candidate1v2(**cfg)
        return p, p.update
    if kind == "c2":
        import candidate2_learned as C
        net, meta = C.load_net("full", seed)
        return C.Cand2Policy(net, **meta["cfg"]), None
    if kind == "truth":
        from baselines_ref import TruthK
        return TruthK(*cfg), None
    raise ValueError(kind)


def physical(log, duration):
    import scenarios as S
    if log is None:
        return dict(failed=1, rmse_N=None, max_abs_error_N=None, peak_force_N=None,
                    no_contact_s=None, longest_no_contact_s=None)
    f = log[log[:, 0] >= 1.0, 3]
    zero = f <= 0.0
    longest = run = 0
    for value in zero:
        run = run + 1 if value else 0
        longest = max(longest, run)
    m = S.metrics(log)
    return dict(failed=0, rmse_N=m["rmse"], max_abs_error_N=m["maxerr"],
                peak_force_N=m["peak"], no_contact_s=float(zero.sum() * 0.008),
                longest_no_contact_s=longest * 0.008,
                force_p95_abs_error_N=m["p95"], sampled_seconds=len(f) * 0.008,
                low_contact_s=float((f < 0.5).sum() * 0.008),
                time_above_7N_s=float((f > 7.0).sum() * 0.008),
                old_tuning_J=m["J"], duration_s=duration)


def replay(job):
    task, seed0, spec, n = job
    os.chdir(ROOT)
    import plant as P
    import scenarios as S
    name, kind, _, trainseed = spec
    # Warm JIT using training data, not the timed test episode.
    pp, _ = S.cohort(task, 1, 1000, CELLS[task])[0]
    ctrl, update = build(spec)
    if hasattr(ctrl, "bind"):
        ctrl.bind(pp)
    P.episode(pp, ctrl, T=0.2, seed=1000, update=update)
    rows = []
    for scene, (p, T) in enumerate(S.cohort(task, n, seed0, CELLS[task])):
        ctrl, update = build(spec)
        if hasattr(ctrl, "bind"):
            ctrl.bind(p)
        timed = TimedPolicy(ctrl) if callable(ctrl) else ctrl
        log = P.episode(p, timed, T=T, seed=seed0 + scene, update=update)
        row = dict(cohort=seed0, task=task, scene=scene, noise_seed=seed0 + scene,
                   method=name, training_seed=trainseed, **physical(log, T))
        tt = np.asarray(timed.times) * 1e3 if callable(ctrl) else np.array([])
        row.update(policy_median_ms=float(np.median(tt)) if len(tt) else None,
                   policy_p95_ms=float(np.percentile(tt, 95)) if len(tt) else None,
                   policy_max_ms=float(tt.max()) if len(tt) else None)
        rows.append(row)
    return rows


METRICS = ["rmse_N", "max_abs_error_N", "peak_force_N", "no_contact_s", "longest_no_contact_s",
           "force_p95_abs_error_N", "low_contact_s", "time_above_7N_s", "old_tuning_J",
           "policy_median_ms", "policy_p95_ms", "policy_max_ms"]


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def validate_rows(rows):
    expected = {(cohort, task, spec[0], spec[3], scene)
                for cohort in (3000, 4000) for task in CELLS
                for spec in specifications(task) for scene in range(40)}
    keys = [(r["cohort"], r["task"], r["method"], r["training_seed"], r["scene"]) for r in rows]
    assert len(keys) == len(set(keys)) and set(keys) == expected, "Missing or duplicated runs"
    for r in rows:
        if r["failed"]:
            assert r["rmse_N"] is None, "A failed scene must not become a finite score"
            continue
        assert all(np.isfinite(r[k]) for k in METRICS[:5])
        assert r["max_abs_error_N"] + 1e-10 >= r["rmse_N"] >= 0
        assert r["max_abs_error_N"] + 1e-10 >= abs(r["peak_force_N"] - 4.0)
        assert 0 <= r["longest_no_contact_s"] <= r["no_contact_s"] + 1e-10
        assert r["no_contact_s"] <= r["low_contact_s"] + 1e-10 <= r["sampled_seconds"] + 1e-10
        if r["no_contact_s"] > 0:
            assert r["max_abs_error_N"] >= 4.0 - 1e-10
    return dict(expected_runs=len(expected), observed_runs=len(rows), unique_runs=len(set(keys)),
                failed_runs=sum(r["failed"] for r in rows), passed=True)


def aggregate(rows):
    out = []
    for cohort, task, method in sorted({(r["cohort"], r["task"], r["method"]) for r in rows}):
        rr = [r for r in rows if (r["cohort"], r["task"], r["method"]) == (cohort, task, method)]
        seeds = sorted({r["training_seed"] for r in rr})
        scenes = sorted({r["scene"] for r in rr})
        row = dict(cohort=cohort, task=task, method=method, n_scenes=len(scenes),
                   n_training_seeds=len(seeds) if seeds != [-1] else 0,
                   failed_runs=sum(r["failed"] for r in rr))
        for metric in METRICS:
            # Average training seeds within each scene before describing scene
            # variation. Seeds are not additional independent physical trials.
            values = []
            for scene in scenes:
                vv = [r.get(metric) for r in rr if r["scene"] == scene]
                values.append(float(np.mean(vv)) if all(v is not None for v in vv) else np.nan)
            valid = np.isfinite(values)
            row[metric + "_mean"] = float(np.mean(values)) if all(valid) else None
            row[metric + "_sd"] = float(np.std(values, ddof=1)) if len(values) > 1 and all(valid) else None
        pk = [r["peak_force_N"] for r in rr]
        row["worst_peak_force_N"] = max(pk) if all(v is not None for v in pk) else None
        out.append(row)
    return out


def replay_checks(rows):
    """Compare replay outcomes with the saved physical outcomes, never timing."""
    out = []
    for cohort, task, name, seed in sorted({(r["cohort"], r["task"], r["method"], r["training_seed"]) for r in rows}):
        old = None
        if cohort == 4000:
            names = {"Candidate1_v1": "C1v1", "Candidate1_v2": "C1v2", "Best_fixed": "B1", "EtR-A": "EtR-A", "TruthK_0.1s": "B3b_0.1"}
            if name in names:
                old = read("candidate1_v2_results.json")["fresh"][task]["methods"][names[name]]["summary"]
        else:
            gate = next(c for c in read("gate_PC125.json") if c["task"] == task and c["dose"] == CELLS[task])
            if name in ("Fixed_low", "Fixed_high", "Best_fixed"):
                key = {"Fixed_low": "B0_low", "Fixed_high": "B0_high", "Best_fixed": "B1"}[name]
                old = gate[key]["test"] if key == "B1" else gate[key]
            elif name.startswith("TruthK_"):
                old = gate["B3b"][name.split("_")[1][:-1]]["test"]
            elif name.startswith("EtR") or name in ("Duan2018", "Song2025"):
                old = next(c for c in read("adaptive_results.json")["cells"] if c["task"] == task)["methods"][name]["test"]
            elif name.startswith("Jin") or name.startswith("Xue"):
                old = next(c for c in read("mpc_results.json")["cells"] if c["task"] == task)[name]["test"]
            elif name == "Learned_VIC":
                old = next(c for c in read("learned_results.json")["cells"] if c["task"] == task)["learned"][str(seed)]["test"]
            elif name == "Candidate1_v1":
                old = read("candidate1_results.json")["test"][task]["methods"]["C1"]["test"]
            elif name == "Candidate1_v2":
                old = read("candidate1_v2_results.json")["oldtest"][task]["methods"]["C1v2"]["summary"]
            elif name == "Candidate2":
                old = next(c for c in read("candidate2_results.json")["cells"].values()
                           if c["task"] == task and c["dose"] == CELLS[task])["methods"]["cand2_full"]["seeds"][str(seed)]["test"]
        if old is None:
            continue
        rr = [r for r in rows if (r["cohort"], r["task"], r["method"], r["training_seed"]) == (cohort, task, name, seed)]
        for metric, key in (("rmse_N", "rmse"), ("max_abs_error_N", "maxerr"), ("peak_force_N", "peak"), ("old_tuning_J", "J")):
            if key not in old:
                continue
            value = float(np.mean([r[metric] for r in rr])) if all(r.get(metric) is not None for r in rr) else None
            delta = value - old[key] if value is not None else None
            out.append(dict(cohort=cohort, task=task, method=name, training_seed=seed, metric=metric,
                            archived=old[key], replay=value, difference=delta,
                            within_numeric_tolerance=value is not None and abs(delta) <= 1e-8 + 1e-5 * abs(old[key])))
    return out


def paired(rows):
    out = []
    rng = np.random.default_rng(2261022)
    for cohort, task in sorted({(r["cohort"], r["task"]) for r in rows}):
        rr = [r for r in rows if r["cohort"] == cohort and r["task"] == task]
        methods = sorted({r["method"] for r in rr})
        data = {}
        for method in methods:
            rm = [r for r in rr if r["method"] == method]
            data[method] = {metric: np.array([np.mean([r[metric] for r in rm if r["scene"] == scene])
                                           for scene in range(40)]) for metric in METRICS[:5]
                            if all(r[metric] is not None for r in rm)}
        for candidate in ("Candidate1_v1", "Candidate1_v2", "Candidate2"):
            for reference in methods:
                if reference.startswith("Candidate"):
                    continue
                for metric in METRICS[:5]:
                    if metric not in data[candidate] or metric not in data[reference]:
                        continue
                    diff = data[candidate][metric] - data[reference][metric]
                    boot = diff[rng.integers(0, len(diff), size=(10000, len(diff)))].mean(axis=1)
                    out.append(dict(cohort=cohort, task=task, candidate=candidate, reference=reference,
                                    metric=metric, difference_mean=float(diff.mean()),
                                    ci95_low=float(np.percentile(boot, 2.5)), ci95_high=float(np.percentile(boot, 97.5)),
                                    n_scenes=len(diff), correction="none; exploratory marginal interval"))
    return out


ORDER = [s[0] for s in specifications("T3")]
ORDER = list(dict.fromkeys(ORDER))
LABELS = {
    "Fixed_low": "Fixed low damping", "Fixed_high": "Fixed high damping",
    "Best_fixed": "Best tuned fixed (B1)", "EtR-A": "Estimate + stiffness rule",
    "EtR-A-Haninger": "Estimate + Haninger damping rule", "EtR-R": "RLS estimate + rule",
    "Duan2018": "Duan (RAS 2018), adapted", "Song2025": "Song (RA-L 2025), adapted",
    "Jin_OS": "Jin (TMECH 2023), one-step", "Jin_MS": "Jin (TMECH 2023), multistep",
    "Jin_MS_literal": "Jin multistep, unanchored diagnostic",
    "Xue_A1": "Xue MPVIC, initial-surface prior", "Xue_nominal_ys0": "Xue MPVIC, nominal initial surface",
    "Learned_VIC": "E27-inspired learned VIC (3 seeds)", "Candidate1_v1": "Candidate 1 v1",
    "Candidate1_v2": "Candidate 1 v2", "Candidate2": "Candidate 2 (3 seeds)",
    "TruthK_0.0s": "Truth-stiffness reference, 0 s delay",
    "TruthK_0.1s": "Truth-stiffness reference, 0.1 s delay",
    "TruthK_0.2s": "Truth-stiffness reference, 0.2 s delay",
}


def render(summary, protocol, checks):
    cols = ["rmse_N", "max_abs_error_N", "peak_force_N", "no_contact_s", "worst_peak_force_N"]
    lines = ["# Physical outcomes of the cobot scanning benchmark", "",
             "The PC125 position-controlled plant runs at 125 Hz with a 4 N force target, "
             "one control-cycle transport delay, a 19 ms servo lag, noisy 1 kHz force sensing "
             "and a 30 Hz filter. T1 scans a stiffness ladder (modulus ratio 30); T3 scans "
             "bone-like ridges with 0.5 s push/retraction ramps. Both main tasks use 1 mm "
             "respiration. Each cohort has 40 paired scenes per task. All controller parameters "
             "and neural weights are frozen from the existing training records; no retuning is performed.", "",
             "Seed 3000 is the original, already inspected test set. Seed 4000 is the previously "
             "used v2 evaluation set, extended here to all controllers. Neither is described as a "
             "new untouched test set. Three neural training seeds are averaged within each scene; "
             "the independent sample size remains 40, not 120. Worst sampled peaks for neural "
             "methods cover all 120 scene/seed runs and are descriptive, not an equal-run-count "
             "comparison against a single fitted controller.", "",
             "Force outcomes use the true simulated force after the original 1 s warm-up. "
             "RMSE, maximum absolute error and peak force are computed per episode, then averaged "
             "over scenes. The worst-peak column additionally gives the largest sampled peak across all "
             "scenes/seeds. Zero-force duration counts F <= 0 as an unloaded-contact proxy, without the old 0.5 N operational "
             "cutoff. Failed integrations remain in the data and give NA physical means, rather "
             "than silently dropping scenes. Samples are 8 ms apart, so all peaks and durations are sampled quantities. "
             "No clinical safety threshold is asserted. Units are N and s.", ""]
    for cohort in protocol["cohorts"]:
        for task in CELLS:
            lines += [f"## {task}, cohort seed {cohort}", "",
                      "| Controller | Force RMSE (N) | Mean max absolute error (N) | Mean peak force (N) | Mean zero-force time (s) | Worst sampled peak (N) | Failed runs |",
                      "|---|---:|---:|---:|---:|---:|---:|"]
            for method in ORDER:
                row = next(r for r in summary if (r["cohort"], r["task"], r["method"]) == (cohort, task, method))
                vv = [row[c if c == "worst_peak_force_N" else c + "_mean"] for c in cols]
                values = [f"{v:.4f}" if v is not None else "NA" for v in vv]
                lines.append("| " + LABELS[method] + " | " + " | ".join(values) + f" | {row['failed_runs']} |")
            lines.append("")
    lines += ["## Policy-update timing", "",
              "Timing includes each policy's estimator, feature extraction and decision call. "
              "It excludes shared plant/servo execution, communication, and JIT compilation. "
              "Warm-up uses training scenes. Measurements are host prototype timings during a "
              "two-process replay, not robot real-time certification. Each table entry averages "
              "per-episode medians or 95th percentiles; it is not the pooled call-distribution percentile. "
              "A dash means a fixed parameter vector, not a claim of zero control computation.", "",
              "| Controller | Decision period (ms) | T1 mean episode-median (ms) | T1 mean episode-p95 (ms) | T3 mean episode-median (ms) | T3 mean episode-p95 (ms) |",
              "|---|---:|---:|---:|---:|---:|"]
    for method in ORDER:
        vals = []
        for task in ("T1", "T3"):
            row = next(r for r in summary if (r["cohort"], r["task"], r["method"]) == (3000, task, method))
            vals += [row["policy_median_ms_mean"], row["policy_p95_ms_mean"]]
        period = "—" if method in ("Fixed_low", "Fixed_high", "Best_fixed", "Duan2018") else "48" if method == "Learned_VIC" else "16" if method.startswith("Jin_") else "8"
        lines.append("| " + LABELS[method] + " | " + period + " | " + " | ".join(f"{v:.5f}" if v is not None else "—" for v in vals) + " |")
    bad = [r for r in checks if not r["within_numeric_tolerance"]]
    groups = sorted({(r["cohort"], r["task"], r["method"]) for r in bad})
    lines += ["", "## Replay agreement", "",
              "Original files are unchanged. Physical outcomes above are newly replayed values. "
              "The comparison tolerance is abs(error) <= 1e-8 + 1e-5 abs(archived); timing is "
              "not compared across machines. See replay_checks.csv for every check. "
              f"{len(bad)} of {len(checks)} archived-metric checks differ beyond this strict tolerance, "
              "including unstable Jin_MS_literal trajectories and the failed low-damping extreme. "
              "This is not a claim of exact cross-environment reproduction. Small differences "
              "also occur in selected Jin, Xue and v1 outputs. The cause of the unstable "
              "variant's larger discrepancy has not been isolated.", "",
              "Flagged cohort/task/controller groups: " + "; ".join(f"{c}/{t}/{m}" for c, t, m in groups) + "."]
    lookup = {(r["cohort"], r["task"], r["method"]): r for r in summary}
    lines += ["", "## Observations on the supplementary seed-4000 cohort", ""]
    for task in ("T1", "T3"):
        b = lookup[4000, task, "Best_fixed"]
        c = lookup[4000, task, "Candidate1_v2"]
        n = lookup[4000, task, "Candidate2"]
        lines.append(f"- {task}: best-fixed mean force RMSE {b['rmse_N_mean']:.4f} N, "
                     f"Candidate1_v2 {c['rmse_N_mean']:.4f} N and Candidate2 {n['rmse_N_mean']:.4f} N. "
                     f"The corresponding mean episode-max absolute errors are {b['max_abs_error_N_mean']:.4f}, "
                     f"{c['max_abs_error_N_mean']:.4f} and {n['max_abs_error_N_mean']:.4f} N.")
    v1 = lookup[4000, "T3", "Candidate1_v1"]
    v2 = lookup[4000, "T3", "Candidate1_v2"]
    h = lookup[4000, "T3", "EtR-A-Haninger"]
    lines += [f"- T3 v1's worst sampled force is {v1['worst_peak_force_N']:.2f} N; v2 reduces it to "
              f"{v2['worst_peak_force_N']:.2f} N. The Haninger-rule comparator has a lower worst sampled peak "
              f"({h['worst_peak_force_N']:.2f} N) despite its larger force RMSE. Mean tracking gains do not imply "
              "dominance on every outcome.",
              "- The 0.25 s push and 4 mm respiration stress cells are outside these two main-cell "
              "replays. Their existing negative findings remain in candidate1_v2_log.md and are not "
              "replaced by these positive main-task results."]
    lines += ["", "## Scope and reproducibility", "",
              "These are physical performance tables, not the old weighted tuning J. The original "
              "J is retained only in per-scene source data for replay checks. The exact scenario, "
              "configuration, source and weight hashes are in protocol.json. Source data and "
              "scene-level SDs are in per_scene.csv and summary.csv. Seed-specific neural "
              "summaries are in training_seed_summary.csv. Pairwise scene-bootstrap intervals "
              "are in paired_differences.csv, without multiplicity correction; no significance "
              "stars or universal superiority claim is attached.", "",
              "Jin and Xue rows retain the existing parameter searches and adaptation choices. "
              "In particular, Xue has kmax=0 and no train-feasible grid setting; those rows are "
              "not evidence that a fully tuned Xue family is inferior. Learned_VIC uses the "
              "existing E27-inspired ES adaptation, with unchanged weights. TruthK rows use "
              "unavailable true stiffness and are information references, not deployable baselines "
              "or globally optimal controllers. The per-scene fixed-grid B2 and isolated "
              "preview/feedforward oracle probes remain historical gate diagnostics; their "
              "missing physical outputs are not fabricated in this report.", "",
              "The lateral axis is position-servoed and not adapted by these policies. No claim "
              "of improved lateral tracking, complete B0 validation or certified near-optimality "
              "follows from these normal-force tables.", "",
              "Reproduce: `python -X utf8 experiments/cobot_scan_benchmark/report_physical_metrics.py`."]
    (OUT / "results.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    os.chdir(ROOT)
    OUT.mkdir(exist_ok=True)
    source_files = sorted([*ROOT.glob("*.py"), *ROOT.glob("*.npz"),
                           *[ROOT / n for n in ("adaptive_results.json", "mpc_results.json", "learned_results.json", "gate_PC125.json",
                                               "candidate1_results.json", "candidate1_v2_results.json", "candidate2_results.json")]])
    import numba
    import scipy
    protocol = dict(cohorts=[3000, 4000], scenes_per_task=40, tasks=CELLS,
                    python=sys.version, system=platform.platform(), numpy=np.__version__,
                    numba=numba.__version__, scipy=scipy.__version__, processor=platform.processor(),
                    git_head=subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
                    workers=args.workers, warmup_s=1.0, reporting_Tc_s=0.008,
                    source_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files},
                    controller_specifications={task: specifications(task) for task in CELLS})
    previous = OUT / "protocol.json"
    if args.resume and previous.exists():
        old = json.loads(previous.read_text())
        if json.dumps(old["controller_specifications"], sort_keys=True) != json.dumps(protocol["controller_specifications"], sort_keys=True):
            raise RuntimeError("Frozen controller specifications changed; do not reuse the replay cache")
        for name, digest in old["source_sha256"].items():
            if name != Path(__file__).name and protocol["source_sha256"].get(name) != digest:
                raise RuntimeError(f"Input {name} changed; do not reuse the replay cache")
    previous.write_text(json.dumps(protocol, indent=2), encoding="utf-8")
    rows = []
    chunks = ROOT.parents[1] / ".sync-runtime" / "physical_reporting"
    chunks.mkdir(exist_ok=True, parents=True)
    jobs = []
    for cohort in protocol["cohorts"]:
        for task in CELLS:
            for spec in specifications(task):
                cache = chunks / f"{cohort}_{task}_{spec[0]}_seed{spec[3]}.json"
                if args.resume and cache.exists():
                    rows += json.loads(cache.read_text())
                else:
                    jobs.append(((task, cohort, spec, 40), cache))
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        pending = {pool.submit(replay, job): (job, path) for job, path in jobs}
        done = 0
        for future in as_completed(pending):
            job, path = pending[future]
            result = future.result()
            path.write_text(json.dumps(result, indent=1, allow_nan=False), encoding="utf-8")
            rows += result
            done += 1
            good = [r["rmse_N"] for r in result if r["rmse_N"] is not None]
            print(f"[{done}/{len(jobs)}] cohort={job[1]} {job[0]} {job[2][0]} seed={job[2][3]} "
                  f"mean RMSE={np.mean(good):.6f} N failures={sum(r['failed'] for r in result)}", flush=True)
    rows.sort(key=lambda r: (r["cohort"], r["task"], r["method"], r["training_seed"], r["scene"]))
    for row in rows:
        for metric in METRICS:
            row.setdefault(metric, 1000.0 if row["failed"] and metric == "old_tuning_J" else None)
    validation = validate_rows(rows)
    (OUT / "validation.json").write_text(json.dumps(validation, indent=2), encoding="utf-8")
    write_csv(OUT / "per_scene.csv", rows)
    summary = aggregate(rows)
    write_csv(OUT / "summary.csv", summary)
    checks = replay_checks(rows)
    write_csv(OUT / "replay_checks.csv", checks)
    seed_summary = []
    for seed in range(3):
        for r in aggregate([r for r in rows if r["training_seed"] == seed]):
            seed_summary.append(dict(training_seed=seed, **r))
    write_csv(OUT / "training_seed_summary.csv", seed_summary)
    write_csv(OUT / "paired_differences.csv", paired(rows))
    render(summary, protocol, checks)
    print(f"Saved {len(rows)} episode records and {len(summary)} summary rows to {OUT}", flush=True)


if __name__ == "__main__":
    main()
