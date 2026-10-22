# Cobot scanning benchmark: task protocol and headroom gate

**Physical-outcome tables (2026-10-22):** [complete report](physical_results/results.md),
[per-scene source data](physical_results/per_scene.csv), and
[replay code](report_physical_metrics.py). All frozen deployable baseline variants,
both Candidate 1 versions, Candidate 2 and truth-stiffness references were replayed
on the same T1/T3 seed-3000 and seed-4000 cohorts, 40 scenes per task and cohort.
The report separates force RMSE, maximum error, peak force, actual no-contact
duration and policy-update timing. It retains failed scenes and documents replay
differences without changing the original result files or retuning controllers.

Purpose: a simulation task set in which the research question can be tested
fairly. The question is whether newly acquired contact information can promptly
update near-optimal interaction parameters with real control benefit. Every
baseline must be reasonable and documented, so no task should make adaptation
look good through an incomplete fixed controller or contrived conditions.

The protocol came out of the setup reverse analysis
(`research/simulation_setup_reverse_analysis.md`). Three independent design
proposals, each with its own probes, were judged by three reviewers (hostile
reviewer, physics, question fit) and then merged. This directory implements the
merged protocol and its first gate.

## What the design round established (probes on 3-8 in-sample scenes; indicative)

- **The fixed controller must be complete.** That means admittance (M, D, K),
  a force integral ki and a direct proportional force path KPF, all tuned
  jointly. Each of the three proposals had left out one term, and that omission
  produced most of their reported headroom. Examples:
  - 125 Hz homogeneous tissue: 0.49 N drops to 0.16 N once the term is added.
  - Stiffness ladder at r = 30: 0.72 N drops to 0.26 N.
- **With the complete family, the loophole holds on fast platforms.** A
  jointly tuned fixed controller is near-optimal on a 500 Hz cobot and on the
  1 kHz torque plant: 0.035-0.039 N fixed vs 0.028-0.037 N oracle on
  homogeneous tissue.
- **The 125 Hz cobot class mostly closes the gap.** On homogeneous tissue and
  on bone-like ladders, fixed and oracle are both inside the 0.5 N band (fixed
  0.08-0.26 N vs oracle 0.03-0.14 N).
- **Only one failure-level gap survives.** It needs the 125 Hz class, a
  bone-like region in the scene, and fast patient motion together:
  - Fixed: 0.69-0.94 N RMSE, peaks 6.9-8.3 N.
  - Truth-stiffness schedule: 0.16-0.21 N, peaks about 5 N.
  - The gap shrinks with slower motion: at a 1 s push ramp it is 0.18-0.23 N
    vs 0.06-0.09 N.
  - It persists with 1 mm breathing, so it is not a breathing-amplitude artefact.
- **The mechanism.** Bone forces a conservative fixed tuning, and that tuning
  then cannot follow fast surface motion on soft tissue.
- **Excluded from prompt-update claims.** First-contact landing is set by
  approach speed before any contact information exists, so it is not a
  prompt-update claim. Friction has no control consequence on a
  position-servoed lateral axis, so shear is reported only.

## Reference platform (`plant.py`)

- **PC125.** Position-controlled cobot of the UR5/xArm6 class at 125 Hz
  (F07:353-357, F08:168-172).
  - One control period of transport delay.
  - First-order servo lag of 19 ms (UR5e TPI 0.1 per 2 ms, F09:353-362).
    Carrying it over to 125 Hz is an assumption.
  - F/T sensor sampled at 1 kHz: sigma 0.08 N, bias random walk, 30 Hz low-pass.
- **Brackets.**
  - PC125dead: mostly dead time instead of a lag. The latency structure must be
    identified on hardware before any claim.
  - PC500: UR5e at 500 Hz, used as a negative control.
- **Controller, shared by every method.**
  - Admittance: M yc'' + D yc' + K (yc - y_ref) = ef + ki ∫ef.
  - Command: y_cmd = yc + KPF ef.
  - Discretisation: exact zero-order hold, with velocity cap and anti-windup.
  - Bounds: M 0.25-3 kg, D 20-2500 N s/m, K 0-1000 N/m, ki 0-30 1/s,
    KPF 0-1e-3 m/N.
- **Contact.** B0 cylindrical finite footprint with SLS relaxation,
  Hunt-Crossley damping and Stribeck friction.
  - Bone-like ridges with partial aperture coverage and an optional bump.
  - Breathing, plus patient push/retract events.
  - Contact stiffness at 4 N: 30 kPa gives 1.15 N/mm, 1 MPa gives 18.4 N/mm.
- **Sources for latency.** Web sources were reached only through search
  excerpts (marked [x]); local texts were read in full.
  - UR servoj: a lag set mainly by lookahead_time (0.03-0.2 s; 0.01 s from SW
    5.26), whose larger values "introduce a higher delay" [x].
  - UR10e forum measurement: servoj at 2 ms with lookahead 0.03 gives about
    34 ms of cross-correlation delay; sync() gives 4-6 ms (UR Client Library
    issue #86).
  - Wang 2024: about 20 ms, first-order in form (F09:353-362).
  - Mathiassen 2016 (UR5 ultrasound): 125 Hz loop, about 16.6 Hz force
    bandwidth [x].
  - No servo identification of xArm was found.
  - Nominal used here: 1 cycle of dead time plus a 19 ms lag, which is at the
    low (conservative) end of these sources.
  - Sweep: dead time 1-2 cycles, lag 10-60 ms.
- **Kinematic arm.** The stiff servo does not deflect under contact force.
- **Policy interface.** A policy may read time, filtered force, commanded and
  encoder positions. Truth columns (contact stiffness) are for oracles only.
- **Cost.** About 0.015 s per 5 s episode.

## Tasks (`scenarios.py`)

| Task | Role | Construction (synthetic, dose variables in bold) |
|---|---|---|
| T1 ladder | Benchmark, timeliness instrument | Soft base 10-40 kPa, 1-2 ridges of **ratio** 10-30 (width 8-25 mm, edge 0.5-2 mm), **speed** 10/20/30 mm/s, **breathing A** 2-5 mm, 8 s |
| T2 push | Benchmark, non-inferiority / false adaptation | Homogeneous soft (10-40 kPa) or stiff (0.5-1.5 MPa), 1-3 push/retract events of **amp** 3-8 mm with **ramp** 0.15-0.5 s, 8 s |
| T3 spine | Headline scanning task | Soft 10-40 kPa, spinous-process ridges at 25-35 mm pitch, width 8-20 mm, **ratio** 10-30, speed 10-20 mm/s, breathing 2-5 mm, 1-2 patient events (**amp**, **ramp**) at random times, 10 s |
| N1 homog | Negative control | Homogeneous 10 kPa-2 MPa with breathing, 6 s |

Sources for patient motion:
- **Breathing.**
  - Chest surface 4-12 mm at 0.2-0.34 Hz (Shafiq & Veluvolu 2014 [x]).
  - Abdominal surface AP 2.5-13.9 mm (Fahmi 2018 [x]).
  - Spine target in prone position AP 3.4 ± 1.6 mm (Mizonobe 2023 [x]).
  - Moreira's 0.2 Hz / 5 mm test case (F10:875-880).
- **Cough.** The contraction phase lasts 0.5 ± 0.2 s; abdominal-wall
  displacement is up to 18 ± 8 mm (Jourdan 2022, dynamic MRI [x]).
  - So our 0.5 s ramp has a source, and the 0.25 s ramp is a fast case.
- **Not sourced.** Flinch or posture-shift speeds were not found; they are
  engineering choices.

Unsourced choices are treated as dose variables and swept:
- ridge geometry and the 10-30x modulus ratio (they come from the research
  question, not from anatomy);
- patient-motion amplitude and speed.

Splits are by seed: train 1000, dev 2000, test 3000.

## Metrics (`scenarios.metrics`)

All metrics use the true contact force after a 1 s warm-up.

- **Accuracy:** RMSE, 95th percentile and maximum |F - Fd|, time inside ±0.5 N.
- **Safety:** peak force, time above Fd + 3 N, contact-loss fraction and
  longest loss, oscillation RMS above 10 Hz, peak shear.
- **Tuning composite.** Declared in advance and used only for tuning:
  J = mean(((F - Fd)/0.5)^2) + 20 loss + mean(max(0, F - 7)^2).
- **Claims.** A claim counts only in physical units. Paired comparisons use
  bootstrap CIs.

## Baselines (all share plant, seeds, score, parameter vector and tuning budget)

- **B0 fixed extremes:** D at its bounds (20 and 2500 N s/m) with ki 10.
- **B1 best single fixed:**
  - Grid M {0.3, 1, 3} × D (12 log steps) × K {0, 300} × ki {0, 3, 10, 30} ×
    KPF {0, 2e-5, 5e-5, 1e-4, 2e-4, 4e-4}.
  - Then Nelder-Mead in log space from the 5 best points.
  - Selection: minimum training J subject to peak ≤ 7 N in at least 95 % of
    training scenes.
- **B2 per-scene oracle:** best grid point on each test scene.
- **B3b truth-stiffness schedule (oracle):** D = max(Dmin, k/wI), KPF = cP/k.
  - Here k is the true footprint stiffness, delayed by 0 / 0.1 / 0.2 s from the
    moment it starts to change under the footprint.
  - The rule parameters are tuned on train.
- **To be added before candidates:** B1-robust (tuned over latency and filter
  mismatch), estimate-then-rule (EKF / lock-in stiffness estimate + rule,
  Haninger D ≥ 2ζ sqrt(M k)), Duan/E04 adaptive damping, Jin E41 MPC-VIC, Xue
  MPVIC, an E27-style learned VIC, and breathing feedforward for all methods.

## Headroom gate (`gate.py`)

Six cells on PC125, 20 train / 40 test scenes each:
- T3 with push ramp {0.25, 0.5} s × breathing {1, 4} mm;
- T1 at ratio 30 × breathing {1, 4} mm.

Pass rule:
- B1 fails at least one physical criterion: peak above 7 N in at least 20 % of
  scenes, mean RMSE above 0.5 N, or time inside ±0.5 N below 70 %;
- B3b at 0.1 s delay meets all of them;
- the paired J-difference CI excludes zero.

Results are in `gate_PC125.json` and `gate_PC125.log` and are summarised
below once complete.

## Results (PC125, 40 test scenes per cell; RMSE N / mean peak N / fraction of scenes with peak > 7 N)

| Cell | B1 best fixed | B2 per-scene fixed | B3b truth-k, 0 s | B3b, 0.1 s | B3b, 0.2 s | J: B1 → B3b 0.1 s (paired 95 % CI of difference) |
|---|---|---|---|---|---|---|
| T3, push ramp 0.25 s, breathing 1 mm | 0.67 / 7.8 / 0.65 | 0.50 / 7.0 / 0.40 | 0.39 / 6.9 / 0.30 | 0.41 / 6.9 / 0.30 | 0.78 / 8.8 / 0.47 | 2.27 → 1.08 (-1.67, -0.79) |
| T3, 0.25 s, 4 mm | 0.74 / 8.1 / 0.68 | 0.51 / 7.0 / 0.38 | 0.39 / 6.9 / 0.30 | 0.41 / 7.0 / 0.30 | 0.71 / 9.6 / 0.53 | 2.93 → 1.09 (-2.89, -1.13) |
| T3, 0.5 s (cough contraction), 1 mm | 0.31 / 5.5 / 0.03 | 0.19 / 4.9 / 0.00 | 0.15 / 4.8 / 0.00 | 0.15 / 4.8 / 0.00 | 0.30 / 5.3 / 0.03 | 0.57 → 0.12 (-0.91, -0.20) |
| T3, 0.5 s, 4 mm | 0.32 / 5.5 / 0.03 | 0.19 / 4.9 / 0.00 | 0.15 / 4.8 / 0.00 | 0.16 / 4.8 / 0.00 | 0.32 / 5.3 / 0.03 | 0.60 → 0.13 (-0.95, -0.21) |
| T1 ladder r = 30, 1 mm | 0.16 / 4.8 / 0.03 | 0.11 / 4.6 / 0.00 | 0.07 / 4.3 / 0.00 | 0.11 / 4.4 / 0.00 | 0.13 / 4.7 / 0.05 | 0.25 → 0.06 (-0.50, -0.02) |
| T1 ladder r = 30, 4 mm | 0.18 / 4.8 / 0.03 | 0.13 / 4.6 / 0.00 | 0.08 / 4.4 / 0.00 | 0.13 / 4.5 / 0.00 | 0.16 / 4.8 / 0.07 | 0.27 → 0.09 (-0.51, -0.01) |

B1 gains are in `gate_PC125.json`. All six cells settle on M ≈ 1-1.7 kg, D 73-180 N s/m,
ki up to its bound, and KPF ≈ 1e-4 m/N.

The strict gate fails in all six cells, for one of two reasons:
- With the sourced 0.5 s ramp and on the ladders, B1 already meets every physical
  criterion.
- With the fast 0.25 s ramp, B3b also has 30 % of scenes above 7 N. Most of
  those peaks occur where bone edges and pushes coincide.

What the gate shows:
1. **Known contact stiffness has real but bounded value.** With the
   gains-only oracle, J falls by a factor of 2.7-4.7 and every CI excludes
   zero. In physical units this halves force RMSE (0.31 → 0.15 N at cough
   speed) and clears peak violations. At sourced patient-motion speed, however,
   both controllers are already inside the 0.5 N band. Failure-level headroom
   appears only with the faster, unsourced 0.25 s pushes, where B1 has 65-68 %
   of scenes above 7 N against 30 % for the oracle.
2. **Timeliness decides it.** With 0.2 s of information delay, the same oracle
   is no better than the best fixed controller, and is worse at the 0.25 s ramp.
   The usable latency budget is about 0.1 s.
3. **On ladders, time variation of the gains adds little.** The per-scene
   fixed controller already matches the oracle at 0.1 s. A stiffness-map
   feedforward on the command is the stronger lever: an extra probe
   (`baselines_ref.TruthKFF`, 40 test scenes) lowered T1 RMSE from
   0.106 N (gains only) to 0.055 N with 0.05 s preview, against 0.155 N for B1.
   The same feedforward does not change the T3 push peaks.
4. **Spatial preview of the gains alone does not help.** Switching to stiff
   gains earlier only hurts soft-tissue tracking: T3 RMSE 0.39 → 0.53 N with
   0.4 s preview.

These are oracle upper bounds; causal methods will recover only part of them.
For the paper, a gains-only claim is bounded by result 1. The data point
towards updating the gains and the command reference together from contact
information (result 3), with a latency budget of about 0.1 s (result 2).

## Reproduce

`python gate.py PC125` (all cells) or `python gate.py PC125 0,1` (selected
cells). Requires numpy, scipy and numba.

## Baselines on the two sourced cells (40 test scenes; J / RMSE N / fraction of scenes with peak > 7 N)

| Method | T3 spine, cough-speed push (0.5 s), 1 mm | T1 ladder r = 30, 1 mm | Notes |
|---|---|---|---|
| B1 best fixed (complete family) | 0.570 / 0.31 / 0.03 | 0.247 / 0.16 / 0.03 | `gate.py` |
| B3b truth-k oracle, 0.1 s | 0.120 / 0.15 / 0 | 0.062 / 0.11 / 0 | upper reference, not causal |
| EtR-A: lock-in stiffness estimate + rule | 0.407 / 0.27 / 0.10 | 0.077 / 0.13 / 0 | `baselines_adaptive.py`; ΔJ vs B1 CI includes 0 in both cells |
| EtR-A with Haninger damping rule | 0.410 / 0.30 / 0.05 | 0.174 / 0.16 / 0.08 | |
| Duan RAS2018 adaptive damping | 6.20 / 1.06 / 0.68 | 2.59 / 0.62 / 0.38 | exactly a fixed PI admittance with K = KPF = 0 |
| Song RA-L2025 NAIPC (mapped) | 0.99 / 0.44 / 0.03 | 0.55 / 0.32 / 0.03 | |
| Jin TMECH2023 MPC-VIC, one-step | 0.549 / 0.32 / 0.03 | 0.200 / 0.15 / 0.03 | `baselines_mpc.py`; anchored to B1 gains |
| Jin MPC-VIC, multi-step | 0.604 / 0.35 / 0.05 | 0.123 / 0.14 / 0.03 | literal (unanchored) version: 16.7 / 4.06 J |
| Xue RAS2025 MPVIC | 7.76 / 1.28 / 0.85 | 1.36 / 0.55 / 0.20 | HC model undefined once pushes exceed the assumed surface |
| E27-style learned VIC (ES, 3 seeds) | 0.47 ± 0.07 J, 0.32 N | 0.57 ± 0.39 J, 0.18 N | `baselines_learned.py`; learned near-constant gains |

What the baselines show:
- These particular field-method adaptations do not establish consistent superiority
  over the complete fixed controller on the two cells. This concerns the mapped
  implementations; it does not establish why the original published results differ.
- A simple causal estimate-then-rule (dither-based stiffness estimate) recovers most
  of the oracle value on the ladder (J 0.077 vs oracle 0.062), but only a small part on
  the spine task with patient motion (0.41 vs 0.12). Its estimate lags by about
  70-100 ms at ridge entry, and patient motion corrupts it.
- So the remaining room for a method is causal, prompt use of contact information
  while the tissue moves. That is the T3 gap between EtR-A and the oracle.

Fidelity caveats for each mapping are in each module's docstring and in the agents'
reports summarised in the JSON files:
- Jin was anchored to B1, with an estimator gain correction.
- Xue's numerical model limited its gain grid.
- The learned policy used evolution strategies instead of DRL.

## Candidate methods on the same cells (40 test scenes; J / RMSE N / fraction of scenes with peak > 7 N)

Both candidates use only causal signals (force, encoder, command, scan position), are tuned or trained on
the train split only, and are run once on test. Details, equations and failed attempts are in each
module's docstring and in `candidate1_log.md` / `candidate2_log.md`; numbers are in the JSON files.

| Method | T3 spine, 0.5 s push, 1 mm | T1 ladder r = 30, 1 mm |
|---|---|---|
| B1 best fixed | 0.570 / 0.31 / 0.03 | 0.247 / 0.16 / 0.03 |
| EtR-A (best causal baseline) | 0.407 / 0.27 / 0.10 | 0.077 / 0.13 / 0 |
| B3b truth-k oracle, 0.1 s | 0.120 / 0.15 / 0 | 0.062 / 0.11 / 0 |
| B3b truth-k oracle, 0 s | 0.109 / 0.15 / 0 | 0.022 / 0.07 / 0 |
| **C1** moving-window estimate of stiffness and surface motion + OCP over D, KPF, u_ff (`candidate1_model.py`) | **0.142 / 0.16 / 0.03** | **0.021 / 0.07 / 0** |
| **C2** student MLP distilled from a truth-stiffness teacher (DAgger, expectile loss) + ES on task cost, 3 seeds (`candidate2_learned.py`) | **0.255 ± 0.012 / 0.23 / 0.03** | **0.039 ± 0.001 / 0.09 / 0** |

Paired 95 % CIs of ΔJ:
- C1 vs B1: T3 −0.43 (−0.87, −0.19), T1 −0.23 (−0.55, −0.05). C1 vs EtR-A: T3 −0.26 (−0.40, −0.15), T1 −0.06 (−0.09, −0.04).
- C1 vs B3b 0.1 s: T3 +0.022 (+0.004, +0.048), T1 −0.042 (−0.064, −0.023).
- C2 vs B1: T3 −0.32 (−0.75, −0.07), T1 −0.21 (−0.53, −0.03). C2 vs EtR-A: T3 −0.15 (−0.26, −0.06), T1 −0.04 (−0.07, −0.02).

Ablations (same tuning or training budget):

| Ablation | T3 J | T1 J |
|---|---|---|
| C1 estimator without surface-motion terms | 10.0 (no train-feasible setting) | 4.27 |
| C1 gains only (u_ff = 0) | 0.341 | 0.076 |
| C1 neighbouring/quasi-Newton update instead of ordinary solve | 0.210 | 0.021 |
| C1 +0.1 s latency on the stiffness estimate | 0.129 | 0.068 |
| C1 whole update every 40 ms | 0.275 | 0.166 |
| C2 no teacher (ES only) | 0.329 | 0.069 |
| C2 teacher only (no task-cost stage) | 0.302 | 0.042 |
| C2 gains only | 0.279 | 0.041 |
| C2 no force history | 0.823 | 0.277 |
| C2 decision every 0.096 s | 0.500 | 0.113 |

What the candidates show:
- Prompt causal use of contact information closes most of the T3 gap at the sourced push speed
  (C1 about 95 %, C2 about 50 % of the EtR-A to oracle gap).
- The identified mechanism is separating patient surface motion from stiffness change. Without it the
  estimated bone stiffness is biased by about ×0.45 (the confound that also limits Jin and EtR-A), and bone
  entry is detected after 368 ms instead of 48 ms. Feedforward of the estimated drift gives about half of
  C1's gain over B1; gains alone still beat EtR-A on T3.
- Parts that do not carry the result: the neighbouring update (no gain in quality or time with four decision
  variables), and C2's task-cost ES stage (small gain; most of C2's value comes from the privileged teacher).
- Open negatives: at the fast 0.25 s push both candidates are worse than B1 in J (C1 25.5 with 2 of 40
  scenes diverging after contact loss, although it wins 35 of 40 scenes; C2 7.0; B1 2.27). On soft
  homogeneous tissue C1 is slightly worse than B1 (RMSE 0.04 vs 0.025 N). C1's stiff start prior
  (v2, k0 = 3e4) fixes a stiff homogeneous failure but has not been re-run on test.

## Candidate 1 v2: safety supervisor, checked on a fresh test cohort (seed 4000, never used before)

`candidate1_v2.py` wraps the unchanged estimator and OCP in a supervisor tied to the traced fast-push failure
(push on bone read as drift → over-retraction → contact loss → stale estimator, wrong contact anchor and
integrator windup → impacts): Simplex-style fallback to the cell's B1 gains on contact loss, bumpless removal
of the loss-time integral, estimator freeze and post-contact reset, |vs| ≤ 0.1 m/s, and a stiff start prior
k0 = 3e4. Supervisor tuned on train (49 configurations) with a rule that T3/T1 train J must not exceed v1.

| Cell (fresh, 40 scenes) | C1v2 | C1v1 | B1 | EtR-A | B3b 0.1 s |
|---|---|---|---|---|---|
| T3 0.5 s, 1 mm | 0.153 / 0.17 / 0.03 | 3.25 (1 runaway, start prior) | 0.638 / 0.32 / 0.05 | 0.397 / 0.27 / 0.03 | 0.141 / 0.17 / 0 |
| T1 r30, 1 mm | 0.022 / 0.07 / 0 | 0.022 | 0.205 / 0.15 / 0.03 | 0.071 / 0.12 / 0 | 0.056 / 0.10 / 0 |
| T3 0.25 s, 1 mm | 2.13 / 0.60 / 0.55 | 5.65 (2 runaways) | 2.26 / 0.67 / 0.73 | 6.49 / 0.95 / 0.68 | 1.41 / 0.47 / 0.43 |
| T3 0.25 s, 4 mm | 2.95 / 0.68 / 0.53 | – | 3.41 / 0.77 / 0.78 | 4.78 / 0.88 / 0.78 | 1.44 / 0.48 / 0.43 |
| T3 0.5 s, 4 mm | 0.739 / 0.33 / 0.10 | – | 0.666 / 0.33 / 0.05 | 0.390 / 0.28 / 0.05 | 0.153 / 0.17 / 0 |
| T1 r30, 4 mm | 0.038 / 0.09 / 0 | – | 0.215 / 0.16 / 0.03 | 0.130 / 0.17 / 0 | 0.089 / 0.13 / 0 |

- Sourced cells reproduce on fresh scenes: ΔJ vs B1 T3 −0.48 (−0.98, −0.16), T1 −0.18 (−0.41, −0.04);
  vs EtR-A T3 −0.24 (−0.38, −0.14). C1v2 is within noise of the 0.1 s oracle on T3.
- Fast push: no runaways; fewer scenes above 7 N than B1 (−0.18 (−0.30, −0.08) at 1 mm, −0.25 at 4 mm), but
  ΔJ vs B1 includes 0. A robustness result, not a performance win.
- Negative: T3 0.5 s with 4 mm breathing is not better than B1 and is worse than EtR-A (+0.35 (+0.03, +0.79)),
  with more contact loss. Soft homogeneous accuracy is still slightly worse than B1. The fallback is B1 itself.
- The v1 fresh runaway shows the original v1 test number depended on the start prior; v2 fixes it.
Details in `candidate1_v2_log.md` and `candidate1_v2_results.json`.
