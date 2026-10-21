# Cobot scanning benchmark: task protocol and headroom gate

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

## Results

(pending)

## Reproduce

`python gate.py PC125` (all cells) or `python gate.py PC125 0,1` (selected
cells). Requires numpy, scipy and numba.
