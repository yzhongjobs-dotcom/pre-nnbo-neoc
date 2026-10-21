# Reverse analysis of simulation and experiment setups (2026-10-21)

Purpose: work out, from how top papers set up their experiments, how this
project should set up its simulation tasks and baselines, so that the task fits
the research question, every baseline is reasonable and documented, and the
evidence can carry either an application story (scanning) or a method story
(prompt near-optimal impedance update from new contact information).

Complements `independent_reverse_analysis_of_top_journal_methods.md`, which
covers how methods are formed, not how experiments are built.

Scope: experiment/simulation sections of 36 local papers (text extractions in
`literature/extracted_text/`). Some table bodies were lost in extraction
(E41 Fig. 4 values, F06 Tab. I, NNBO Tabs. II-VI, E27/E04 Tabs. II-III,
F08 Tab. II, R13 Tab. I); only numbers stated in the body text are used.
Line references below are to those text files and are approximate.
Corrections to earlier labels: E16 is Raiola et al. RA-L2018 (safety/energy
aware impedance), E47 is Cao, Cheng, Li TCDS2024 (passive MPC impedance for pHRI).

## 1. How the experiments derive from the claim

Strong papers name each claimed capability and build one scenario in which a
method lacking that capability visibly fails. Weak papers run generic demos.

| Claim type | Scenario built to need it | Examples |
|---|---|---|
| Precision-compliance trade-off needs online K | Reference runs into an obstacle: too soft stalls (accuracy fails), too stiff gives large force (compliance fails) | E41 (Franka, curved obstacle, fixed K 100 vs 1500 N/m = the method's own bounds); E38 grinding: Qf=0 gives ~40 N and damage, Qf=0.5 oscillates, tuned Qf=0.2 works |
| Stiff contact needs damping matched to stiffness | Soft / uncertain-location / stiff environment grid | Haninger RAS2023 (15, 45, 126 N/mm, ±1 cm; constraint D >= 2 zeta sqrt(M Ke)); low-impedance MPC loses contact on stiff |
| Uncertainty makes high gains costly | Unexpected contact earlier than planned | Hammoud (3 cm block, DDP vs risk-sensitive with identical cost, 500 Monte Carlo runs, trot success 47% -> 65%) |
| Force tracking on unknown environment | Ladder: stiffness step -> slope -> sine -> composite, free space to contact | Duan RAS2018 (1000->3000->5000 N/m, fd 5-50 N), E04 (3000->5000 N/m at 5 s, fd 50 N), E18, E27 |
| Contact keeping / safety on moving tissue | Tissue retracts (contact loss) and pushes (force spike) at random times | Xue RAS2025 (silicone 5-18 kPa, 100 Hz ROS loop; RMSE with/without disturbance) |
| Physiological motion rejection | 0.2 Hz / 5 mm breathing, recorded heartbeat trace | Moreira CMPB2014 (peak error 0.9 N -> 0.1 N, compensation 87%) |
| Each adapted quantity has its own role | Phases that each need a different component (feedforward, then K, then D, then disturbance removed) | R11 (constant field, position-dependent and velocity-dependent divergent fields; ILC-only ablation fails) |
| Robustness to the competitor's weak point | Scenario forcing the competitor's assumption to fail | R17 (saturation tightened so violation must occur), R19 (rotation/noise/bias sweeps) |
| Speed / latency | Difficulty-variable sweep plus latency breakdown | YOPO (speed 2-10 m/s, 10 repeats; latency 1.6 ms vs 5-62 ms) |

Typical paper: 2-5 scenarios, each tied to one claim, plus 1-3 ablations
("same method, component off") and often one sweep over the variable that
creates the difficulty (stiffness, speed, noise, friction).

## 2. How fixed and other baselines are set

| Practice | Papers | Assessment |
|---|---|---|
| Fixed K at the proposed method's own bounds (high = precise, low = compliant) | E41 (100 / 1500 N/m), E38 (index extremes) | Shows both sides fail; standard for trade-off claims |
| Fixed parameters "optimal by trial and error per method and per case", values listed | E04 (md 2, bd 200, kd 0 in sim; md 26, bd 2500 real) | Fair and documented, but per-case tuning |
| Fixed damping hand-picked as "ideal", then shown failing elsewhere | Duan RAS2018 (b=42 plane, b=60 slope, no b works on sine) | One-sided |
| Same impedance vector for all methods | Xue, Wang TASE2024, E47 | Isolates the method's added mechanism, not the value of adaptation |
| Iso-accuracy: tune all to the same tracking accuracy, compare force | F06 (Jin TMECH2025) | Good fairness protocol for trade-off claims |
| Sweep of fixed set-points, show the adaptive method gets high-setting quality at low-setting cost | F07 (breast US: matches 8/6/3 N constant force quality with less force), R21 (4 fixed thrusts; adaptive beats best fixed by about 20 points sim, 11 real) | Strongest, reviewer-proof form |
| Values not reported | E27 constant impedance, E41 MPC-VIC/RL-VIC tuning | Weakness to avoid |
| Same cost, weights and reference; only the mechanism differs | Hammoud, F04, NNBO, F01, R13 (same optimiser) | Correct for method ablations |

No paper in the set grid-searches fixed K/D and picks the best single fixed
setting on a training set and per scenario. Doing so (plus the two extremes)
is the most defensible fixed baseline and is a gap this project can fill.

Where fixed impedance is shown to fail in the field:
- slope and sinusoidal surfaces, as a steady-state force bias (E04, Duan, E18, E27);
- obstacles and unexpected contact, as a force spike (stiff side) or a stall or tracking loss (soft side) (E41, Hammoud);
- stiff environments with low damping, through oscillation or contact loss (Haninger, E47 with zeta 0.01);
- moving tissue, through contact loss on retraction (Xue).

On a flat plane with a stiffness step, fixed impedance with force feedback
generally copes (E27: "all three methods coped with the stiffness switches";
E04 plane: all similar). Our earlier headroom result (fixed + force PI within
0.3 % of the oracle on 2-8x stiffness steps) is therefore consistent with the
field, not a contradiction of it.

## 3. Metrics

- Accuracy: force RMSE / MAE / max error [N] or % of fd; overshoot [N];
  settling time [s]; position/tracking RMSE and max [mm].
- Compliance / safety: max and RMS contact force [N]; impact or push-induced
  force rise; contact-loss indicator; oscillation (energy above ~15 Hz, Haninger);
  passivity / dissipated energy [J], power [W] (E16, E47); stiffness level.
- Cost of the method: computation time per update vs control period (Xue, F04,
  Haninger, YOPO); training samples / time (E27, E40, R19).
- Two-sided reporting (accuracy and compliance in separate physical units) is
  standard for trade-off claims (E41, F06, E38, F07, R05).

## 4. Magnitudes used in the field

| Quantity | Range in the papers |
|---|---|
| Robot task-space K | 100-1500 N/m (E41), 500-5000 N/m (E37, E40), K = 0 with damping-only adaptation (E04, E27, Duan) |
| Damping | 20-500 Ns/m; damping ratio 0.1-0.9 |
| Environment | soft tissue 1-18 kPa (Xue); linear 40-5000 N/m; soft/hard 500 vs 20000 N/m (E44, 40x); 15-126 N/mm (Haninger); OOD 250-25000 N/m (E27) |
| Desired force | 1-3 N tissue (Xue, Moreira, Wang), 3-12 N phantoms (F07, F08; >5 N no image gain), 5-50 N industrial |
| Scan speed | 1-15 mm/s typical (F07, E04), up to 0.15 m/s (Wang) |
| Surface motion | breathing 0.2 Hz / 5 mm (Moreira), random retract/push (Xue) |
| Loops | force/impedance 1-5 ms; gain update via MPC/learning 5-100 Hz; actuation lag ~20 ms modelled (Wang TPI 0.1) |
| Sensor noise | rarely quantified; E27 uses 30 % of the tolerance band |

## 5. Statistics and the split between simulation and hardware

- Statistics are weak across the field: mostly single runs; best are Hammoud
  (500 Monte Carlo), R19 (10 seeds, 95 % CI), E40 (20 trials), Haninger (3).
- Simulation carries the controlled evidence: ground truth, sweeps, forced
  assumption failure, Monte Carlo. Hardware carries plausibility and transfer
  with 3-20 trials; the best papers repeat the same group design on hardware (R19, R21).
- Learning papers split train/test by environment family and add an
  out-of-distribution sweep (E27 random polynomial terrains for training,
  K_env 250-25000 N/m x fd 10-50 N for testing).

## 6. How the story type changes the design

| | Application story (Xue, AuSoScan F08, Wang F09, Zhang F07, R01, R05) | Method story (E41, F06, Haninger, Hammoud, NNBO, E04, E27, R11, R17, R19, YOPO) |
|---|---|---|
| Main evidence | Real platform, phantoms or volunteers, task-quality metric (image confidence, reconstruction, Cobb angle) | Benchmark-like scenarios each isolating one capability, plus ablations and sweeps |
| Baselines | Same controller without the new model or mode (F08), same impedance vector (Xue), fixed set-point sweep (F07) | Fixed extremes, tuned representative methods, same-cost ablations, oracle optimum (E45, F03) |
| Theory | Light | Conditions checked by diagnostics (R03 excitation, F01 certificate eigenvalue) |
| Risk | Single-domain significance; few trials | Benchmark looks artificial unless tied to a real task |

A method story with scanning as the representative application can combine
both: field-standard benchmark scenarios (which reviewers already accept as
the places VIC matters) plus scanning scenarios that motivate the problem.

## 7. Implications for this project

1. Each scenario must contain a mechanism under which the best fixed
   interaction controller is wrong somewhere: soft contact wants low damping and
   high force-loop bandwidth; stiff (bone-like) contact wants high damping or
   lower force gain; free-space or obstacle segments want high K for precision
   and low K at contact. Stiffness steps of 2-8x on a flat soft plane with a force
   loop do not create this (field and our headroom check agree).
2. Fixed baselines: two extremes at our own gain bounds, plus the best single
   fixed setting grid-searched on the training cohort, with K/D and force-loop
   gains tuned jointly. Report the per-scenario oracle fixed setting as an upper
   reference. All values are listed.
3. Representative adaptive baselines are re-tuned by the same budget on the same
   training cohort and score: Duan/E04 adaptive force tracking; Jin E41 MPC-VIC;
   Xue MPVIC (online HC estimate + MPC); a learned VIC in the E27 style;
   force PI / admittance (Wang type) as the strong simple method.
4. Metrics are two-sided in physical units: force RMSE, max and overshoot
   [N], lateral tracking [mm], contact loss, peak force, oscillation energy; plus
   computation time per update and training cost.
5. Monte Carlo cohorts (>= 20 seeded scenes per family), paired comparisons
   with bootstrap CIs, train/dev/test split by seed, an OOD family.
6. Research questions are stated up front and mapped one to one to experiment
   families (R19 style); ties and losses are reported (R13, YOPO style).

First physics probe in the existing plant (`experiments/realistic_scan_task`,
fixed Kz=10 N/m, force PI kp 1.5, ki 10, 4 mm breathing, 5 s): with the force
loop, the best Dz on 30 kPa tissue is 5 Ns/m (force RMSE 0.029 N; 0.117 N at
Dz 150), whereas on a 1 MPa bone-like region Dz <= 20 gives oscillation
(RMSE 1.7-2.1 N, peaks 10-19 N) and Dz 150 is needed (0.119 N). So the damping
needed depends strongly on the environment once a bone-like region is present,
which matches Haninger's well-damped condition. However, when the force-loop
gains are tuned jointly with Dz (kp 0.3-3, ki 5-40, Dz 5-150), one fixed setting
(kp 0.3, ki 20, Dz 60) gives 0.035 / 0.034 / 0.038 N on 30 kPa / 150 kPa / 1 MPa
against per-environment optima of 0.025 / 0.030 / 0.038 N: about 15 % headroom,
with all errors far inside a 0.5 N band. So with a clean 1 kHz force loop and
2 ms delay, a jointly tuned fixed controller is close to optimal across a 30x
homogeneous stiffness range. Any task set must be checked against this
loophole. The field shows variable impedance mattering under conditions this
probe lacks: transitions and unexpected contact, actuation and sensing lag of
position-controlled platforms (about 20 ms, Wang TASE2024; 100-125 Hz loops,
Xue, Zhang, AuSoScan), and tissue push/retract. Those conditions must come from
a stated reference platform, not be inflated to favour adaptation.
