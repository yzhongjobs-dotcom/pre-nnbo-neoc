# Draft: Simulation Experiments (paper section, IEEE T-RO / TMECH style)

Status: working draft written from existing results only. No new experiments were run.
Every number below is copied from the result files named in brackets. Paths are relative
to `experiments/cobot_scan_benchmark/` unless they start with `research/`.

Source key used in the brackets:

| Key | File |
|---|---|
| [R] | `README.md` |
| [G] | `gate_PC125.json` (headroom gate, B0/B1/B2/B3b, seed 3000) |
| [A] | `adaptive_results.json` (EtR-A, EtR-R, Haninger rule, Duan, Song) |
| [M] | `mpc_results.json` (Jin one-step / multi-step / literal, Xue) |
| [L] | `learned_results.json` (E27-style learned VIC) |
| [C1] | `candidate1_results.json`, `candidate1_log.md`, `candidate1_model.py` docstring |
| [V2] | `candidate1_v2_results.json`, `candidate1_v2_log.md` |
| [V3] | `candidate1_v3_results.json`, `candidate1_v3_log.md` |
| [C2] | `candidate2_results.json`, `candidate2_log.md`, `candidate2_learned.py` docstring |
| [P] | `physical_results/results.md`, `summary.csv`, `paired_differences.csv` |
| [S] | `research/simulation_setup_reverse_analysis.md` |

Notation used throughout: C1 is the model/optimisation method (moving-window estimator +
short-horizon optimal choice + supervisor); C1v1, C1v2 and C1v3 are its tested instances; C2 is the
privileged-teacher distilled network. "J" is the pre-declared tuning composite; physical claims are
made in newtons. Paired intervals are 95 % bootstrap intervals of the method-minus-reference
difference over scenes.

---

## V. SIMULATION EXPERIMENTS

The experiments answer five questions, each tied to one experiment family:

- **Q1 (benefit).** Does prompt, causal use of contact information improve normal-force regulation on
  moving, heterogeneous tissue over the best *complete* fixed interaction controller? (Sec. V-C,
  six cells, fresh cohort.)
- **Q2 (information).** Which information carries the benefit? (Sec. V-D1: surface-motion terms,
  estimator bias, detection delay, feedforward versus gains.)
- **Q3 (timeliness).** How fast must the information arrive? (Sec. V-D2: oracle delay, estimate
  latency, update period.)
- **Q4 (computation).** Is solver speed a bottleneck? (Sec. V-D3: neighbouring update; negative.)
- **Q5 (generality of the principle).** Does a learned policy that receives the same information
  through a different route show the same effect? (Secs. V-C4 and V-D4, C2.)

Where the method fails or does not help is reported in the same sections and collected in Sec. V-E.

### A. Experimental Setup

#### 1) Platform model (PC125)

The reference platform is a position-controlled cobot of the UR5 / xArm6 class with a 125 Hz external
interaction loop [R; `plant.py`]. Its constants and their sources are listed in Table I. The robot is
modelled kinematically: the stiff joint servo is not deflected by contact force, as in the
simulation model of Wang (TASE 2024, F09:338-372). The lateral axis follows the scan path under
position control and is not adapted by any method; only the normal axis is admittance and force
controlled. Lateral and shear quantities are therefore reported only, never claimed [R; P].

**Table I. PC125 platform constants and sources.**

| Quantity | Value | Source / status [R; `plant.py`] |
|---|---|---|
| Control period Tc | 8 ms (125 Hz) | UR5/xArm 125 Hz ultrasound systems (F07:353-357, F08:168-172) |
| Transport delay | 1 control period | engineering choice; low end of the sources below |
| Servo lag | first order, τa = 19 ms | UR5e TPI 0.1 per 2 ms (Wang 2024, F09:353-362); carrying it over to 125 Hz is an assumption |
| F/T sensor | 1 kHz, σ = 0.08 N, bias random walk, 30 Hz first-order low-pass | engineering values |
| Velocity cap | 0.25 m/s | engineering value |
| Desired force FD | 4 N | within the 3-12 N phantom range of F07/F08 and above the 1-3 N tissue range of Xue/Moreira/Wang [S, Sec. 4] |

Latency sources considered [R]: the UR `servoj` lag set mainly by `lookahead_time` (0.03-0.2 s; 0.01 s
from SW 5.26); a UR10e forum measurement of about 34 ms cross-correlation delay with `servoj` at 2 ms
and lookahead 0.03 (4-6 ms with `sync()`, UR Client Library issue #86); about 20 ms first-order lag
in Wang 2024; a 125 Hz loop with about 16.6 Hz force bandwidth in Mathiassen 2016 (UR5 ultrasound).
The web sources were reached only through search excerpts, and no servo identification of the xArm
was found [R]. The nominal (1 cycle dead time + 19 ms lag) is at the low, conservative end of these
sources, i.e. it does not inflate the latency that favours adaptation. A latency bracket (PC125dead,
mostly dead time) and a 500 Hz negative control (PC500) are defined in `plant.py`; the declared sweep
(dead time 1-2 cycles, lag 10-60 ms) has no result file and is listed as a gap (Sec. V-G).

**Shared interaction controller.** Every method, including every baseline, acts through the same
normal-axis law [R; `plant.py`]:

  M ÿc + D ẏc + K (yc − yref) = ef + ki ∫ef dt,  y_cmd = yc + KPF ef + u_ff,  ef = FD − F_filt,

discretised by exact zero-order hold with a velocity cap and anti-windup (±20 N s). Bounds:
M 0.25-3 kg, D 20-2500 N s/m, K 0-1000 N/m, ki 0-30 s⁻¹, KPF 0-1e-3 m/N, |u_ff| ≤ 20 mm. A policy may
read only time, the filtered measured force, the commanded and encoder positions and its own state;
the truth columns (contact stiffness under the footprint) are available to oracles only [R; C1; C2].

#### 2) Contact model

The contact is the B0 cylindrical finite footprint (R = 30 mm, L = 40 mm) with a spatially varying
modulus, standard-linear-solid relaxation, Hunt-Crossley damping and Stribeck friction; bone-like
ridges have partial aperture coverage and an optional bump [`plant.py`; R]. At 4 N the contact
stiffness is 1.15 N/mm for 30 kPa tissue and 18.4 N/mm for 1 MPa [R]. Patient motion enters as a
surface trajectory: breathing plus push/retract events. All values are synthetic engineering values,
not a calibrated robot or tissue model [`plant.py`].

#### 3) Tasks and the six evaluation cells

Tasks [R; `scenarios.py`]:

- **T3 spine (headline scanning task).** Soft base 10-40 kPa; spinous-process ridges at 25-35 mm
  pitch, width 8-20 mm, modulus ratio 10-30; scan speed 10-20 mm/s; breathing; 1-2 patient push or
  retract events of amplitude 3-8 mm at random times; 10 s.
- **T1 stiffness ladder (timeliness instrument).** Soft base 10-40 kPa; 1-2 ridges of width 8-25 mm and
  edge 0.5-2 mm; scan speed 10, 20 or 30 mm/s; breathing; 8 s.
- **N1 homogeneous tissue (negative control).** 10 kPa-2 MPa with breathing; 6 s. Used here only as a
  known-behaviour check (3 soft + 3 stiff scenes per configuration).
- **T2 push (non-inferiority / false adaptation)** is defined in the protocol but has no reported
  result; it is listed as a gap.

The evaluation uses six dose cells on PC125 [G; V3]:

| Cell | Task | Fixed dose | Role |
|---|---|---|---|
| T3 0.5/1 | T3 | push ramp 0.5 s, breathing 1 mm | sourced push speed (cough contraction 0.5 ± 0.2 s, Jourdan 2022 [x]) |
| T3 0.5/4 | T3 | ramp 0.5 s, breathing 4 mm | sourced push speed; breathing inside the cited range |
| T3 0.25/1 | T3 | ramp 0.25 s, breathing 1 mm | fast, unsourced push (stress) |
| T3 0.25/4 | T3 | ramp 0.25 s, breathing 4 mm | fast, unsourced push (stress) |
| T1/1 | T1 | ratio 30, breathing 1 mm | ladder |
| T1/4 | T1 | ratio 30, breathing 4 mm | ladder |

Breathing sources [R]: chest surface 4-12 mm at 0.2-0.34 Hz (Shafiq & Veluvolu 2014 [x]); abdominal AP
2.5-13.9 mm (Fahmi 2018 [x]); spine target in prone position 3.4 ± 1.6 mm (Mizonobe 2023 [x]); 0.2 Hz /
5 mm test case of Moreira (F10:875-880). The 4 mm cells therefore match the cited prone-spine value,
while 1 mm is a low dose below the cited ranges, used to show that the gap is not a breathing-amplitude
artefact [R]. Ridge geometry, the 10-30× modulus ratio and flinch or posture-shift speeds are not
sourced; they are treated as dose variables [R].

The cells were chosen by a design round that ran probes on 3-8 in-sample scenes; those probes are
indicative only [R]. They established that (i) every one of three independent proposals had omitted
one controller term, and that omission produced most of their reported headroom (125 Hz homogeneous
tissue 0.49 → 0.16 N; stiffness ladder at r = 30 0.72 → 0.26 N once the term is added); (ii) on fast
platforms (500 Hz cobot, 1 kHz torque plant) a jointly tuned complete fixed controller is near the
oracle (0.035-0.039 N fixed vs 0.028-0.037 N oracle on homogeneous tissue); (iii) only the
combination of the 125 Hz class, a bone-like region and patient motion leaves a gap [R]. The
benchmark is thus deliberately set where the evidence says adaptation can matter, and the same
evidence says where it cannot.

#### 4) Scenario cohorts and splits

Splits are by scene seed [R; V2; V3]:

| Cohort | Seed | Scenes per cell | Use |
|---|---|---|---|
| train | 1000 | 20 (40 for the v2 fast-push supervisor tuning and all v3 work) | tuning and training |
| dev | 2000 | 20-60 | development and checkpoint selection |
| original test | 3000 | 40 | gate, all baselines, C1v1 and C2, run once each; later inspected (see Sec. V-E6) |
| fresh test 1 | 4000 | 40 | C1v2 final test, all six cells; later used read-only for diagnosis of v3 |
| fresh test 2 | 5000 | 40 | C1v3 final test, all six cells, run once after every choice was fixed |

On every cohort all methods are run on the same scenes with the same sensor-noise seeds, so all
comparisons are paired [V2; V3; C2]. In T3 cells the push amplitude (3-8 mm), the number of events
and the ridge ratio are random; in T1 cells the scan speed is random [`scenarios.py`].

#### 5) Metrics

All metrics use the true simulated contact force after a 1 s warm-up [R; `scenarios.metrics`]:

- **Accuracy:** force RMSE, 95th percentile and maximum of |F − FD|, fraction of time inside ±0.5 N.
- **Safety:** peak force; fraction of scenes with peak above FD + 3 N = 7 N ("> 7 N"); contact loss
  (F < 0.5 N) as a time fraction, as the fraction of scenes with any loss, and as the longest loss;
  oscillation RMS above 10 Hz; peak shear (reported only).
- **Tuning composite (declared in advance, used only for tuning and selection):**
  J = mean(((F − FD)/0.5)²) + 20 · loss + mean(max(0, F − 7)²).
- **Computation:** policy call time per update on the host prototype (not real-time certification).

A claim counts only in physical units [R]. J is reported because it is the selection objective, but
every headline difference is also given in force RMSE, and the replayed physical tables [P] add
maximum error, peak force and zero-force time.

#### 6) Statistics

Paired differences use a nonparametric bootstrap over scenes (10 000 resamples, 2.5/97.5 percentiles;
`tune.paired_ci`) [V3]. A difference is called significant only when its interval excludes zero. No
multiplicity correction is applied [P]; with six cells and several references per cell the reader
should treat isolated marginal intervals accordingly. For C2, three training seeds are averaged within
each scene before the bootstrap; the independent sample size remains 40 scenes, not 120 [P; C2].
"Wins" counts scenes in which the method has lower J than the reference [V3].

#### 7) Tuning budgets and the selection rule

All methods are tuned or trained on the train split only, selected under one rule, and run once on
each test cohort. The rule is: **minimum mean train J subject to peak ≤ 7 N in at least 95 % of
train scenes**; if no setting is feasible, the plain minimum J is used and reported as infeasible
[R; `tune.best_fixed`]. In the fast-push cells no method met the peak rule on train: B1 [G], C1v1
T3fast (0 of 49 feasible) [C1], the C1v2 supervisor search [V2] and the per-cell EtR-A tuning [V2].
Comparisons in those cells are therefore among controllers that all violate the train peak rule.

Budgets are summarised with the baselines in Table II. The development budget is not equal across
methods: before its declared 49-configuration tuning, C1 used about 45 exploratory train
configurations [C1]; v2 added about 38 supervisor variants on dev plus a 49-configuration supervisor
search [V2]; v3 used 96 further configurations on train/dev/N1-dev [V3]. B1 has the largest single
search (below). This asymmetry is stated as a limitation (Sec. V-G).

#### 8) How the fixed baseline B1 is set

B1 is the best single fixed controller of the **complete** family (M, D, K, ki, KPF), tuned per cell
[R; `tune.best_fixed`; G]:

1. A grid M ∈ {0.3, 1, 3} kg × D ∈ 12 log-spaced values in [20, 2500] N s/m × K ∈ {0, 300} N/m ×
   ki ∈ {0, 3, 10, 30} s⁻¹ × KPF ∈ {0, 2e-5, 5e-5, 1e-4, 2e-4, 4e-4} m/N, i.e. 1728 settings, each
   evaluated on the 20 train scenes of the cell.
2. Nelder-Mead in log space from each of the 5 best grid points (up to 200 evaluations each).
3. Selection by the rule above.

The resulting gains (M, D, K, ki, KPF) are [G]:

| Cell | M (kg) | D (N s/m) | K (N/m) | ki (s⁻¹) | KPF (m/N) | train-feasible |
|---|---|---|---|---|---|---|
| T3 0.25/1 | 1.589 | 73.5 | 197.5 | 4.56 | 1.45e-4 | no |
| T3 0.25/4 | 1.001 | 152.4 | 677.3 | 30.0 | 1.24e-4 | no |
| T3 0.5/1 | 1.000 | 180.5 | 320.4 | 30.0 | 1.05e-4 | yes |
| T3 0.5/4 | 1.000 | 179.6 | 247.5 | 30.0 | 1.06e-4 | yes |
| T1/1 | 1.688 | 104.3 | 0.0007 | 11.3 | 1.16e-4 | yes |
| T1/4 | 1.681 | 114.9 | 0.0009 | 14.5 | 1.16e-4 | yes |

This is stronger than common practice in the field. In the 36 papers analysed, fixed baselines are set
at the proposed method's own gain bounds (E41: K 100 / 1500 N/m), by per-case trial and error (E04),
by a hand-picked "ideal" damping (Duan RAS 2018), or not reported (E27); none grid-searches the fixed
family and selects the best single setting on a training set [S, Sec. 2]. Two observations from our
data show why this matters. First, removing one term from the fixed family produced most of the
headroom reported by the design proposals (above). Second, the Duan adaptive-damping law maps exactly
onto a fixed PI admittance with K = KPF = 0 [A docstring]; tuned with its own budget it reaches
J 6.20 / 2.59 (T3 / T1) against 0.570 / 0.247 for B1 [A; G], i.e. the missing direct force path, not
adaptation, decides that comparison. B0 extremes (D at 20 and 2500 N s/m, ki 10) are also run and fail
by orders of magnitude (T3 0.5/1: J 1514 and 171) [G].

### B. Baselines

Table II lists every comparison method, how it was mapped onto the shared interface, its tuning budget
and the main fidelity caveat. All share the plant, seeds, score and selection rule.

**Table II. Baselines.**

| Method | Information | Mapping and fidelity (one line) | Budget |
|---|---|---|---|
| B0 fixed extremes | none | D at its bounds, ki 10 | – |
| **B1** best fixed | none | complete family (M, D, K, ki, KPF), Sec. V-A8 | 1728 grid + 5 × ≤200 NM, 20 train scenes [G] |
| B2 per-scene fixed (oracle) | scene identity | best of the same grid on each test scene; upper reference for any fixed controller | per test scene [G] |
| B3b truth-stiffness schedule (oracle) | true footprint stiffness k, delayed 0 / 0.1 / 0.2 s | D = max(Dmin, k/wI), KPF = cP/k, M 0.3, K 0, ki fixed; rule tuned on train per delay. Gains only: an information reference, not an optimal controller or an upper bound for methods with feedforward | rule grid on train [G] |
| EtR-A estimate-then-rule | causal | lock-in stiffness estimate from a force-targeted 20.8 Hz command dither (first-differenced force and position, known F/T filter) + the B3b rule; best causal baseline | 120 configurations [A] |
| EtR-A + Haninger rule | causal | same estimate, D = max(Dmin, 2ζ√(M k̂)) (Haninger RAS 2023, Eq. 20) | 120 [A] |
| EtR-R | causal | RLS on increments; biased by unmeasured surface motion by construction | 45 [A] |
| Duan RAS 2018 | causal | exactly a fixed PI admittance (K = KPF = 0, ki = σ/λ); σ capped at 0.24 by the plant ki bound (paper σ = 0.5 at 5 ms ≡ ki 100). Reproduces the paper's slope-drift correction (mid-ramp bias 0.85 N constant impedance vs 0.0026 N tuned Duan) | 135 grid + refinement; no train-feasible setting [A] |
| Song RA-L 2025 NAIPC | causal | algebraic reduction to admittance + lagged direct force path; mapping error ≤ 0.0007 N; paper constraints KPF·D·ki ≤ 1/4 and M ki < D enforced | 228 grid [A] |
| Jin TMECH 2023 MPC-VIC, one-step / multi-step (anchored) | causal | (log D, KPF) as MPC inputs; RLS stiffness with a train-calibrated gain correction (×1.7 on T1); compliance cost centred at the cell's B1 gains (a documented deviation); LM instead of iLQR; energy tank | 16 / 24 configurations [M] |
| Jin multi-step, literal | causal | compliance cost centred at zero impedance as in the paper; fails (J 16.7 / 4.06) | 8 [M] |
| Xue RAS 2025 MPVIC | causal (+ initial surface position, the paper's Assumption 1) | HC estimator + MPC; explicit-Euler model limits (M, dmin, dmax) to two pairs; ki ≤ 10; HC model undefined once pushes exceed the assumed surface; no train-feasible setting | 24 (+1 fully causal variant) [M] |
| E27-style learned VIC | causal | one-hidden-layer MLP (126 parameters) adjusting log D, log KPF every 0.05 s; trained by separable CMA-ES instead of DRL; learned near-constant gains | 3 seeds × 14 040 episodes [L] |

The mapped field methods are fair implementations on this interface, not reproductions of the original
experiments; why their published results differ from ours is not established by these data [R].
"To be added" baselines that were declared but not implemented are B1-robust (tuned over latency and
filter mismatch) and a breathing feedforward for all methods [R]; both are gaps (Sec. V-G).

**Proposed methods.** C1 (`candidate1_model.py`, `candidate1_v2.py`, `candidate1_v3.py`) estimates,
over a moving window of the last Nw control intervals (16 intervals = 128 ms by default), the local
log-stiffness k, an equivalent surface trajectory h(τ) = vs τ + as τ²/2 and a force level g0, by
arrival-cost-regularised constrained Gauss-Newton on the one-step-inverted F/T filter output with the
1 kHz robot position reconstructed from the known servo model [C1]. A finite-horizon problem
(H = 16 ticks) chooses (log D, KPF, u0, u1) — admittance damping, direct force gain and a ramp-then-hold
reference feedforward u_ff — by one warm-started projected Gauss-Newton iteration every 8 ms, through the
exact robot-side predictor, with equilibrium robustness probes and a rule-form gain prior [C1]. The
supervisor (v2) adds a Simplex-style fallback to the cell's B1 admittance on contact loss, a bumpless
integrator, estimator freeze and post-contact reset, and |vs| ≤ 0.1 m/s, with a stiff start prior
k0 = 3e4 N/m [V2]. v3 changes only configuration selection (each 4 mm cell reuses its 1 mm
configuration) [V3]. C2 (`candidate2_learned.py`) is a causal MLP student around the EtR-A lock-in
front end whose stiffness head is distilled by DAgger from a privileged truth-stiffness teacher with
an expectile loss (τ = 0.95), then fine-tuned by ES on the task cost; one policy serves T1 and T3 [C2].

### C. Main Results

#### 1) Headroom: what knowing the stiffness is worth (seed 3000)

**Table III. Gate cells, 40 test scenes per cell: RMSE (N) / mean peak (N) / fraction of scenes > 7 N;
last column J and paired 95 % CI of the B3b(0.1 s) − B1 difference [R; G].**

| Cell | B1 | B2 per-scene fixed | B3b 0 s | B3b 0.1 s | B3b 0.2 s | J: B1 → B3b 0.1 s |
|---|---|---|---|---|---|---|
| T3 0.25/1 | 0.67 / 7.8 / 0.65 | 0.50 / 7.0 / 0.40 | 0.39 / 6.9 / 0.30 | 0.41 / 6.9 / 0.30 | 0.78 / 8.8 / 0.47 | 2.27 → 1.08 (−1.67, −0.79) |
| T3 0.25/4 | 0.74 / 8.1 / 0.68 | 0.51 / 7.0 / 0.38 | 0.39 / 6.9 / 0.30 | 0.41 / 7.0 / 0.30 | 0.71 / 9.6 / 0.53 | 2.93 → 1.09 (−2.89, −1.13) |
| T3 0.5/1 | 0.31 / 5.5 / 0.03 | 0.19 / 4.9 / 0.00 | 0.15 / 4.8 / 0.00 | 0.15 / 4.8 / 0.00 | 0.30 / 5.3 / 0.03 | 0.57 → 0.12 (−0.91, −0.20) |
| T3 0.5/4 | 0.32 / 5.5 / 0.03 | 0.19 / 4.9 / 0.00 | 0.15 / 4.8 / 0.00 | 0.16 / 4.8 / 0.00 | 0.32 / 5.3 / 0.03 | 0.60 → 0.13 (−0.95, −0.21) |
| T1/1 | 0.16 / 4.8 / 0.03 | 0.11 / 4.6 / 0.00 | 0.07 / 4.3 / 0.00 | 0.11 / 4.4 / 0.00 | 0.13 / 4.7 / 0.05 | 0.25 → 0.06 (−0.50, −0.02) |
| T1/4 | 0.18 / 4.8 / 0.03 | 0.13 / 4.6 / 0.00 | 0.08 / 4.4 / 0.00 | 0.13 / 4.5 / 0.00 | 0.16 / 4.8 / 0.07 | 0.27 → 0.09 (−0.51, −0.01) |

The headroom is real but bounded. Known stiffness at 0.1 s delay lowers J by a factor of 2.7-4.7 in
every cell with every interval excluding zero, and halves RMSE at cough speed (0.31 → 0.15 N) [R].
But the pre-declared strict gate (B1 fails a physical criterion while B3b at 0.1 s meets all) fails
in all six cells: at the sourced 0.5 s ramp and on the ladders B1 already meets every physical
criterion (RMSE < 0.5 N, ≥ 70 % in band, < 20 % of scenes > 7 N), and at the 0.25 s ramp the oracle
itself has 30 % of scenes above 7 N [R; G]. With 0.2 s information delay the same oracle is no better
than B1 (and worse at the 0.25 s ramp), so the usable latency budget is about 0.1 s [R]. On ladders a
per-scene fixed controller nearly matches the 0.1 s oracle, and a stiffness-map feedforward is the
stronger lever (`baselines_ref.TruthKFF`, T1 RMSE 0.055 N with 0.05 s preview vs 0.106 N gains-only
and 0.155 N for B1) [R]. The paper's claim is therefore framed as a force-accuracy and contact-keeping
improvement inside a regime where the best fixed controller is already acceptable at sourced speeds,
not as rescuing a failing fixed controller.

#### 2) Sourced cells, all baselines and candidates (seed 3000)

**Table IV. T3 0.5/1 and T1/1, 40 test scenes: J / RMSE (N) / fraction of scenes > 7 N [R; A; M; L;
C1; C2]. Paired ΔJ CIs vs B1 in brackets where reported in the result files.**

| Method | T3 0.5/1 | T1/1 |
|---|---|---|
| B1 best fixed | 0.570 / 0.31 / 0.03 | 0.247 / 0.16 / 0.03 |
| B3b 0.1 s (oracle) | 0.120 / 0.15 / 0 | 0.062 / 0.11 / 0 |
| B3b 0 s (oracle) | 0.109 / 0.15 / 0 | 0.022 / 0.07 / 0 |
| EtR-A | 0.407 / 0.27 / 0.10 [−0.164 (−0.552, +0.091)] | 0.077 / 0.13 / 0 [−0.169 (−0.491, +0.005)] |
| EtR-A + Haninger | 0.410 / 0.30 / 0.05 [−0.160 (−0.598, +0.094)] | 0.174 / 0.16 / 0.08 [−0.072 (−0.447, +0.200)] |
| EtR-R (RLS) | 1.91 / 0.56 / 0.55 | 0.296 / 0.25 / 0.08 |
| Duan RAS 2018 | 6.20 / 1.06 / 0.68 | 2.59 / 0.62 / 0.38 |
| Song RA-L 2025 | 0.99 / 0.44 / 0.03 | 0.55 / 0.32 / 0.03 |
| Jin one-step (anchored to B1) | 0.549 / 0.32 / 0.03 [−0.021 (−0.147, +0.066)] | 0.200 / 0.15 / 0.03 [−0.046 (−0.131, −0.002)] |
| Jin multi-step (anchored) | 0.604 / 0.35 / 0.05 [+0.034 (−0.279, +0.238)] | 0.123 / 0.14 / 0.03 [−0.124 (−0.447, +0.055)] |
| Jin multi-step, literal | 16.7 / 1.61 / 0.73 | 4.06 / 0.71 / 0.55 |
| Xue MPVIC (initial-surface prior) | 7.76 / 1.28 / 0.85 | 1.36 / 0.55 / 0.20 |
| Learned VIC (E27-style, 3 seeds) | 0.47 ± 0.07 J, 0.32 N | 0.57 ± 0.39 J, 0.18 N |
| **C1v1** | **0.142 / 0.16 / 0.03** [−0.428 (−0.869, −0.185)] | **0.021 / 0.07 / 0** [−0.226 (−0.549, −0.049)] |
| **C2** (3 seeds, mean ± sd) | **0.255 ± 0.012 / 0.23 / 0.03** [−0.315 (−0.746, −0.072)] | **0.039 ± 0.001 / 0.09 / 0** [−0.207 (−0.528, −0.034)] |

Among the mapped field methods none is significantly better than B1 on T3; the only interval that
excludes zero, Jin one-step on T1 (−0.046 (−0.131, −0.002)), is marginal and uncorrected for
multiplicity [A; M]. The learned VIC's seed 0 is significantly *worse* than B1 on T1
(+0.116 (+0.015, +0.260)) [L]. The best causal baseline, EtR-A, recovers most of the oracle value on the
ladder (J 0.077 vs 0.062) but little on the spine task with patient motion (0.407 vs 0.120): its
estimate lags by about 70-100 ms at ridge entry and patient motion corrupts it [R; A]. The remaining
room for a method is therefore causal, prompt use of contact information *while the tissue moves*.
C1v1 recovers 92 % of the EtR-A-to-oracle(0.1 s) J gap on T3 ((0.407 − 0.142)/(0.407 − 0.120)) and C2
53 % ((0.407 − 0.255)/(0.407 − 0.120)) (computed from [C1; C2; A]). C1v1 vs EtR-A: T3
−0.264 (−0.397, −0.154), T1 −0.057 (−0.087, −0.035); C2 vs EtR-A: T3 −0.152 (−0.261, −0.062), T1
−0.038 (−0.065, −0.018) [C1; C2].

Against the oracles the picture is two-sided. On T3, C1v1 is slightly but significantly worse than
B3b at 0.1 s in J (+0.022 (+0.004, +0.048)) [C1]; on T1 it is better than B3b at 0.1 s
(−0.042 (−0.064, −0.023)) and equal to B3b at 0 s (−0.001 (−0.004, +0.002)) [C1]. Because B3b is a
gains-only reference, "better than the 0.1 s oracle" on T1 means that prompt information plus
feedforward beats delayed information without feedforward, not that C1 exceeds an optimum.

**Physical-unit replay (seeds 3000 and 4000) [P].** The replayed tables report RMSE, mean maximum
absolute error and peak force for all frozen controllers on the two sourced cells. Key paired
differences (N):

| Cohort / task | C1v2 − B1, RMSE | C1v2 − B1, max error | C1v2 − B3b 0.1 s, RMSE | C1v2 − B3b 0.1 s, max error | C2 − B1, RMSE | C2 − B1, max error |
|---|---|---|---|---|---|---|
| 3000 / T3 | −0.148 (−0.201, −0.110) | −0.350 (−0.663, −0.113) | +0.009 (+0.002, +0.016) | +0.290 (+0.211, +0.379) | −0.085 (−0.137, −0.047) | −0.178 (−0.486, +0.072) |
| 4000 / T3 | −0.151 (−0.224, −0.097) | −0.344 (−0.657, −0.068) | +0.007 (+0.001, +0.014) | +0.318 (+0.200, +0.442) | −0.093 (−0.161, −0.042) | −0.211 (−0.467, +0.020) |
| 3000 / T1 | −0.088 (−0.150, −0.047) | −0.443 (−0.754, −0.224) | −0.039 (−0.053, −0.027) | −0.172 (−0.249, −0.102) | −0.062 (−0.125, −0.023) | −0.312 (−0.599, −0.106) |
| 4000 / T1 | −0.080 (−0.134, −0.042) | −0.396 (−0.712, −0.178) | −0.032 (−0.043, −0.022) | −0.115 (−0.171, −0.062) | −0.058 (−0.110, −0.021) | −0.283 (−0.578, −0.072) |

So on T3 the oracle's advantage over C1 is small in RMSE (< 0.01 N) but about 0.3 N in the maximum
error and peak force; "C1 matches the oracle" holds for J and RMSE, not for the worst sample of an
episode. C2's reduction of the maximum error on T3 is not significant [P].

#### 3) Six cells, fresh cohort: headline result (C1v3, seed 5000)

C1v3 was fixed on train/dev only and run once on seed 5000; B1, EtR-A (per-cell tuned) and B3b at
0.1 s were recomputed on the same scenes and noise seeds [V3].

**Table V. Seed 5000, 40 scenes per cell: J / RMSE (N) / mean max |F − FD| (N) / fraction > 7 N /
fraction of scenes with contact loss / number of scenes with J > 20 [V3, `final_seed5000`].**

| Cell | C1v3 | B1 | EtR-A | B3b 0.1 s (oracle) |
|---|---|---|---|---|
| T3 0.5/1 | 0.184 / 0.185 / 1.49 / 0.025 / 0 / 0 | 0.379 / 0.290 / 1.51 / 0.025 / 0.025 / 0 | 0.489 / 0.301 / 1.88 / 0.150 / 0.100 / 0 | 0.187 / 0.185 / 1.07 / 0 / 0.025 / 0 |
| T3 0.5/4 | 0.262 / 0.203 / 1.50 / 0.100 / 0.025 / 0 | 0.398 / 0.299 / 1.51 / 0.025 / 0.025 / 0 | 0.450 / 0.295 / 1.71 / 0.075 / 0.050 / 0 | 0.200 / 0.192 / 1.10 / 0.025 / 0 / 0 |
| T3 0.25/1 | 2.468 / 0.610 / 5.46 / 0.600 / 0.400 / 1 | 2.288 / 0.682 / 4.52 / 0.725 / 0.550 / 0 | 7.976 / 0.998 / 9.20 / 0.725 / 0.600 / 4 | 1.995 / 0.529 / 4.22 / 0.400 / 0.375 / 0 |
| T3 0.25/4 | 3.545 / 0.694 / 6.15 / 0.625 / 0.400 / 2 | 2.926 / 0.744 / 4.79 / 0.775 / 0.575 / 0 | 7.487 / 1.059 / 8.33 / 0.750 / 0.675 / 4 | 2.053 / 0.535 / 4.26 / 0.400 / 0.375 / 0 |
| T1/1 | 0.033 / 0.080 / 0.50 / 0 / 0 / 0 | 0.224 / 0.179 / 0.97 / 0.075 / 0.025 / 0 | 0.128 / 0.153 / 0.93 / 0.025 / 0 / 0 | 0.098 / 0.129 / 0.68 / 0 / 0 / 0 |
| T1/4 | 0.047 / 0.098 / 0.61 / 0 / 0 / 0 | 0.274 / 0.199 / 1.02 / 0.125 / 0.025 / 0 | 0.202 / 0.202 / 1.26 / 0 / 0 / 0 | 0.156 / 0.163 / 0.85 / 0.025 / 0 / 0 |

**Table VI. Seed 5000, paired 95 % CIs of C1v3 minus reference; wins vs B1 out of 40 [V3].**

| Cell | ΔJ vs B1 | ΔRMSE vs B1 (N) | Δfrac > 7 N vs B1 | Δloss-scenes vs B1 | ΔJ vs EtR-A | ΔJ vs B3b 0.1 s | wins vs B1 |
|---|---|---|---|---|---|---|---|
| T3 0.5/1 | **−0.195 (−0.250, −0.145)** | **−0.105 (−0.130, −0.083)** | 0 (0, 0) | −0.025 (−0.075, 0) | **−0.305 (−0.455, −0.184)** | −0.003 (−0.035, +0.024) | 38 |
| T3 0.5/4 | −0.136 (−0.264, +0.053) | **−0.095 (−0.129, −0.054)** | +0.075 (0, +0.175) | 0 (−0.075, +0.075) | **−0.188 (−0.248, −0.127)** | +0.062 (−0.041, +0.227) | 38 |
| T3 0.25/1 | +0.180 (−0.635, +1.215) | −0.072 (−0.152, +0.014) | −0.125 (−0.275, +0.025) | **−0.150 (−0.275, −0.050)** | **−5.51 (−10.53, −2.02)** | +0.474 (−0.562, +1.646) | 29 |
| T3 0.25/4 | +0.619 (−0.622, +2.109) | −0.050 (−0.160, +0.068) | −0.150 (−0.325, 0) | **−0.175 (−0.300, −0.075)** | **−3.94 (−7.26, −1.40)** | **+1.49 (+0.10, +3.23)** | 27 |
| T1/1 | **−0.191 (−0.335, −0.083)** | **−0.099 (−0.143, −0.063)** | −0.075 (−0.175, 0) | −0.025 (−0.075, 0) | **−0.096 (−0.147, −0.057)** | **−0.065 (−0.098, −0.036)** | 38 |
| T1/4 | **−0.228 (−0.414, −0.090)** | **−0.101 (−0.148, −0.063)** | **−0.125 (−0.225, −0.025)** | −0.025 (−0.075, 0) | **−0.155 (−0.207, −0.109)** | **−0.109 (−0.167, −0.059)** | 38 |

(Bold: interval excludes zero.)

Reading of the headline cohort:

- **Slow pushes and ladders (four cells).** C1v3 lowers force RMSE against B1 by 0.095-0.105 N in all
  four cells, with every interval excluding zero; J is significantly lower in three of the four. On
  T3 0.5/4 the J difference is not significant (−0.136 (−0.264, +0.053)) and C1v3 has more scenes above
  7 N than B1 (4 vs 1; three of the four at 7.03-7.08 N, one at 10.9 N, a scene that is also 9.95 N for
  EtR-A) [V3]. C1v3 is better than EtR-A in J in all four cells and within noise of the 0.1 s oracle on
  both T3 cells.
- **Fast pushes (two cells).** C1v3 has significantly fewer scenes with contact loss than B1
  (−0.150 and −0.175) but no J gain (intervals include zero), the reduction in scenes above 7 N is not
  significant on this cohort (intervals touch or include zero), and 1-2 of 40 scenes per cell have
  J > 20 where B1 has none [V3]. This is a robustness result, not a performance win (Sec. V-E1).

#### 4) Consistency across cohorts (earlier cohorts, labelled)

The same method family was evaluated on three cohorts as it evolved. The instances differ (v1 → v2 adds
the supervisor and start prior; v2 → v3 changes only the 4 mm configuration choice), and the earlier
cohorts were later inspected, so only seed 5000 is a clean final test for C1v3.

**Table VII. ΔJ vs B1 (paired 95 % CI) across cohorts [C1; V2; V3].**

| Cell | Seed 3000, C1v1 (original test; later inspected) | Seed 4000, C1v2 (fresh for v2; read-only diagnosis for v3) | Seed 5000, C1v3 (fresh, final) |
|---|---|---|---|
| T3 0.5/1 | −0.428 (−0.869, −0.185) | −0.48 (−0.98, −0.16) | −0.195 (−0.250, −0.145) |
| T1/1 | −0.226 (−0.549, −0.049) | −0.18 (−0.41, −0.04) | −0.191 (−0.335, −0.083) |
| T3 0.5/4 | – | **+0.07 (−0.57, +0.69)** (C1v2 negative) | −0.136 (−0.264, +0.053) |
| T1/4 | – | −0.18 (−0.41, −0.05) | −0.228 (−0.414, −0.090) |
| T3 0.25/1 | +23.19 (−0.75, +65.18) (C1v1, 2 diverging scenes) | −0.13 (−0.87, +0.67) | +0.180 (−0.635, +1.215) |
| T3 0.25/4 | – | −0.46 (−1.75, +0.73) | +0.619 (−0.622, +2.109) |

Absolute J values on seed 4000 (C1v2 / B1 / EtR-A / B3b 0.1 s): T3 0.5/1 0.153 / 0.638 / 0.397 / 0.141;
T1/1 0.022 / 0.205 / 0.071 / 0.056; T3 0.5/4 0.739 / 0.666 / 0.390 / 0.153 [V2]. On seed 3000 the re-run
C1v2 gives T3 0.5/1 J 0.139 and T1/1 0.021 [V2, `oldtest`].

The sourced slow-push and ladder results reproduce on three cohorts. The one cell that changed sign is
T3 0.5/4: C1v2 was not better than B1 and significantly worse than EtR-A on seed 4000
(+0.35 (+0.03, +0.79)) [V2]. The cause was traced on train/dev and read-only seed-4000 traces to
configuration selection, not to the breathing model: the 49-configuration search for that cell on 20
train scenes had only 8 peak-feasible settings, the selected one had no excitation dither and a
twice-faster stiffness random walk, and under it stiffness and surface drift were not separable while
the probe followed the breathing surface (median k̂/k 2.12-2.53 on soft tissue across train, dev and
seed 4000, against 1.11-1.14 with the 1 mm configuration) [V3]. The same setting is equally bad on the
1 mm cell (dev J 0.58 vs 0.12), and the constant-velocity-plus-acceleration drift model fits 4 mm,
0.33 Hz breathing to a cubic remainder of about 1e-8 m over the 0.13 s window [V3]. v3 therefore drops
per-cell re-tuning for breathing amplitude. The old negative reproduces for C1v2 on seed 5000
(ΔJ vs B1 +0.624 (+0.007, +1.573)), and C1v3 improves on it there (ΔJ vs C1v2 −0.760 (−1.555, −0.213))
[V3]. This history is reported because the v3 change was motivated by a test-cohort result; seed 5000
is the cohort that controls for it.

#### 5) Candidate 2: a learned instance of the same information principle

**Table VIII. C2 on seed 3000, 3 training seeds (mean ± sd of J; seed-averaged paired CIs) [C2].**

| Cell | C2 J | RMSE (N) | frac > 7 N | ΔJ vs B1 | ΔJ vs EtR-A | ΔJ vs B3b 0.1 s |
|---|---|---|---|---|---|---|
| T3 0.5/1 | 0.255 ± 0.012 | 0.226 | 0.033 | −0.315 (−0.746, −0.072) | −0.152 (−0.261, −0.062) | +0.135 (+0.092, +0.184) |
| T1/1 | 0.039 ± 0.001 | 0.093 | 0 | −0.207 (−0.528, −0.034) | −0.038 (−0.065, −0.018) | −0.023 (−0.041, −0.008) |
| T3 0.25/1 (secondary) | 7.03 ± 2.69 | 0.877 | 0.60 | **+4.76 (+1.07, +10.46)** | −3.02 (−7.60, +1.58) | +5.95 (+2.16, +12.00) |

C2 is one shared policy for both tasks (the baselines are tuned per cell) and is consistent across
training seeds on the sourced cells (sd 0.012 and 0.001). It is better than B1 and EtR-A on both sourced
cells but clearly behind the oracle on T3, and significantly worse than B1 on the fast-push cell (the
EtR-A row there uses the 0.5 s tuning because no fast-cell tuning existed when C2 was tested) [C2]. In
the physical replay C2's RMSE advantage over B1 reproduces on seed 4000 (T3 −0.093 (−0.161, −0.042) N,
T1 −0.058 (−0.110, −0.021) N) [P]. C2 was not run on seed 5000 or on the 4 mm cells. Known-behaviour N1
checks: soft 15 kPa J 0.0025-0.0028 vs B1(T3 gains) 0.0033; stiff 1 MPa J 0.051-0.068 vs B1(T3 gains)
50.2 and B1(T1 gains) 0.051 [C2].

**Computation.** Host-prototype policy-call times, single core, Python, mean of per-episode medians
(p95) [P]: C1v2 0.088-0.090 ms (0.101-0.106 ms), C2 0.015 ms (0.016-0.017 ms), EtR-A 0.025-0.032 ms,
Jin multi-step 0.248-0.276 ms, Xue 0.096-0.101 ms, against an 8 ms control period. Mean per-episode
maxima reach 0.31 ms for C1v2 [P], and the largest single C1v1 call recorded on the seed-3000 test was
8.3 ms [C1, `t_call_max`], i.e. about one control period. These are prototype timings during a
two-process replay, excluding plant execution and communication; they are not real-time certification.

### D. Ablations and Mechanism

All C1 ablations were re-tuned with the final solver under the same 49-configuration budget on train
and run once on the seed-3000 test (C1v1 design, k0 = 5000) [C1]. C2 ablations used the same training
budget, initialisation and seeds [C2].

#### 1) Information: what has to be separated

**Table IX. C1 ablations, seed 3000, J (RMSE N; fraction > 7 N) [C1].**

| Variant | T3 0.5/1 | T1/1 |
|---|---|---|
| C1v1 (full) | 0.142 (0.163; 0.025) | 0.021 (0.067; 0) |
| (a) no surface-motion terms (vs = as = 0) | 10.0 (1.24; 0.925), no train-feasible setting | 4.27 (0.756; 0.70), no train-feasible setting |
| (b) gains only (u_ff = 0; drift still used in prediction) | 0.341 (0.250; 0.075) | 0.076 (0.109; 0.025) |
| (c) neighbouring / quasi-Newton update | 0.210 (0.198; 0.025) | 0.021 (0.067; 0) |
| (d1) +0.1 s latency on the stiffness used by the OCP | 0.129 (0.161; 0) | 0.068 (0.115; 0) |
| (d2) whole update every 40 ms | 0.275 (0.230; 0.05) | 0.166 (0.131; 0.075) |
| B1 / EtR-A / B3b 0.1 s, for reference | 0.570 / 0.407 / 0.120 | 0.247 / 0.077 / 0.062 |

**Surface motion versus stiffness change.** Removing the surface-trajectory terms from the window
model is the only ablation that destroys the method: no setting is train-feasible, and on test J rises
to 10.0 (T3) and 4.27 (T1), worse than B1 by +9.44 (+5.40, +14.55) and +4.02 (+1.50, +7.54) [C1]. On dev
scenes (n = 10) the estimate without these terms is biased low on bone (median k̂/k 0.450 on T3, 0.469
on T1, against 1.026 and 1.049 with them) and on soft tissue (0.335 and 0.427 against 1.084 and
1.229), and the first bone entry is detected later (median 144 ms vs 0 ms on T3; 368 ms vs 48 ms on T1)
[C1, `estimator_bias_dev`]. The same confound — scanning onto bone raises force while the robot
retracts, and surface motion adds force change at a fixed robot position — biases the RLS estimate used
for Jin low by ×0.4-0.6 [M docstring] and corrupts EtR-A during patient motion [R]. This is the
mechanism the paper claims: separating patient surface motion from stiffness change is what makes the
contact information usable.

**Excitation.** The window never identifies k confidently: the posterior-to-prior std ratio of log k
is 0.92-0.99, and the data-only std is 0.34-0.47 with the dither and 0.6-1.0 without it [V3]. Removing
the 20.8 Hz dither raises T3 dev J from 0.135 to 0.205, because bone-entry detection relies on it, and
an excitation-aware dither was harmful in closed loop (T3 dev-12 J 0.31 vs 0.14) [V3]. The dither is
therefore part of the method and part of its cost on soft tissue (Sec. V-E2).

**Feedforward versus gains.** With u_ff = 0 (drift still used for prediction), C1 remains better than
B1 on T3 (−0.230 (−0.608, −0.001)) and than EtR-A (−0.066 (−0.122, −0.018)), but loses about half of
its T3 gain over B1 (gains-only recovers 0.229 of the 0.428 J reduction) and about a quarter on T1
(0.171 of 0.226) (computed from [C1]). On T1 the gains-only variant is not significantly different
from EtR-A (−0.001 (−0.027, +0.042)) [C1]. Both levers are therefore used: the gains carry the
stiffness information, the feedforward carries the drift estimate. This matches the oracle probes
(feedforward is the stronger lever on ladders) [R].

#### 2) Timeliness

- **Oracle delay (Table III).** The truth-stiffness schedule loses its value between 0.1 s and 0.2 s of
  delay in every cell [G].
- **Estimate latency.** Delaying only the stiffness used by the OCP by 0.1 s raises T1 J from 0.021 to
  0.068, i.e. to the level of the 0.1 s oracle (ΔJ vs B3b 0.1 s +0.006 (−0.015, +0.027)); on T3 it does
  not hurt (0.129 vs 0.142; vs B3b 0.1 s +0.009 (−0.0003, +0.019)) [C1]. The T1 advantage of C1 over
  the 0.1 s oracle is therefore a timeliness effect; on T3 at the 0.5 s ramp, timeliness of the
  stiffness alone is not the binding factor (the drift estimate is not delayed in this ablation).
- **Update period.** Running the whole update (estimator and decision) every 40 ms instead of 8 ms
  raises J to 0.275 (T3) and 0.166 (T1); on T1 it is then no longer significantly better than B1
  (−0.081 (−0.462, +0.187)) [C1]. For C2, a decision every 0.096 s raises J to 0.500 ± 0.054 (T3) and
  0.113 ± 0.072 (T1), and two of three seeds never met the train peak rule [C2].

#### 3) Neighbouring update (negative ablation)

The neighbouring / quasi-Newton update was part of the earlier candidate design as a way to update the
decision cheaply when the context changes. It gave no benefit, and the paper claims no solver
acceleration:

- On train it was selected at J 0.210 on T3 against 0.144 for the ordinary one-iteration solve
  (0.025 vs 0.025 on T1), so the ordinary warm-started Gauss-Newton solve was made the final solver
  before any test run [C1, `candidate1_log.md`].
- On test: T3 0.210 vs 0.142; T1 0.0208 vs 0.0205 [C1]. It used fewer rollout sets per update
  (14.2 vs 18.0 on T3) but the call time was the same (mean 126.5 vs 125.4 µs) [C1].
- In a controlled solver study on 6 dev scenes per cell (same hyperparameters, solver swapped), the
  closed-loop J of the three NU variants (T3 0.066-0.068, T1 0.024-0.028) is not better than one
  ordinary iteration (0.0655 / 0.0242), and converging the ordinary solve to 5 iterations changes
  closed-loop J only marginally (0.0649 / 0.0229) while raising the mean update time from 112 to 293 µs
  on T3 [C1, `solver_study`]. With four decision variables and about 0.1 ms per update against an 8 ms
  period, computation is not the bottleneck in this problem; information is.

#### 4) Candidate 2: which parts carry the result

**Table X. C2 ablations, seed 3000, J mean ± sd over 3 seeds, with the seed-averaged paired ΔJ vs
EtR-A (a ΔJ vs the full model is not stored in the results file) [C2].**

| Variant | T3 0.5/1 J | ΔJ vs EtR-A (T3) | T1/1 J | ΔJ vs EtR-A (T1) | train-feasible seeds |
|---|---|---|---|---|---|
| full (teacher + ES) | 0.255 ± 0.012 | −0.152 (−0.261, −0.062) | 0.039 ± 0.001 | −0.038 (−0.065, −0.018) | 3/3 |
| (a) no teacher (ES only, from EtR-A) | 0.329 ± 0.038 | −0.078 (−0.139, −0.027) | 0.069 ± 0.007 | −0.008 (−0.022, +0.004) | 3/3 |
| (b) teacher only (no ES stage) | 0.302 ± 0.054 | −0.104 (−0.171, −0.050) | 0.042 ± 0.001 | −0.035 (−0.061, −0.016) | 3/3 |
| (c) gains only (u_ff residual = 0) | 0.279 ± 0.025 | −0.127 (−0.239, −0.037) | 0.041 ± 0.002 | −0.036 (−0.062, −0.017) | 3/3 |
| (d) no force history | 0.823 ± 0.175 | +0.417 (+0.212, +0.650) | 0.277 ± 0.059 | +0.199 (+0.068, +0.363) | 0/3 |
| (e) decision every 0.096 s | 0.500 ± 0.054 | +0.093 (−0.013, +0.225) | 0.113 ± 0.072 | +0.036 (−0.010, +0.103) | 1/3 |

(sd is the population standard deviation over the three seeds as stored in [C2]; variants d and e
that never met the train peak rule report their last DAgger iterate and are included anyway [C2].)

- The privileged teacher carries most of C2's value: without it (a), T1 is no longer distinguishable from
  EtR-A; the task-cost ES stage adds little on top of the teacher (b vs full) [C2]. A diagnostic
  regression confirmed that the causal force history carries more stiffness information than the
  lock-in estimate (held-out log-RMSE 0.36 vs 0.54) [C2, `candidate2_log.md`].
- Removing the force history (d) makes the policy worse than EtR-A on both cells; a slower decision
  period (e) loses the advantage. Both repeat the C1 conclusions (information and timeliness).
- The feedforward residual contributes little (c vs full); the attempt to distil a stiffness-map
  feedforward failed (the u_ff head formed a high-gain force-to-position loop; dev J 100-400) [C2]. C2
  is therefore mainly a gains-only instance of the principle, whereas C1 also uses drift feedforward.
- Note: c seed 0 is identical to b seed 0, because ES did not improve dev for that seed [C2].

### E. Robustness and Limitations (main text)

#### 1) Fast 0.25 s pushes: robustness only

The 0.25 s ramp is a fast, unsourced case (twice the speed of the sourced cough contraction). There
C1v3 does not improve J over B1 on seed 5000 (+0.180 (−0.635, +1.215) and +0.619 (−0.622, +2.109)),
reduces the fraction of scenes with contact loss (−0.150 (−0.275, −0.050) and −0.175 (−0.300, −0.075)),
does not significantly reduce scenes above 7 N on this cohort, and has 1 (1 mm) and 2 (4 mm) scenes with
J > 20 out of 40 where B1 has none [V3]. On seed 4000 C1v2 had no such scenes and significantly fewer
scenes above 7 N than B1 (−0.175 (−0.30, −0.075), −0.25 (−0.40, −0.125)) [V2], so the >7 N result does
not replicate as significant on the final cohort. C1 stays well behind the 0.1 s oracle at 4 mm
(+1.49 (+0.10, +3.23)) [V3]. C1v1 without the supervisor diverged after contact loss in 2 of 40 test
scenes (J 241.1 and 731.6) [C1]; the v2 supervisor removed those (post-hoc diagnostic on seed 4000:
16.99 J and 4 runaways without the supervisor vs 2.13 with it, T3 0.25/1) [V2].

The residual failure mode is known: a push on bone right after a short contact loss, after which the B1
fallback gains bounce at 0-12 N (15.6 Hz) on 25-29 kN/m bone and the band exit (12 consecutive ticks in
[1, 6] N) cannot be met [V3]. Five stiffness-aware fallback variants were tried on train/dev and
rejected, because a stiffer fallback slows the in-loss approach (about 12 vs 60 mm/s) and each loss
tick costs more than 20 in J, while a lower KPF removes the fast re-contact path [V3]. Every
approach-speed law tried in v2 also failed [V2]. C2 is significantly worse than B1 in this regime
(Table VIII). In short, neither candidate offers a performance gain at fast pushes.

#### 2) Soft homogeneous tissue: slightly worse than B1

On N1 15 kPa scenes, C1v3 has RMSE 0.046 N against 0.026 N for B1 with the T3 configuration and 0.039
vs 0.023 N with the T1 configuration (3 scenes each) [V3]. The decomposition on N1 dev scenes attributes
about a quarter of the gap to the dither (0.046 → 0.040 N without it) and most of the rest to the
soft-side gain regime (D at the 20 N s/m floor, KPF at its 1e-3 cap); less aggressive soft-side settings
recover at most 0.004 N and do not help T3 [V3]. The gap is inside the 0.5 N band by a factor of about
10 and is reported, not fixed. Conversely, on stiff 1 MPa homogeneous tissue the T3-tuned B1 fails
(J 50.2, RMSE 3.34 N) while C1v3 (T3 configuration) has J 0.014 and RMSE 0.060 N [V3];
this shows that B1 is cell-specific, not that C1 is uniformly better (3 scenes, known-behaviour check).

#### 3) No single configuration covers both push regimes

C1v3 uses one configuration per push regime and task (T3 slow, T3 fast, T1). A declared pooled search
for one shared configuration (51 candidates, 40 train scenes per cell) selected the T3 configuration
with the T1 dither level [V3]. On seed 5000 it is as good as or better than the regime configurations
on the four slow/ladder cells (T3 0.5/4 ΔJ vs B1 −0.202 (−0.267, −0.144)) but significantly worse than
B1 on both fast-push cells (+3.86 (+1.06, +7.60) and +1.79 (+0.43, +3.34)), with 4 and 2 scenes above
J = 20 [V3]. The slow-push configurations lose contact in fast pushes, and the fast configuration gives
J 0.35 on T3 and 0.13 on T1 on train [V3]. Deployment therefore needs knowledge of the push regime, or
a regime-switching layer that has not been built.

#### 4) The fallback is B1

The supervisor that removed the fast-push runaways falls back to the cell's B1 admittance on contact
loss [V2]. C1v2/C1v3 is therefore C1 plus the best fixed controller, and B1's tuning is part of its
cost. Its residual failure (above) occurs under these fallback gains.

#### 5) Simulation only, one synthetic plant

All results come from one synthetic 1-D normal-axis plant with a kinematic arm, one latency model at
the low end of the sources, synthetic tissue and unsourced ridge geometry and push speeds [R;
`plant.py`]. The design-round probes indicate that on faster platforms (500 Hz, 1 kHz torque) a jointly
tuned fixed controller is near-optimal [R], so the benefit shown here is specific to the 125 Hz
position-controlled class. The replay of archived metrics is not exact across environments: 19 of 232
checks differ beyond a strict tolerance, including unstable Jin-literal trajectories and C1v1 on seed
4000 [P]; the C1v1 fast-push seed-3000 value was 25.5 originally and 17.2 on re-run, because runaway
scenes are chaotic [V2].

#### 6) Test-contamination history and how fresh cohorts address it

- Seed 3000 was the original test cohort for the gate, all baselines, C1v1 and C2. After the C1v1 test,
  two diverging test traces (scenes 6 and 36) were used to diagnose the fast-push failure (no tuning on
  them) [V2], and the start-prior issue was found on N1 [C1]. Seed 3000 numbers are therefore
  "already known" for every later instance [V2].
- Seed 4000 was fresh for C1v2 and all recomputed baselines [V2]. Its traces were then re-simulated
  read-only to confirm the T3 0.5/4 mechanism for v3 [V3], and the physical replay extended it to all
  controllers [P]; it is not a clean test for C1v3.
- Seed 5000 was used once, after every v3 choice was fixed, for C1v3, C1v3-shared, C1v2, B1, EtR-A and
  B3b [V3]. It is the headline cohort.
- Residual concerns that fresh cohorts do not remove: the decision to make v3 was triggered by a
  test-cohort negative; the cells and doses were chosen after in-sample design probes [R]; C2 has only
  the seed-3000 test (plus the seed-4000 physical replay) [C2; P].

### F. Figures Needed

| Fig. | Content | Data source |
|---|---|---|
| 1 | Block diagram: window estimator (k, vs, as, g0) → short-horizon choice of (D, KPF, u0, u1) through the robot-side predictor → shared admittance; supervisor and B1 fallback; C2 as teacher/student with the same inputs | method description, no data |
| 2 | Example T3 scene: surface trajectory (breathing + push), stiffness under the footprint, force for B1 / EtR-A / C1v3 / B3b | requires a replay of one stored seed-5000 scene with `plant.episode` (no new tuning) |
| 3 | Headroom: J and RMSE of B1, B2, B3b at 0 / 0.1 / 0.2 s per cell; shows the bounded headroom and the 0.1 s budget | `gate_PC125.json` |
| 4 | Sourced cells, all baselines: RMSE, mean max error and peak force with scene-bootstrap intervals, seeds 3000 and 4000 | `physical_results/summary.csv`, `paired_differences.csv` |
| 5 | Headline forest plot: paired ΔJ and ΔRMSE of C1v3 vs B1, EtR-A, B3b 0.1 s for the six cells (seed 5000) | `candidate1_v3_results.json` → `final_seed5000` |
| 6 | Cross-cohort consistency: ΔJ vs B1 per cell for seeds 3000 / 4000 / 5000 with instance labels | `candidate1_results.json` (`test`), `candidate1_v2_results.json` (`fresh`, `oldtest`), `candidate1_v3_results.json` |
| 7 | Mechanism: k̂/k (soft, bone) and bone-entry detection delay with and without surface-motion terms; k̂/k under the v2 vs 1 mm configuration at 4 mm breathing | `candidate1_results.json` → `estimator_bias_dev`; `candidate1_v3_log.md` table (values in `candidate1_v3_results.json` → `development_notes`) |
| 8 | Ablation bars for C1 and C2 (J and RMSE, both sourced cells) with B1 / EtR-A / oracle lines | `candidate1_results.json` → `test`; `candidate2_results.json` → `cells` |
| 9 | Timeliness: oracle delay curve (0 / 0.1 / 0.2 s) overlaid with C1 latency (+0.1 s) and update-period (40 ms) ablations and C2 decision-period ablation | `gate_PC125.json`, `candidate1_results.json`, `candidate2_results.json` |
| 10 | Solver study: closed-loop J and update time vs rollout sets for NU and ordinary solves (negative result) | `candidate1_results.json` → `solver_study` |
| 11 | Fast-push robustness: per-scene J of C1v3 vs B1 (log axes) and contact-loss scene fractions, both fast cells; marks the 1-2 runaway scenes | `candidate1_v3_results.json` → `J_per_scene`, `loss_per_scene` |
| 12 | Residual failure trace: loss on bone, B1-fallback bounce at 15.6 Hz | requires a replay of T3_4 train scene 36 or fast-push dev scenes 6/43 (named in `candidate1_v3_log.md`) |
| 13 | Computation: policy-call median and p95 per method vs the 8 ms period | `physical_results/summary.csv` |

### G. Gaps a Reviewer Will Raise

1. **No hardware validation.** All evidence is from one synthetic plant; no phantom or volunteer
   experiment, and no image-quality or task-level outcome (only normal force).
2. **Sensitivity to platform latency.** The declared sweep (dead time 1-2 cycles, lag 10-60 ms) and the
   PC125dead bracket have no results; the latency structure "must be identified on hardware before any
   claim" [R]. The benefit depends on the 125 Hz class.
3. **Other robots.** Torque-controlled and 500 Hz-1 kHz platforms appear only in in-sample design
   probes, which suggest the fixed controller is near-optimal there [R]. A negative-control run on
   PC500 with the final methods is missing.
4. **Missing declared baselines.** B1-robust and breathing feedforward for all methods were declared but
   not implemented [R]. Jin was anchored to B1 and Xue's grid was limited by its numerical model; a
   reviewer may ask for closer reproductions.
5. **Missing task.** T2 (non-inferiority / false adaptation on homogeneous tissue with pushes) was
   defined but not run; N1 has only 3 + 3 scenes per configuration.
6. **Fast pushes.** No method, candidate or baseline, meets the peak rule; C1 gives robustness (fewer
   loss scenes) but no J gain and occasional runaways; the fallback-bounce failure is unresolved.
7. **Regime dependence.** No single configuration covers both push regimes; regime detection is not
   implemented.
8. **Unequal development budget.** C1 used more exploratory configurations than the mapped baselines
   (Sec. V-A7), though B1 had the largest single search. A budget-matched comparison or a sensitivity
   to C1's hyperparameters is not reported.
9. **Statistics.** 40 scenes per cohort, many intervals, no multiplicity correction; several
   significant differences are marginal.
10. **C2 coverage.** C2 has not been evaluated on seed 5000 or on the 4 mm cells.
11. **Computation claims.** Timings are host Python prototypes; one recorded single call reached
    8.3 ms; no real-time implementation exists.
12. **Unsourced doses.** Ridge geometry, modulus ratio and push speeds are engineering choices; the
    sourced push speed is where the gain is bounded (both controllers inside the 0.5 N band).

---

### Author notes: inconsistencies found between files (resolve before submission)

1. `README.md` states that without surface-motion terms "the estimated bone stiffness is biased by
   about ×0.45 ... and bone entry is detected after 368 ms instead of 48 ms". The ×0.45 is the T3 value
   (T1: 0.469), while 368 vs 48 ms is the T1 value (T3: 144 vs 0 ms); all are dev medians over 10 scenes
   [C1, `estimator_bias_dev`]. This draft reports both cells.
2. `README.md` says C1 closes "about 95 %" of the EtR-A-to-oracle gap on T3. On seed 3000 the computed
   value is 92 % ((0.407 − 0.142)/(0.407 − 0.120)); about 95 % holds on seed 4000 for C1v2
   ((0.397 − 0.153)/(0.397 − 0.141)). C2 is 53 % (README "about 50 %").
3. `README.md` (gate) says all six B1 solutions have "ki up to its bound". In `gate_PC125.json` ki is 30
   only in three cells; it is 4.56 (T3 0.25/1), 11.3 (T1/1) and 14.5 (T1/4). K also ranges 0-677 N/m.
   The README also does not mention that B1 is not train-feasible in the two fast-push cells (train J
   includes the 1000 infeasibility penalty).
4. `README.md` (v3 summary) lists "fewer contact-loss scenes than B1 (−0.15 and −0.18)"; the v3 log and
   JSON give −0.150 and −0.175. The seed-4000 claim of fewer scenes above 7 N (CIs excluding 0) does not
   replicate as significant on seed 5000 (CIs touch or include 0).
5. "Within noise of the 0.1 s oracle on T3" (README, v2 and v3) holds for J and RMSE, but the physical
   replay shows the oracle significantly better in RMSE by +0.007 to +0.009 N and in mean maximum error
   by about +0.29 to +0.32 N on seeds 3000 and 4000 [P].
6. The README calls T3 0.5/1 and T1/1 the "sourced" cells. The 0.5 s ramp is sourced, but 1 mm breathing
   is below every cited breathing range, while 4 mm matches the prone-spine value. The sourced label
   should refer to push speed only.
7. The EtR-A row for the fast-push cell differs by file: `candidate2_results.json` uses the T3 0.5 s
   tuning (seed 3000 J 10.05), whereas `candidate1_v2_results.json` re-tuned EtR-A for that cell
   (seed 3000 J 6.03). Label the tuning wherever EtR-A appears in fast-push tables.
8. The C1v1 fast-push value on seed 3000 is 25.5 in `candidate1_results.json` but 17.2 on the re-run in
   `candidate1_v2_results.json` (`oldtest`); the physical replay also flags 4000/T3/Candidate1_v1 as
   differing beyond tolerance [P]. Runaway scenes are not reproducible to tolerance.
9. `baselines_learned.py` says "~130 parameters"; `learned_results.json` records 126. Minor.
10. `physical_results/results.md` says Xue "has kmax=0"; `mpc_results.json` specs do not list kmax, so
    this could not be checked from the results file.
