# Redesigned scan task: representative tasks and headroom check

Purpose: replace the synthetic 2DOF task of
[dual_candidate_task_simulation](../dual_candidate_task_simulation/README.md),
whose optimum under the shared physical score was a gain-box corner (fixed
(Kx,Kz,Dx,Dz)=(1500,30,100,4) beat Jin and both candidates in 24/24 test scenes) with a plant whose difficulty is comparable to
representative field work (E41, E27, Xue RAS2025, Haninger RAS2023, Hammoud),
and first check whether environment-dependent impedance has any value there
(step 0 of the redesign plan). No candidate method is run yet.

## Plant (`plant.py`)
- Planar 3R arm (x, z, pitch), joint-torque interface, 5 kHz integration, 1 kHz
  Cartesian impedance loop, 2 ms transport delay, 100 Hz torque bandwidth.
- Robot-side model error: 1% mass error, unknown payload 0-0.1 kg, uncompensated
  joint Coulomb/viscous friction (scale 0.5-1.5).
- B0 cylindrical finite footprint (R 30 mm, L 40 mm, B0 psi relation) with
  spatially varying modulus E(x) and friction mu(x); standard-linear-solid
  relaxation, Hunt-Crossley damping, Stribeck friction, horizontal pressure
  component on the curved face (-m1/R).
- Respiratory surface motion 2-6 mm at 0.2-0.33 Hz.
- Wrist F/T noise (0.05/0.08 N, 3 mNm) with bias random walk; encoder and
  velocity noise. Target force 4 N.
- Adapted gains Kx, Kz, Dx, Dz with bounds and rate limits; pitch impedance fixed.
- `mode=1` adds an explicit force PI outer loop on Fz (kp 1.5, ki 10 1/s),
  a standard strong simple method.

## Tasks (`scenarios.py`)
- T1 boundary: soft<->stiff boundary (ratio 2-6, width 1-6 mm) with friction change.
- T2 inclusion: stiff inclusion (ratio 3-8, width 8-25 mm).
- T3 homogeneous: E 20-300 kPa, respiration and robot-side error only.
Score: tolerance-normalised force error (0.5 N) and lateral error (1 mm), plus
contact-loss, over-force (F_d+3 N) and shear (2.5 N) penalties.

## Results (6 scenes per task, seed 2026)
`headroom.py`: 72-point fixed-gain grid. `switched_oracle.py`: oracle schedule
that switches (Kz, Dz) at the known boundary, with 0/0.3 s lead.

| Task | Force PI | Best single fixed | Per-scene oracle gain | Boundary-switched oracle gain |
|---|---|---:|---:|---:|
| T1 | no | 3.090 | 0.2% | 9.0% |
| T1 | yes | 0.861 | 0.3% | 0.3% |
| T2 | no | 2.995 | 0.9% | 15.2% |
| T2 | yes | 0.772 | 0.2% | 0.3% |
| T3 | no | 2.460 | 9.5% | n/a |
| T3 | yes | 0.558 | 0.3% | n/a |

With force PI the force RMSE is about 0.04-0.06 N; the remaining score is the
lateral error (about 0.9 mm), with Kx and Dx at their upper bounds in every scene.

## What this shows
- Without force feedback, knowing the environment is worth 9-15% (time-switched
  Kz/Dz at material changes), consistent with the 1D pilot.
- With a standard force PI loop, the same knowledge is worth about 0.3%:
  force feedback absorbs boundary transients, respiration and robot-side bias.
- Kx/Dx again sit at their bounds: nothing in these tasks penalises tangential
  stiffness; lateral error would also be removed by integral action.

Interpretation: these three tasks are benign (slow changes relative to a force
loop, clean 1 kHz F/T, 2 ms delay, at most about 8 N/mm contact stiffness, no
cost on tangential stiffness). In that regime fixed impedance plus a force
loop is known to suffice (cf. E04/E27 planar cases), so the result shows the
tasks still lack the conditions under which variable impedance matters; it is
not evidence against the research question. Kx/Dx at their bounds is the same
missing-tradeoff symptom as the old task. Two lessons carry forward: any
comparison must include fixed impedance plus a force loop, and force-loop
gains/feedforward are part of the interaction controller a method may adapt.
Limits: six scenes per task, coarse grids, synthetic parameters, two-period
delay only, piecewise oracle rather than full time-varying OCP.

Next representative tasks to add: sharp stiffness edges (bone/hard inclusion,
10-30x) at higher speed with realistic force-sensing bandwidth and latency
(e.g. a 100 Hz ROS loop as in Xue RAS2025); sudden tissue push/retraction
(Xue); a tangential tolerance so that high tangential stiffness has a cost.

## Reproduce
`python headroom.py 6` (~2.5 min on 4 cores), `python switched_oracle.py` (~7 min).
Requires numpy, scipy, numba.
