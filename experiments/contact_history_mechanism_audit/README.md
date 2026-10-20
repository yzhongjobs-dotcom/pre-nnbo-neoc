# Does contact history change future response and impedance selection?

2026-10-18. **Exploratory mechanism check** after literature reading, compared with `57e08ca`. This is not a new scanning controller, strong-prior reproduction or unknown-parameter identification experiment, and does not rank NEOC/PPO.

## Motivation and exact questions

Tiezzi/Kao T-RO2007 studies finite-area relaxation, creep and loads; Moreira et al. CMPB2014 compares tissue models and compensates physiological motion. Contact memory is established knowledge, not our discovery. NNBO/Xue HC parameterizations, AuSoScan sliding contact and Matschek static GP force maps each have capabilities/conditions; none should be dismissed as lacking contact models.

Ask only: **Can identical current indentation, velocity and force produce different future forces, and does that difference materially improve impedance selection?** The first does not prove the second.

Use a known passive generalized Maxwell model with two relaxation branches:

\[
 F=k_\infty\delta+z_1+z_2,\qquad
 \dot z_i=k_i\dot\delta-z_i/\tau_i.
\]

Stored energy is \(S=\frac12k_\infty\delta^2+\sum_i z_i^2/(2k_i)\), satisfying
\(\dot S=F\dot\delta-\sum_i z_i^2/(k_i\tau_i)\). Under sustained compression/fixed parameters this is an established passive rheological element, not new contact theory or full-robot passivity. Parameters are explicit mechanism examples, **not calibrated to real tissue**.

Two realizable monotone smooth loading histories last 20ms, ending now or 150ms earlier. Adjust loading amplitude rather than arbitrary hidden states so both now have 4mm,0m/s,3N, then hold indentation. Material parameters agree; internal states differ. Fixed-parameter HC \(F=(k+c\dot\delta)\delta^n\) without extra states is single-valued in current motion and cannot exactly represent both holding futures. **This does not deny good fits to some measured trajectories or tracking relaxation through online-varying parameters.** Moreira reports good HC fits; its experiments cannot be rewritten as proving HC failure.

## Results and limits

- Holding-force curves differ by at most about **0.2806N** at 128ms, and about 0.2752N at 100ms. Any single-valued holding prediction based only on the identical snapshot has error at least 0.1403N on one history at 128ms. This is model information aliasing, not a general tissue-error lower bound.
- With true branch parameters/full known displacement history, a standard two-state linear observer recovers/predicts internal states. At 100Hz force sampling and independent Gaussian noise 0.01/0.05N,100 repetitions each give mean holding RMSE about 0.0008/0.0040N. Parameters, displacement and preparation are ideal: **not unknown-tissue identification or out-of-domain generalization**. Zero-noise reconstruction error is about 4.4e−16N.
- With the same 1kg normal plant, initial actual K/D=500/40, gain/rate bounds and target-following dynamics, compare 36 fixed K/D targets over 0.3 seconds. Best candidates are 100/15 and 30/15, but **using the same 30/15 gives mean excess cost only 0.4334% across the two cases**. This has no population-statistical meaning and does not establish important impedance benefits. No below-1% rejection rule is imposed; significant capability evidence is simply missing.

Cost explicitly weights force error, indentation-reference error and velocity squares, without preferred-K/D labels; see source/JSON for coefficients/all candidates. Against the known-state 36-candidate reference, the observer selects the best candidate throughout these noise samples, not a continuous-action global optimum. Energy filters, actuator saturation, real delay, contact loss, tangential scanning and sufficient direct force feedback are absent. All candidates retain positive indentation/force, so compression equations are not extrapolated into separation.

Independent DOP853 checks on both best candidates give maximum relative cost difference about 7.92e−10; halved-step difference about 7.37e−10. Holding also has an analytical exponential solution. Future noise is not used for prediction. The 100 repeats condition on two histories, not 100 independent tissues.

## Retained use

Reuse the known-model observer, realizable history construction and task-cost checks to distinguish prediction differences from impedance benefits. No memory network or mandatory new experiment is selected. Ordinary observers are legitimate tools; sufficiency does not prohibit research in this direction.

Reproduce with `OPENBLAS_NUM_THREADS=1 python experiments/contact_history_mechanism_audit/check.py`, using NumPy/SciPy. [result.json](result.json) retains hashes, parameters, all 36 costs and independent integration errors. This does not inherit B0 contact, HJB, energy or certification proofs.
