# Candidate 2 design log (train/dev only; the test split was run once on the final design)

Dev numbers are mean J over the 20-scene dev cohorts (seed 2000) at the two cell doses.
Reference on dev: EtR-A (= the zero-output student) T3 0.210, T1 0.119, selection score 0.713
(score = mean over cells of J_dev / J_dev(B1)).

1. **Teacher tuning** (train cohorts, 72 rules). Best: TruthKFF, 0 s latency, (wI 80, cP 2, ki 30,
   Dmin 100, preview 0.05 s): train J T3 0.135 / T1 0.007 (B1 0.384 / 0.116). Same rule without
   the feed-forward: T3 0.136 / T1 0.011 (preview 0). A surface-motion feed-forward oracle
   (true surface height, 0.03 s preview) reached T3 0.012 but was not used: it bypasses the
   force loop's delay and noise with an unobservable signal, so its labels are not inferable.
2. **v1: absolute action head, fixed prior (M 1.5, D 120, K 0, ki 10, KPF 1.15e-4), u_ff
   up to 8 mm, teacher = TruthKFF.** DAgger students diverged in closed loop (dev J 100-400).
   Hybrid tests: with the teacher's u_ff and student gains J ~0.2-1.4; with the student's u_ff
   J ~100. The u_ff head, regressed from force history, forms a high-gain force-to-position loop
   (u_ff errors of ~1 mm are ~35 N on bone). Dropped the absolute u_ff head.
3. **v2: residual around the stiffness rule on the lock-in estimate (zero output = EtR-A,
   verified to 1e-11), independent D/KPF residuals, structured k_hat feed-forward
   c_ff*FD*(1/k_hold - 1/k_ref).** The structured feed-forward was harmful (R^2 of teacher
   u_ff on phi = -0.15; dev J 2-5 with it vs 0.1-1 without): 1/k_hat on soft tissue is noisy.
   Disabled (USE_FF = False). The residual-gain student was still unstable (T1 scenes J 40-70):
   gain chatter from tick-to-tick noisy features. A first-order output filter (tau 0.03 / 0.1 s)
   reduced but did not remove the failures.
4. **Copycat leak.** Features 26 (y_cmd - y) and 33/34 (applied D, KPF) copy the executed action;
   when the teacher executes in DAgger they reveal k_true. Masked permanently. DAgger dev then
   T3 0.4-2.0, T1 0.06-1.3, still worse than EtR-A on T3. Residual spans 0.35 and teacher latency
   0.1 s were worse (latency-0.1 teacher: T1 J up to 5.5, degrading with iterations).
5. **Diagnosis: can history infer the hidden state?** Regressing log k_true on the causal
   features (prior-policy rollouts, 48 train episodes, 16 held out): held-out log-RMSE MLP 0.36
   vs lock-in raw 0.54, peak hold 0.58 (bias +0.31), affine-corrected 0.49; at transitions 0.34
   vs 0.46. So the history does carry more stiffness information than the lock-in estimate.
6. **v3: distil the teacher's privileged state.** The teacher is the truth-stiffness rule, so
   its unsaturated action is a known function of k_true; the student's stiffness head
   (residual of log k around log k_hold, +-2 e-folds) is regressed on log k_true, the rule maps it
   to D and KPF. MSE regression: dev T3 0.54-3.4 (worse): it removed the conservative bias of the
   peak-held lock-in, and under-estimates of k on bone destabilise contact.
7. **v3 + expectile loss.** tau 0.8: dev T3 0.25-0.80, T1 0.06. tau 0.95: after 6 iterations
   T3 0.211, T1 0.047, score 0.458 (EtR-A 0.713). Adopted tau 0.95.
8. ES sigma: 0.03 on the 325 last-layer weights produced exploding population members (mean
   fitness 5-30 x B1); final runs use 0.01 and 50 generations.
9. **Final runs** (full, a, c, d, e x seeds 0-2; ~5 min each on 2 processes). Dev selection
   scores: full 0.429 / 0.428 / 0.455; a 0.651 / 0.543 / 0.673; c 0.460 / 0.428 / 0.448;
   d never feasible on the train peak rule; e feasible only for seed 0 (0.760). An earlier queue
   was aborted after ~10 min because numpy/OpenBLAS threads in the parent pushed CPU use above 2
   cores; rerun with single-threaded BLAS (identical code and seeds).
10. **Test (run once).** See candidate2_results.json and the final report.
