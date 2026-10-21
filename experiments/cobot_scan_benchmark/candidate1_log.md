# Candidate 1 design log (train/dev only until the single test run)

Each step below was checked on the T3 train cohort (6 or 20 scenes, seed 1000), with a few T1 train checks.
About 45 exploratory configurations were tried this way before the declared tuning. None of this used the test split.

1. **First version.** Estimator plus a 4-variable OCP (log D, KPF, u0, u1).
   - The feedforward was a linear ramp over the whole horizon. The cost was nominal only, from the current
     state, with a weak gain prior.
   - Result: train J 0.8-58, and u_ff saturated at 20 mm.
   - The stiffness estimate itself was good from the start: k_hat/k_true median about 1.0, p10-p90
     0.85-1.25, on T3.
2. **Diagnosis 1 (ki = 30).** The force-integral state wound up during a push. The OCP then held the force
   with a ramping u_ff while yc ramped the other way: a feedforward/integrator null direction that no
   short-horizon cost restrains.
   - ki = 0 removed the drift but the rule became poor (J 1.16).
   - Adding a terminal admittance-velocity term (w_v) made it worse.
3. **Isolate the information.** OCP with a very strong rule prior, gains only, which is effectively
   MHE estimate + B3b rule:
   - J 0.12-0.17 on 6 scenes, 0.205 on 20 train scenes.
   - For comparison on train: B1 0.384, EtR-A 0.370, oracle B3b(0.1 s) 0.137.
   - So the estimator carries most of the value.
4. **Feedforward still hurt** (0.15-0.26), even with the true surface velocity given to the predictor as an
   oracle diagnostic (J 0.87-15).
   - The predictor itself was checked open-loop against the plant and is accurate.
   - Cause: a horizon-long u_ff ramp extrapolated past the predicted surface motion, giving a slow
     (~0.75 s) null-mode oscillation of u_ff against yc.
   - Fix: ramp-then-hold feedforward, u0 + u1 min(i Tc, t_u). On a clean push scene J went 0.149 -> 0.009,
     with estimated drift (0.0088) about as good as oracle drift (0.0075).
   - T3 train: 0.205 -> 0.169.
5. **Free gains hurt** (pD = 10: 0.32; pD = 1: 0.45). The deterministic nominal prediction rewards
   aggressive error correction, and u_ff can cancel the predicted ringing.
   - A probe that adds a disturbance at the current state made it worse: it rewards fast rejection.
   - Fix: equilibrium robustness probes, an unpredicted surface ramp at k_hat and at 2 k_hat with no
     feedforward. With these the gain freedom gives 0.146-0.157, and even weak priors (pD = 0.1) work.
6. **Declared tuning.** 49 configurations per method variant and cell: the default plus 32 random draws
   plus 16 one-parameter perturbations, same draws for every variant.
   - The neighbouring update (as specified) was selected at train J 0.210 on T3 versus 0.144 for the
     ordinary 1-iteration solve; on T1 both were 0.025.
   - **Decision before any test run:** final solver = ordinary warm-started GN, and the neighbouring update
     becomes ablation (c). Ablations a, b and d were re-tuned with the final solver; the earlier NU-based
     ablation tunings are kept in the JSON as superseded and were never run on test.
7. **Test** was run once per variant on T3 and T1, and once for C1 on the secondary fast-push cell.
8. **Post-test checks** (off the test split):
   - Known-behaviour N1 found a start-up failure on 1 MPa tissue: the soft k0 prior gives soft gains at
     t = 0. Setting k0 = 3e4 fixes it and is J-identical on T3/T1 train and dev, so it is recommended
     (v2); it was not re-tested.
   - Fast push: 2 of 40 test scenes diverge after contact loss. A contact-loss fix (lossfix) made dev
     worse (5.8 -> 92) and was rejected. This failure is unresolved.

**Things that did not help:** ki = 0; terminal velocity cost; longer horizons for the feedforward
(H = 32-48 were unstable); current-state disturbance probes; free gains without probes; the neighbouring
update (no gain in decision quality, rollouts or time); lossfix.
