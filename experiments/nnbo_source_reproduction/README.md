> Existing local validation. Obsolete topics/plans were removed; original configurations/results remain. Current plans:[HANDOFF](../../project_records/HANDOFF.md).

# Minimal original NNBO mechanism reproduction

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent; PDF reading; nature-figure(Python plotting)
- Origin Mode: run + source-equation checks + reproducibility
- Origin Date: 2026-10-12
- Verification Status: REPRODUCED; four original trajectory/metric groups rerun identically, not exact reproduction of every figure or HJB-optimality proof
- Version Label: nnbo_source_minimal_v1
- Status: PASS_MINIMAL_MECHANISM_NOT_EXACT_FIGURE_REPLICATION
- Python 3.10.11，NumPy 1.26.4，SciPy 1.13.1
- Latest main run:16.54 seconds, exit0.

## 1. Completed scope

Implemented/executed the original 1D Hunt–Crossley force-input augmented system, three-term single critic, online weights and exploration. Both weights complete0–7s without contact loss/nonfinite states or weights. Larger position-error weighting produces positions nearer target and larger force.

Original(52) Newton Downhill inversion of recorded force/critic to position references is independently checked pointwise by bracketing roots. **It is not connected into a measured-force→reference→ideal-position-robot closed loop.**

Thus scope is **minimal §III-A–C/§IV augmented-system simulation and §III-D inversion checks**, not parameter identification, hardware,[44]/[45] comparisons or all figures. This is not current R15 five-state HJB–KKT root reproduction.

## 2. Source correspondence and parameters

Kong et al.,*Neural-Network-Based Optimal Impedance Control for Robots in Physical Interaction With Soft Environments*, IEEE TSMC: Systems55(10),6463–6475(2025), DOI10.1109/TSMC.2025.3579017.

Local source:[NNBO PDF](../../literature/papers/NNBO.pdf).

`results.json` records PDF/original-code SHA256; extracted text is in literature. Duplicate equation screenshots were removed; check the PDF directly.

| Source | Implementation |
|---|---|
| §II-D,(21),§IV |r_dot=-r+0.03, r0=0.08, xd=r|
| §III-A,(23)–(25) |x_dot=-Ke/Be+f/(Be(x-xe)^n)+ε; X=(e, r), e=x-r|
| §III-C,(35)–(37) |Jhat=W1 e²+W2 er+W3 r²; fhat=-(2W1 e+W2 r)/(2RBe(x-xe)^n)|
| (38)–(41) |Hhat=W^Tν+Qe²+Rfhat²+ζ²; W_dot=-κνHhat/(1+ν^Tν)²|
| §IV（PDF p7—8） | Ke=1149，Be=567，n=1.2，xe=0.01 m，x0=0.011 m；W0=(20,20,20)，κ=5；Q=500000/1000000，R=1 |
| §IV exploration/disturbance |Complete original trigonometric exploration added to force;ε=1000 e sin(r)|
| §III-D,(52)/Algorithm1 |Independent Newton Downhill inversion at fixed measured force, r, W, checked with brentq|

ν uses original nominal F+Gfhat; physical dynamics use explored fhat+N and stated disturbance. True disturbances are not inserted into nominal critic dynamics.

## 3. Unreported settings and reproduction assumptions

1. **Original ζ is unspecified.** Main runs use tight stated-disturbance bound|1000e sin(r)|; also run conservative1000|e|.Both satisfy|ε|≤ζ. Retain both, without selecting by visual similarity.
2. **Exploration end time is unspecified.** Apply the full original signal throughout7s, without truncation/amplitude reduction. Terminal oscillations preclude strict-steady-state claims.
3. The 7 s duration comes from the displayed range of Fig.2–3; outputs follow the original 0.001 s interval. Adaptive RK45 integrates the continuous learning law with maximum internal step 0.001 s. This differs from explicit Euler weight updates every 1 ms followed by holding the weights constant.
4. Shallow initial contact is stiff; fixed1msRK4 crossed the contact domain during development. Adaptive integration fixes steps without clipping force, weights, position or altering learning laws. Independent Radau agrees.
5. Force-input augmented integration instantiates §III-A–C, not the original inner robot-position loop. Retain this distinction in R15 transfer.

## 4. Measured results

Values are actual outputs at t=7s; applied force includes exploration. At that instant xd≈30.046mm and tissue surface10mm.

|Q|Disturbance-bound assumption|Finalx/mm|Finalappliedforce/N|Peakforce/N|FinalHhat|
|---:|---|---:|---:|---:|---:|
| 500000 | tight_known_disturbance | 22.468 | 6.636 | 14.330 | 62.938 |
| 1000000 | tight_known_disturbance | 24.708 | 7.849 | 18.126 | 79.754 |
| 500000 | uniform_sine_bound | 26.063 | 8.562 | 20.619 | 87.214 |
| 1000000 | uniform_sine_bound | 27.047 | 9.068 | 22.529 | 91.160 |

Under main assumptions, increasingQ changesx≈22.47→24.71mm and force≈6.64→7.85N, matching qualitative accuracy/force tradeoffs, not exact plotted curves.

**Hamiltonian residual does not approach zero; weights still change during the last second.** Stable execution/plausible interaction do not establish exact HJB critics or R15 root thresholds. Bounds substantially change policies, showing sensitivity to missing settings.

## 5. Numerical/interface checks

- Scalar policies/learning laws agree termwise with independent matrix forms.
- Maximum nominal HC-plus-disturbance/applied-force identity error≈7.11×10^-15N.
- Four runs produce7001 outputs each with positive indentation, without artificial saturation.
-4×701 inversions agree with independent brentq, residual<10^-10N, at most5 Newton steps.
- Independent Radau versusRK45 maximum position difference2.58e-10m, weight difference2.51e-07.
- Four originalt/y/data arrays and metrics rerun identically; runtime need not.
- Record finite7s excitation Gram matrices; positive eigenvalues there do not prove PE on arbitrary windows.

## 6. Costs, optimality and R15

CSV separates runningQe²+Rfhat², finite0–t accumulation, critic values andHhat. Fig.2(d)–3(d) show decreasing cost, but forward integrals of nonnegative running costs cannot decrease. Missing details prevent pointwise reconstruction; do not rename quantities to manufacture agreement.

Nonzero terminal tracking error/force prevent claims of finite infinite-horizon cost or solved infinite-horizon optimality. Finite simulation does not independently prove original bounded-weight arguments/PE assumptions.

This establishes a runnable reference for the original single-critic learning mechanism. Earlier R15 failure concerned initial solution of another fixed-collocation root; **running the original simulation does not automatically prove R15 root solvability**. The original learning mechanism and B0 fixed root differ in states, actions and representations, so its three-dimensional weights cannot be transferred directly.

## 7. Files and reruns

- [Driver](run_reproduction.py),[frozen configuration](config.json),[complete results](results.json).
- Retain four complete`Q*_*.npz`trajectories/results.json. OldCSV/inversion tables remain recoverable from Git or driver.
- [Reproduction check](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/NNBO%E5%8E%9F%E6%96%87%E6%9C%80%E5%B0%8F%E5%A4%8D%E7%8E%B0_2026-10-12/reproducibility.json),[first run](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/NNBO%E5%8E%9F%E6%96%87%E6%9C%80%E5%B0%8F%E5%A4%8D%E7%8E%B0_2026-10-12/run1/results.json),[run log](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/NNBO%E5%8E%9F%E6%96%87%E6%9C%80%E5%B0%8F%E5%A4%8D%E7%8E%B0_2026-10-12/run_log.txt).
- [Plot script](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/NNBO%E5%8E%9F%E6%96%87%E6%9C%80%E5%B0%8F%E5%A4%8D%E7%8E%B0_2026-10-12/plot_diagnostics.py): Matplotlib was absent and installation incomplete; no diagnostic plots were generated or inspected. Numerical reproduction does not require it; CSV/NPZ trajectories remain.
- [Summary/rerun comparison script](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/NNBO%E5%8E%9F%E6%96%87%E6%9C%80%E5%B0%8F%E5%A4%8D%E7%8E%B0_2026-10-12/summarize_reproduction.py).

Run in this directory:

```powershell
python -X utf8 run_reproduction.py
```

First-round duplicate copies were removed; Git retains history. Numerical runs require NumPy/SciPy; rerun in a separate copy/output location.
