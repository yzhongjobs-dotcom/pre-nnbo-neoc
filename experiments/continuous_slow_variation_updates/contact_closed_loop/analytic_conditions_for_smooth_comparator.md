> Existing local validation record. Obsolete topics/plans were removed; original configurations/results remain. Current plans: [HANDOFF](../../../project_records/HANDOFF.md).

# Smooth-contact comparator: analytical conditions for direct NEOC

This simplified validation proof is neither a new full-dimensional R15 proof nor machine certification. It isolates slow parameters, direct prediction and error guarantees; scripts additionally check floating-point implementation, finite-time surface motion and sampled control.

## 1. Smooth critic root with an exact reference

With local normal stiffness θ, mass m and environmental damping b, freeze θ/reference. State x=(e, v) and force increment u satisfy

\[
\dot x=A(\theta)x+Bu,\quad
A=\begin{bmatrix}0&1\\-\theta/m&-b/m\end{bmatrix},\quad
B=\begin{bmatrix}0\\1/m\end{bmatrix}.
\]

The stationary regulation cost is ∫(xᵀQx+r u²)dt, Q=diag(q_e, q_v)>0, r>0. It penalizes **force increments** about equilibrium, not a finite infinite-horizon integral of sustained nonzero force squared. The quadratic single-critic V=xᵀPx root is

\[
F(P,\theta)=A^TP+PA-PBB^TP/r+Q=0.
\]

Feedback u=−ke−dv is directly normal PD impedance, where

\[
s=\sqrt{\theta^2+q_e/r},\quad
k=s-\theta,\quad
w=\sqrt{b^2+2mk+q_v/r},\quad d=w-b,
\]

\[
P=\begin{bmatrix}
r[bk+(\theta+k)d]&rmk\\rmk&rmd
\end{bmatrix}.
\]

Substitution checks Riccati. k, d>0 and characteristic polynomial mλ²+(b+d)λ+(θ+k) ensure frozen-loop stability. Here P₁₁>0 and det(P)/(r²m)=bkd+sd²−mk² ≥ s_min d_min²−m k_max²>0, so P≻0. Two stabilizing roots would differ by a homogeneous Sylvester equation between stable loops, whose unique solution is zero. The stabilizing root is unique/smooth, with sensitivity

\[
A_c^TC+CA_c=-(A_\theta^TP+PA_\theta),\quad
A_c=A-BB^TP/r.
\]

The left Lyapunov operator is invertible because A_c is stable. Compute C and check independent analytical-solution differences; predicted points need no new root solve.

## 2. Predeclare radius and release directly within it

Let c=q_e/r,θ∈[θ_l,θ_h]. Then

\[
k'=\theta/\sqrt{\theta^2+c}-1,\quad
0<k''=c/(\theta^2+c)^{3/2}\le M_k.
\]

Take L_k=1−θ_l/√(θ_l²+c), w_min=√(b²+2mk(θ_h)+q_v/r). Then

\[
|d'|\le mL_k/w_{\min},\quad
|d''|\le M_d:=mM_k/w_{\min}+m^2L_k^2/w_{\min}^3.
\]

Gains recover linearly from P, so P_p=P_a+C_a(θ−θ_a) exactly gives their first-order predictors. For ideal exact anchors/sensitivities,

\[
|k_p-k^*|\le M_k|\Delta\theta|^2/2,\quad
|d_p-d^*|\le M_d|\Delta\theta|^2/2.
\]

With ε_k=0.1 N/m and ε_d=0.002 Ns/m, preselect

\[
r_\theta=\min\{\sqrt{2\varepsilon_k/M_k},\sqrt{2\varepsilon_d/M_d}\}\approx13.0\ \mathrm{N/m}.
\]

Within this radius predict directly from anchors; create new anchors only outside it. Candidate residuals/true errors do not determine acceptance. Inexact anchors/sensitivities need additional budgets; this floating-point check does not interval-enclose numerical errors.

## 3. Slow-variation stability and limits

Use V=xᵀP*(θ(t)) from exact frozen roots, with gain deviation Δg=(Δk,Δd). In homogeneous error dynamics,

\[
\dot V\le-[1-\gamma_\theta-\gamma_g]x^TQx,
\]

\[
\gamma_\theta=\|Q^{-1/2}P_\theta Q^{-1/2}\|\,|\dot\theta|,
\quad
\gamma_g=2\|Q^{-1/2}PB\|\,\|\Delta g Q^{-1/2}\|.
\]

Analytical interval bounds include material sampling, causal interpolation over a release cycle and low-level half-step holds in Δg. γ_θ+γ_g<1 gives uniform homogeneous decay. Moving surfaces, references and sampled-reference deviations are separately bounded external inputs; excluding them cannot establish vanishing full-task error.

**The proof shows how to predeclare direct-NEOC validity without per-candidate rollout checks.** It proves neither global time-varying task optimality nor continuous force/energy certification, and excludes full R15 nonlinear contact, gain-rate OCP, MHE and preview. Scripts check those physical metrics retrospectively at samples.
