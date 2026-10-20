> Existing local validation. Obsolete topics/plans were removed; original configurations/results remain. See [HANDOFF](../../project_records/HANDOFF.md).

# Candidate revision1: value-root tracking under constraint activation

2026-10-13. **Untested candidate revision, not a replacement main method or established leading-journal/independent novelty.** Motivated by actual smooth-KKT condition failure at weakly active constraints; see [minimal test](results.json).

## 1. Why revise the question/method?

Slow parameters do not imply fixed active sets. Optimal gain rates may leave/enter rate or energy bounds. Differentiating the complete critic/action/multiplier root fails at potentially singular simultaneous zero multipliers/slacks; original implicit-sensitivity/C² neighborhoods need not cover these points.

Existing ENE handles piecewise active changes, and semismooth KKT/critical-cone methods remain strong candidates, not failed alternatives. Ask: **Can strongly convex impedance Hamiltonians support a more regular value-coefficient root, retain prediction accuracy/action feasibility through some activations, and reduce maintenance cost against these strong methods?**

## 2. Minimal mechanism revision

Retain gain-rate actions, actual gains and energy accounting. Eliminate actions/multipliers from the outer root:

\[
m_x(p)=\min_{u\in\mathcal U(t,x)}(u^TR_gu+p^Tu),\qquad
p=G_g^T\nabla_xV_W.
\]

Use minimized HJB collocation/terminal equations to define critic-only \(\widetilde F(W,\eta)=0\), and apply NEOC to W.

Low-level action recovery is

\[
u_W=\Pi_{\mathcal U(t,x)}^{R_g}\!left(-\tfrac12R_g^{-1}G_g^T\nabla_xV_W\right).
\]

This is original Hamiltonian recovery, not an added prediction-acceptance pipeline. Current state/gains/energy determine allowed actions each cycle. Covered parameter neighborhoods require no candidate-residual/rollout acceptance; predeclared errors/ranges manage anchors.

Initially restrict changes to material/model parameters, with references, action cost and same-state allowed sets independent of them. State/time variation remains in execution. Parameter-dependent references, sets or Rg require separate set-perturbation/Lipschitz arguments.

## 3. Direct conditional theory

Assume nonempty closed convex U, positive-definite Rg, bounded sets over collocation domains and sufficiently smooth coefficients/bases. Strong convexity gives unique minimizers,

\[
\nabla_pm_x(p)=u_W,
\qquad
\|u(p_1)-u(p_2)\|_{R_g}
\le\tfrac12\|R_g^{-1/2}(p_1-p_2)\|.
\]

so minimized values are C¹ in p with Lipschitz derivatives, without differentiable actions at every switch. Roots may be C¹,¹ rather than C². **These established smoothness/projection tools are not themselves new contributions.**

If \(\|\widetilde F_W^{-1}\|\le\kappa\),\(\|\widetilde F_\eta\|\le B\) throughout the neighborhood and joint derivatives have Lipschitz L, root sensitivity admits

\[
C=-\widetilde F_W^{-1}\widetilde F_\eta,\quad
\|C\|\le c=\kappa B,\quad
\operatorname{Lip}(C)\le M=\kappa L(1+c)^2.
\]

This requires the root graph/connecting paths inside the neighborhood and follows local C¹ implicit functions/derivative differences. For ideal anchors/sensitivities,

\[
\|W^*(\eta_a+\Delta\eta)-W_a-C_a\Delta\eta\|
\le\tfrac12 M\|\Delta\eta\|^2.
\]

Add anchor, sensitivity, numerical and estimation errors with original multistep accounting. Projection Lipschitz bounds propagate action errors through active changes. Every recovered action lies in U; original continuous dynamics/(31) preserve gain/energy boundaries.

**Not automatic:** critic domain accuracy, outer-root regularity, physical/contact feasibility, time-varying global optimality, numerical accuracy or total speed. Elimination cannot repair bad bases, missing nodes or model error. Original state-boundary issues remain.

### 3.1 Hamiltonian performance through active changes

At common actual state/U, let linear terms be p, p̂ and recovered actions u(p), u(p̂). Strong-convex quadratic/closed-convex conjugacy gives

\[
0\le H_p(u(\hat p))-H_p(u(p))
\le\tfrac14\|R_g^{-1/2}(\hat p-p)\|^2.
\]

Proof:φ(u)=uᵀRg u+I_U(u) has conjugate with Lipschitz derivative in the Rg metric; the left side is its Bregman remainder between−p,−p̂. Fixed active sets are unnecessary.

If \(\|R_g^{-1/2}G_g^T\nabla\phi\|\le L_p\), critic error εW gives Hamiltonian gap≤Lp²εW²/4. Ideal εW≤M||Δη||²/2 yields gap≤Lp²M²||Δη||⁴/16. **Fourth order concerns local prediction-induced excess only, not configuration/estimation error or total changing-environment task cost.**

For a common frozen problem with feasible in-domain trajectories and HJB/terminal residual boundsρ,τ, verification bounds continuous-optimal gap by2(Tρ+τ)+∫extra Hamiltonian gap. Time-varying η(t) additionally needs ∫∂ηV·ηdot plus model/estimation/execution errors; frozen fourth-order remainders are not whole-task fourth-order optimality.

## 4. Completed falsifiable tests

Single-node Hamiltonian root

\[
F(w,\theta)=1-2\theta w+\min_{|u|\le0.4}(u^2+2wu)=0.
\]

At critical θ=1.05, the complete3-variable KKT Jacobian has rank2 and minimum singular value≈9.9×10⁻¹⁷; reduced Fw=−2.9, with absolute lower bound2.6 over the test domain. First w derivatives are continuous; second derivatives jump, showing C¹,¹ without C².

| Path | Queries/anchors | Maximum critic-root error | Preset bound | Unprojected action violations | Recovered value-update violations |
|---|---:|---:|---:|---:|---:|
| Bound to release | 301/4 | 0.001873 | 0.002 | 0 | 0 |
| Free to constrained | 301/4 | 0.001619 | 0.002 | 23 | 0 |

Projected direct-action NEOC has release-path error≈0.00680, versus≈0.00153 after value updates. Entering bounds gives equivalent projected predictors, a retained no-advantage result. Four initial states versus independently analytical constrained scalar-policy cost give worst excess≤2.63×10⁻⁶. Single-node validation is not full R15 proof.

The [multinode test](multinode.py) uses16 nodes, a smooth-derivative cardinal-cubic critic,121 queries/9 anchors, and maximum derivative error≤0.002 across the same critical parameter. Off-node independent true-policy cost differences reach2.79×10⁻⁸, configuration differences4.01×10⁻⁹; reduced condition number≤14.34, while full complementarity remains singular. **Projected-action and value predictors have equal maximum action error here.** This supports regularity/feasibility, not general performance advantage. It remains scalar force control, not four-gain-rate/energy HJB.

## 5. Candidate status

Retain conditional theory/completed checks only. This is not merged into B0, a new main research line or default continuation task. Future question/method review determines adoption.
