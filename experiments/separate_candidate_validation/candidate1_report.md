# Candidate 1: contact-model-driven gain-rate optimization and neighboring updates

**Current judgment: a complete leading-journal method/evidence is not established.** The original question remains; this is a reduced implementation/iteration, not full B0 or the entire candidate concept.

Reduced MHE → same-model finite-horizon OCP → directional-residual/cached-curvature prediction with limited correction → constrained actual gain rates. Direct OCP temporarily replaces an unimplemented critic root without inheriting node-root certificates. Sources are predictive VIC, ENE and quasi-Newton methods; no new general solver is claimed.

Expensive full sensitivity refresh was replaced by cached curvature/directional updates. Distant-target dead zones motivated B0 gain-rate actions; prediction/execution now align at 100Hz sampling/smoothing. Compatible-model averaging and minimum-residual single models were separately checked; extra cost with local changes did not justify default adoption.

## Independent results

One new 24-scenario cohort is reported at three tradeoffs, not 72 independent samples. Positive values improve composite cost, not necessarily each force/motion component. Unverified nonzero-force NNBO inversion comparisons are excluded from contribution evidence. E41 is only the declared OS scan adaptation.

| Weight | Improvement over adapted E41/% (95% interval) | Improvement over ordinary warm-start MPC/% (95% interval) |
|---|---|---|
| 1 | -3.959 [-6.883, -1.397] | -2.614 [-5.243, -0.491] |
| 12.5 | 3.681 [-3.480, 10.547] | -3.259 [-7.494, 0.117] |
| 25 | 2.271 [0.438, 4.374] | -0.039 [-0.624, 0.517] |

Main independent results assume zero release wait at design frequencies: model12.5Hz, cheap feedback/network100Hz. Ordinary MPC is not assumed limited to12.5Hz. Six validation cases separately check single-task100Hz requests with measured waiting, including front-end/projection expenses, not robot certification.

- Prototype w=1: median 1.399ms, p95 2.872ms;41 pending updates skipped among3600 requests. Front-end cost recorded separately.
- Prototype w=12.5: median 1.400ms, p95 3.101ms;33 pending updates skipped among3600 requests. Front-end cost recorded separately.
- Prototype w=25: median 1.411ms, p95 3.310ms;33 pending updates skipped among3600 requests. Front-end cost recorded separately.

## Limits and next revision

This carrier omits B0 positive-definite Qg/Rg, terminal costs, full trajectory constraints, strictly positive gain lower bounds and two-region Hertz front ends, with no platform calibration/cross-model/hardware generalization. Projection checks or budgeted QP references do not prove real safety/global near-optimality.

Restore complete task relationships and study how new contact information changes dynamic impedance/timely execution. Retain ordinary solving when fast enough, without invented bottlenecks. Establish reliable model-to-action effects rather than framing finer local solvers as the main contribution.
