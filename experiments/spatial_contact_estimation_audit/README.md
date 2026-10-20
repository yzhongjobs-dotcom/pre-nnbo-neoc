# Conditional checks of spatial contact estimation and task scales

Two independent tools extracted from the old receding-policy branch. Source/results were moved unchanged; no new experiment was run. They check geometry/material confounding, noise propagation and simplified task scales in the shared front end, not a full MHE or paper method.

- `transport_probe.py`: a finite elastic foundation with known initial surface and tangential lever arm, an affine material field, and updates of material/global displacement from force and normal first moment. Includes sharp boundaries, convective-viscous mismatch and raw wrench noise. Raw Fn/Ft noise is 0.03N each and My noise 0.0001Nm; moment noise after lever-arm subtraction is about 0.000757Nm. The 5mm/0.2Hz condition comes from Moreira's original experiments, not project measurements.
- `task_geometry.py`: analytical cost checks for a frozen linear spring, quasistatic behavior and a shared world reference, explaining why weights make fixed gains nearly optimal. It is not an upper bound on nonlinear dynamic task benefits.

Under affine fields with raw wrench noise, displacement-estimation RMS is about 0.021–0.030mm; sharp-boundary/viscous mismatch gives about 0.058/0.086mm. These are estimates at synchronized sample instants, without compensation for 20ms sensor delay, closed-loop validation, demonstrated novelty or real-tissue validity. Full parameters/errors are in [results](transport_probe.json); per-step arrays are in `transport_probe.npz`.

If needed, run `python experiments/spatial_contact_estimation_audit/transport_probe.py` or `python experiments/spatial_contact_estimation_audit/task_geometry.py` with NumPy/SciPy. Runs overwrite matching outputs. This consolidation checks transfer integrity without rerunning them.

Pre-extraction version: `db28365f33dc314ce8bfe83a23796dd18cf27f5b`; original location: `实验/闭环策略滚动阻抗/`. Relevant sources are in the [literature index](../../literature/literature_index.md).
