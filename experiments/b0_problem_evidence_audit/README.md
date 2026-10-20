# Conditional audit of B0 problem premises

This is a theoretical/scale check, not a new controller-performance experiment. Interpretation and research-level judgments are in [field state](../../research/field_state.md), section 5.

`preview_scale_audit.py` uses B0 appendix B.2, D.1 and G.2 and retains three results: partial-coverage signals and local condition numbers with known/unknown boundaries; geometry before center arrival; and contact sizes required for fixed synthetic geometry under different loads. Independent quadrature checks the closed-form signals.

Reproduction: install root NumPy/SciPy dependencies and `matplotlib==3.9.4`, then run:

```bash
python experiments/b0_problem_evidence_audit/preview_scale_audit.py
```

The script needs only NumPy/SciPy/Matplotlib and does not depend on the removed old PPO directory.

Outputs are `results.json` and `footprint_information.png`. Plot moments are divided by half-width. Information condition numbers use equally normalized independent noise to illustrate structure; actual moment/force noise ratios require calibration.

The second source/JSON signal preserves sensor-coordinate signs: it is the normal-force contribution to M_y, the negative of appendix first moment m_1. The plot displays its magnitude.

![Footprint signals and local condition numbers](footprint_information.png)

The horizontal interval 0 to 1 runs from leading-edge contact to probe-center arrival at the material boundary; it does not mean the whole patch has entered the new material.

This example assumes known half-width and original-region material and separability of the normal first moment from total moment. It does not implement full MHE with unknown friction/surface motion/bias, calibrate force/moment sensors or measure tissue response. The load table concerns applicability of synthetic scan geometry, not clinical force limits or validation of the spatial-pressure approximation. B0 already states the degeneracies; they are not new contributions or research questions here.
