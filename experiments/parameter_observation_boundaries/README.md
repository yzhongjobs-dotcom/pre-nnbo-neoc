# Parameter and observation boundary checks

Two independent tools check the candidates' shared front end and update assumptions; neither implements full MHE.

- `state_model_discriminator.py` and `state_model_discriminator/results.json`: distinguish known spatial-state motion, true material changes and estimator/model revisions, showing that slow physical changes do not automatically imply small descriptor updates.
- `contact_equivalence_audit.py` and `contact_equivalence_audit/`: recheck the material/surface equivalence degeneracy already present in B0 under the complete wrench model used here, retaining independent integration and identical-history arrays.

These check model conditions, not a new field-level failure, principal method or closed-loop advantage. Parameters/results are unchanged. Rerunning scripts overwrites matching outputs; use a copy or the correct Git version. Main dependencies are NumPy/SciPy; see source and B0 appendix for equations.
