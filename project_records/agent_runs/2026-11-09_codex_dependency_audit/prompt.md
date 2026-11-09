Read-only dependency audit in E:/github_project/scan-impedance-update. Do NOT edit, delete, run simulations, commit or push. Do not read research conclusions beyond what is listed.

Folder: experiments/rossi_failure_study (call it S). Use grep/AST scripts; run nothing that simulates.

Live roots (code that must keep working): S/corrected_comparison.py, S/keff_rescore.py, S/simple_ls_validation.py, S/rib_transition_eval.py, S/memory_segment_eval.py, S/overlay_compare.py, S/baseline_compare.py, S/id_baselines.py, S/ctrl_baselines.py, S/id_layer_root_cause.py, S/low_estimate_forensics.py, S/plant_arm_check.py, S/plant_arm_effect.py.

Task 1. Transitive import closure of the live roots within the repo (include experiments/motion_aware_scan and others on sys.path; check how sys.path is set). List the closure's .py files in S, and the .py files in S NOT in the closure.

Task 2. For every entry directly under S/results (file or folder), find every reference from any .py in the repo, including indirect forms: string literals, f-strings, glob patterns, B.load/ORIGINAL_LOAD/load/dump names, `OUT/...`, `HERE/'results'/...`, stored `log_file` and `samples_file` paths inside JSON. Mark whether each reference is (a) from a closure file on a path that runs in normal use of the live roots, (b) from a closure file only in a stage or function not used by corrected_comparison/keff_rescore/simple_ls_validation, or (c) only from non-closure scripts. Pay attention to: overlay_compare.chosen() → baseline_test_r2.json; corrected_comparison.legacy_check → overlay_test.json + a baseline_logs npz; memory_closed_loop, memory_filters, memory_filters2; h_fix_eval. Give sizes.

Task 3. Non-Python references. Find which results files, and which of these specific stale items, are cited by these entry docs, with line numbers: research/identification_module_scheme_20261108.md, project_records/HANDOFF.md, experiments/experiment_index.md, research/disturbance_model.md, research/identification_field_state.md, research/method_reassessment_20261108.md, research/paper_draft_zh.md, research/paper_experiments_final_zh.md, S/baseline_interpretation.md, S/overlay_interpretation.md, S/baseline_protocol.md, README files. Stale items: baseline_interpretation.md, overlay_interpretation.md, baseline_summary*, baseline_module*, baseline_test*, overlay_*, memory_closed_loop_*, h_fix_eval.json, memory_segment_eval.json, rib_transition_eval.json, memory_root_causes.json, low_estimate_forensics.json, gain_truth_ceiling*, id_layer_root_cause.json. Also check whether any .py or doc outside S references S/results.

Task 4. Which files in S/results are tracked by git vs ignored (`git ls-files`, `git check-ignore`), and the git-tracked size total.

Output: tables only, under 2 pages, no recommendations beyond a final column "no live dependency: yes/no". Do not re-derive science. Skip anything not listed.
