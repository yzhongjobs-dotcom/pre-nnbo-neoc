# Contact real data — collection note

All files here are local research inputs. `literature/manual_inbox/` is ignored by Git, so raw downloads are not submitted.

## Obtained

- `static/Acquisition informations.tab` — Dahmani et al., Borealis dataset 10.5683/SP3/ASTGWY. 36 indentation levels from 0 to 7 mm in 0.2 mm steps, six acquisitions, voltage and calculated probe–phantom force (N). The accompanying `docs/Dahmani_README.txt` describes a gelatin phantom with a cylindrical inclusion; ultrasound image and RF files are separate and were not downloaded. Force at the nominal 0 mm row is already about 6.10–6.54 N; across the table it ranges about 5.23–16.00 N. Resolve this offset/preload before treating it as a zero-contact force curve. No time or velocity series is present.
- `static/Data Sheet 1.CSV` and `static/Data Sheet 3.PDF` — Remus et al. 2024 S1 coefficient/fit-quality table and S3 plotted experimental force–displacement curves. Filenames are preserved exactly as received. S3 is not a raw time series. S4 was not present in the folder.
- `dynamic/cyclic_relaxation_corona_radiata.csv` and `dynamic/cyclic_relaxation_visual_cortex.csv` — Greiner et al. files received; not yet inspected beyond filename/source match. The record describes human brain white/gray matter, cyclic loading and compression/tension relaxation.
- `dynamic/repository_test_indentation.xlsx` — Callera et al. file received; temporal, force and displacement contents have not yet been checked.
- `dynamic/Rheometer data.zip` — retained without extraction as requested. Categorized as dynamic from its filename; dataset provenance has not yet been confirmed.
- `scanning/USForce_episode_000000.parquet`, `USForce_info.json`, `USForce_episodes.jsonl`, `docs/USForce_README.md` — one representative 29.8 MB parquet episode plus metadata. The dataset reports 90 episodes, 15,527 frames, 30 FPS. Fields include ultrasound/base/wrist images, robot pose, force/torque, Cartesian action, task and timestamp. Episode 0 is left-kidney long-axis scanning (125 frames); sample timestamps advance by 1/30 s. The dataset card does not state force/pose units or calibration. State force/torque is zero at the first frames and nonzero later. The 30 FPS stream cannot resolve 20.8 Hz contact phase, and robot motion is not tissue indentation without surface/contact registration.
- `docs/Dahmani_README.txt` and `docs/USForce_README.md` preserve the source descriptions.

## Still needs manual download

- Neumann et al. 2024: `Figure 14 - Dynamic mechanical analysis.zip` (19.3 MB), `Figure 13 - Cyclic loading test.zip` (15.8 MB), `Figure 7 & 8 - Loading rate sensitivity test.zip` (14.9 MB). https://zenodo.org/records/14536174 . Zenodo record metadata was readable; file downloads returned 403 through both its content API and standard file links. Dataset is spinal-cord indentation, not abdomen.
- Remus S4: `DataSheet4.PDF` is optional and was not present in the folder. Publisher supplement page: https://www.frontiersin.org/articles/10.3389/fbioe.2024.1384062/full#supplementary-material . It contains surface-displacement means/ranges.
- MULTIS in-vivo indentation: user DOI https://doi.org/10.18735/S5R97D returned 404. The SimTK project https://simtk.org/frs/?group_id=1032 lists `002_MasterList_indentation.csv` and `001_MasterList_indentation.csv`; both download links redirect to SimTK login. Do not combine these with the separate postmortem tissue paper https://pmc.ncbi.nlm.nih.gov/articles/PMC6962198/.
- Boyer et al. 2009 full text: https://doi.org/10.1111/j.1600-0846.2008.00324.x . Wiley returned 403 and the university-hosted PDF link found during source check timed out. PubMed abstract is available at https://pubmed.ncbi.nlm.nih.gov/19152580/ . It reports in-vivo skin tests on 46 people, 10–60 Hz and 1–10 μm displacement amplitudes, Kelvin–Voigt behavior, and a 10 Hz complex-modulus range of 7.2±2.1 to 10.7±2.6 kPa across age groups; this is a boundary on small-amplitude skin response, not abdominal data or a raw series.

No complete imaging/RF archive or 2.8 GB Neumann archive was downloaded. US-PPNR was optional and left out; USForce already supplies a small robot-scanning sample.

## Reuse of repository literature and implications for the current plant

Existing PDFs were found and reused: Beber, Wijayarathne, AuSoScan and Moreira; nothing was downloaded again for these papers. They are model/method sources, not independent measurements.

`experiments/rossi_failure_study/scenes3.py` embeds Remus Table 4 mean force–depth polynomial coefficients for relaxed/activated abdomen measured with a 10 mm hemisphere. `physics_plant3.py` transfers that curve through `probe_curve`: equal `depth/sqrt(area)` state and a stiffness area exponent `gamma` (0.5 for half-space, 1 for a thin layer), then maps the reference response via `curve_table`. Remus S1/S3 can check the fit and region/state variation; S3 is plotted experimental data, not a raw time series. The half-space/thin-layer geometry transfer remains a modeling assumption until checked against actual-probe data. Dahmani is the closest direct probe geometry check, but its phantom and unexplained ~6 N zero-depth force prevent direct abdominal-parameter validation. Greiner/Neumann/Callera can support independent dynamic-mechanism tests; spinal cord and brain/pancreas properties must not be transferred as abdomen parameters. Boyer and USForce constrain context/ranges only, with the limitations above.
