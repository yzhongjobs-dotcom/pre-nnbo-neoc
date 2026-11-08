# Task: targeted fact check of the real scanning scene (no method work)

Project: robot ultrasound scanning at constant contact force on the human back or trunk, with only force and probe position available. We need to know how the scene physically happens, before any modelling. Read nothing else in the repo except what is listed here. Do not read project_records/HANDOFF.md.

Already known (do not repeat, only correct if contradicted): `research/scan_scene_evidence.md`. It covers abdominal surface motion (Quirk: AP 8±3 mm peak to peak, 4.2±1.1 s), rib covering 7–32 mm and rib gap 6–19 mm (Patterson, cadaver), and Remus S1 abdominal indentation.

Answer these questions. Each fact needs: value with spread, population and posture, how it was measured, source with DOI or URL and location in the source (page, table or figure), and grade [full text], [abstract only] or [not verified]. Say plainly when you found nothing. Do not infer values without saying so.

Q1. Breathing motion of the back surface. In prone, sitting or side-lying subjects: peak-to-peak displacement and direction, how much it varies from cycle to cycle (amplitude, period, baseline drift), waveform shape (inhale/exhale asymmetry, harmonics), quiet versus deep. How do amplitude and phase vary with position on the back (upper versus lower thorax, paraspinal versus lateral, over ribs versus between ribs)? Are there measured phase differences between positions? Useful sources: optoelectronic plethysmography, structured-light plethysmography, surface imaging in radiotherapy, prone imaging studies.

Q2. How the back or spine is actually scanned by robots and sonographers: posture, probe type and orientation, path (along the spine, across the ribs), speed, contact force, whether subjects breathe freely. Search local files first with grep only: `literature/extracted_text/` (spine, scoliosis, back), then the web.

Q3. Mechanics across a rib on the back. Posterior rib width and spacing and overlying tissue thickness by level in living subjects (imaging). Measured contact force-displacement or stiffness over a rib versus an intercostal space versus paraspinal muscle, in N/mm with indenter size, in vivo if possible. How stiffness depends on indentation depth. What is known about the shape and length of the stiffness transition along the skin.

Q4. Coupling between breathing and the mechanics. Does measured stiffness of the back or thorax change with breathing phase or lung volume (for example posteroanterior spinal stiffness studies, by how much, and in which direction)? Do ribs and overlying tissue move relative to a skin-contacting probe during breathing (rib kinematics, tissue sliding), by how many mm? Does probe force itself change breathing or the surface motion?

Q5. Muscle activation under the probe. Time course and size of stiffness change in back or trunk muscles when activated (cough, startle, voluntary contraction), from shear-wave elastography, indentation or EMG onset timing. Only measured values.

Q6. Any public dataset with force and probe position recorded during ultrasound scanning of a living back or trunk.

Constraints:
- Scene facts only. No method suggestions, no comparison of algorithms, no statement about what our simulation should do.
- Do not download paywalled papers. Do not ask for login. Do not edit any repo file except the output file. No commit, no push. Leave `00_基线快照_2026-10-14/` alone.
- Run no simulations.
- Stop when each question has its best available answer or a plain "not found"; about 20 targeted searches in total.

Deliverable: a report of at most 3 pages, in Chinese with English paper titles kept. For each question: a short table of facts, then one line on what remains unknown. End with a list of the 5 facts most uncertain or most surprising relative to the known list above.
