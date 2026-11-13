**The requested signal is not clean enough in this harness: L1–L4 failed; L5 passed.** Exactly 20 episodes ran, with unchanged estimator/controller and no new feedback signal.

Repository writes were refused. Outputs are in `D:\pip_env\temp\lateral_dither_20261113_144358`: [script](D:/pip_env/temp/lateral_dither_20261113_144358/experiments/rossi_failure_study/lateral_dither_check.py), [complete JSON](D:/pip_env/temp/lateral_dither_20261113_144358/experiments/rossi_failure_study/results/lateral_dither_check.json), and figures for [36](D:/pip_env/temp/lateral_dither_20261113_144358/research/figures/lateral_dither/case_36_a0.10.png), [39](D:/pip_env/temp/lateral_dither_20261113_144358/research/figures/lateral_dither/case_39_a0.10.png), [12](D:/pip_env/temp/lateral_dither_20261113_144358/research/figures/lateral_dither/case_12_a0.10.png). Existing repository files were untouched; nothing was committed or pushed.

The plant interpolates `xcmd` on its 1 ms grid and drives probe position `s[1]` through a first-order lag. Contact position follows through tangential spring–mass dynamics, damping \(2ζ_t\sqrt{K_tm_t}\), friction, and axial/lateral compliance coupling; it does not rigidly follow the encoder. The commanded velocity bound is violated in the slow cases, although filtering prevents contact reversal at 0.1 mm; **0.2 mm causes contact reversals in cases 38 and 12**.

Throughout the tables, entries follow amplitudes **[0, 0.05, 0.1, 0.2] mm**, separated by semicolons. **P/C** means plateau/changing ticks; **—** means unidentifiable or no qualifying ticks.

**1. Achieved motion.** Each entry is probe/contact. Commanded velocity/scan-speed ratios are `[0,.213,.426,.851]` for 36/39/0, `[0,.574,1.149,2.297]` for 38, and `[0,.561,1.122,2.243]` for 12.

| Case | Achieved amplitude, mm | Achieved velocity/scan-speed ratio |
|---|---|---|
|36|0/.0005; .022/.019; .045/.038; .089/.077|0/.002; .095/.081; .190/.163; .380/.328|
|39|0/.0003; .022/.019; .045/.038; .089/.077|0/.001; .095/.081; .190/.163; .380/.328|
|38|0/.0002; .027/.032; .053/.067; .106/.107|0/.003; .305/.370; .609/.773; 1.219/1.227|
|0|0/.0003; .022/.019; .045/.038; .089/.077|0/.001; .095/.082; .190/.164; .380/.329|
|12|0/.00008; .023/.028; .046/.056; .091/.109|0/.001; .255/.315; .511/.631; 1.022/1.223|

**2. Force ripple and cost.** Ripple is N rms, P/C. Both frequencies were fitted jointly in local 24-tick windows with quadratic detrending.

|Case|15.625 Hz P/C|20.833 Hz P/C|L|ΔL|
|---|---|---|---|---|
|36|.140/.093; .143/.235; .129/.379; .198/.764|.135/.161; .145/.171; .108/.143; .159/.213|.790;1.158;1.369;1.725|0;.368;.578;.935|
|39|.033/.100; .044/.130; .055/.149; .066/.156|.076/.091; .084/.089; .086/.085; .079/.064|.082;.106;.141;.146|0;.024;.059;.064|
|38|.052/.072; .082/.185; .135/.331; .206/.518|.101/.122; .094/.062; .087/.057; .073/.115|.106;.201;.448;.903|0;.095;.342;.797|
|0|.029/—; .034/—; .041/—; .062/—|.091/—; .093/—; .093/—; .089/—|.056;.062;.064;.082|0;.006;.008;.026|
|12|.031/—; .034/—; .041/—; .057/—|.109/—; .109/—; .107/—; .102/—|.0955;.0959;.0961;.0983|0;.0004;.0006;.0028|

**3. Gradient agreement.** Zero-lag entries are correlation/κ; optimized entries are correlation/κ/lag in ms. Positive lag means the signal follows truth.

|Case|Zero lag|Best lag within ±200 ms|
|---|---|---|
|36|—; .592/2.200; .533/2.412; .759/4.027|—; .750/2.822/80; .556/2.493/72; .769/4.103/16|
|39|—; .105/.191; .165/.190; .325/.204|—; .121/.223/72; .173/.201/32; .360/.229/72|
|38|—; .858/3.259; .872/2.266; .868/1.783|—; .866/3.279/40; .873/2.260/16; .869/1.774/24|
|0|—; .317/3.560; .501/3.682; .601/3.407|—; .324/3.617/16; .519/3.778/48; .653/3.664/64|
|12|—; .126/1.220; .368/2.023; .592/2.335|—; .190/1.854/−104; .430/2.382/−104; .636/2.535/−112|

**4. Plateau noise.** First column: empirical/predicted sd of \(c\), N/m. Second: implied lateral/V2 rate noise, /s, using the offline best-lag κ.

|Case|sd(c): empirical/predicted|Rate noise: lateral/V2|
|---|---|---|
|36|—;1465/697;860/344;525/172|—/2.65;1.51/1.86;1.00/.71;.37/.46|
|39|—;642/693;396/344;248/178|—/3.15;8.29/3.14;5.67/3.02;3.13/3.25|
|38|—;1040/590;759/296;520/151|—/2.78;.34/2.82;.36/2.37;.31/1.79|
|0|—;699/694;413/348;319/177|—/3.24;.56/3.27;.32/3.42;.25/3.54|
|12|—;516/678;292/339;210/170|—/3.30;.30/3.30;.13/3.33;.09/3.50|

**5. Stiffness error.** Median \(|\ln(\hat k/k_{\mathrm{eff}})|\), joint/single-regressor fit. At zero dither, the identifiable axial fit is retained while \(c\) is undefined.

|Case|Plateau joint/single|Changing joint/single|
|---|---|---|
|36|.064/.064; .046/.060; .103/.106; .116/.115|.205/.205; .289/.352; .301/.432; .300/.308|
|39|.098/.098; .122/.103; .100/.080; .098/.088|.146/.146; .151/.143; .154/.139; .215/.207|
|38|.094/.094; .099/.114; .117/.151; .166/.214|.112/.112; .285/.482; .283/.430; .197/.458|
|0|.096/.096; .106/.102; .110/.110; .117/.145|—;—;—;—|
|12|.087/.087; .086/.088; .086/.081; .087/.094|—;—;—;—|

**6. Case 12 body motion.** \(c/F_{24}\), units 1/m; event window is 2.000–2.993 s.

|Amplitude, mm|Inside mean ± sd|Outside mean ± sd|Mean difference / plateau sd|
|---|---|---|---|
|0|Unidentifiable|Unidentifiable|—|
|.05|12.91 ± 44.94|20.61 ± 64.46|.119|
|.10|.64 ± 22.85|18.69 ± 36.49|.495|
|.20|−7.46 ± 13.32|16.02 ± 26.21|.896|

The predictions were recorded before episodes:

- **L1 failed:** at 0.1 mm, best correlations were .556/.173 and κ was 2.493/.201 in cases 36/39.
- **L2 failed:** lateral noise was .999/5.670 /s versus V2’s .705/3.016 /s; neither achieved threefold improvement.
- **L3 failed:** added changing-tick ripple reached .286 N in 36 and .259 N in 38; ΔL exceeded .03 in 36, 39 and 38.
- **L4 failed:** overall error increases were .0447, .0042, .0138, .0211 and −.0011 for 36,39,38,0,12.
- **L5 passed:** the 0.1 mm event mean shift was .495 plateau sd; all nonzero amplitudes passed.

Important limits: without dither, a linear scan lies in D4’s nullspace—direct matrix subtraction initially manufactured spurious lateral information, corrected using the mathematically equivalent spectral form of M. Five joint stiffness fits were nonpositive in case 36 at 0.1 mm and are counted/excluded from logarithms. Case 38 has only 20 changing ticks at that amplitude; cases 0/12 have none outside events. Finally, \(d\ln k_{\mathrm{eff}}/dt\,/v\) includes depth and velocity effects, so these results do not isolate a fixed-depth spatial derivative or establish body-motion immunity.