**The memory fixes case 7, but the full hypothesis fails.** K2 and K4 pass; K1, K3, K5, K6 and K7 fail. No tuning, extra episodes, existing-file edits, commits or pushes.

Sandbox writes required this mirror: `D:\pip_env\temp\memory_transport_mirror_99ou0lkl`. [Code](D:/pip_env/temp/memory_transport_mirror_99ou0lkl/experiments/rossi_failure_study/memory_transport_loop.py) · [Results JSON](D:/pip_env/temp/memory_transport_mirror_99ou0lkl/experiments/rossi_failure_study/results/memory_transport_loop.json). Verified **168 new runs, 84 reused rows and exactly 24 requested traces**; runtime 315 s.

Throughout, columns are **H=huber, D=delag2, K=kf, KL=kf_long, NT=kf_noT, BK=band_kf**. `As` denotes A_slow_thin; `Lo/Lf` original/fresh A_large.

Mean/median L:

```text
Group       H           D           K           KL          NT          BK
A       .106/.078   .106/.078   .097/.070   .098/.072   .097/.070   .490/.066
As      .591/.158   .591/.156   .146/.155   .167/.158   .151/.155   .141/.150
B       .050/.038   .050/.038   .049/.032   .054/.041   .049/.032   .226/.103
A+B     .077/.067   .077/.067   .076/.063   .080/.066   .076/.063   .098/.094
Lo      .253/.124   .191/.115   .306/.107   .297/.112   .397/.107   .278/.124
Lf      .439/.199   .261/.192   .222/.175   .420/.182   .905/.174   .645/.177
```

Per-case L; N/V identify neural-white/viscous scenes by their original `i`:

```text
Case      H      D      K     KL     NT     BK
 6      .212   .212   .192   .197   .192   .184
 7     2.881  2.882   .211   .214   .211   .201
 8      .196   .193   .225   .329   .254   .195
 9      .119   .119   .119   .120   .119   .114
10      .035   .035   .033   .034   .034   .032
11      .104   .104   .097   .107   .097   .117
36      .790   .494   .494   .534  1.668   .553
37      .316   .260   .994   .892   .344   .749
38      .106   .106   .120   .119   .120   .107
39      .082   .075   .057   .053   .069   .044
40      .142   .124   .075   .076   .085   .142
41      .085   .085   .094   .105   .094   .076
N1      .711   .484   .478  2.273   .834   .536
N3      .081   .081   .077   .096   .076   .072
N4      .209   .196   .269   .299  6.791  5.325
N5     1.551   .355   .335   .611  1.274   .351
N7      .189   .189   .178   .177   .176   .162
N8      .435   .281   .415   .426   .393   .377
V1      .130   .102   .100   .125   .154   .139
V3      .048   .048   .052   .061   .052   .048
V4      .915   .915   .172   .188   .172   .319
V5      .829   .307   .402   .538   .750   .191
V7      .116   .116   .106   .106   .106   .100
V8      .060   .058   .085   .136   .085   .122
```

Identification uses the earlier own-arm masks and pooling. Entries are **median |e| / mean e**; U is sustained-underestimate **median/worst**. `—` means no eligible samples. Original A_large is split into fast (`OF`: 36,37,39,40) and slow (`OS`: 38,41).

```text
Group Arm   Plateau       Rise          Fall        Reversal         U
A     H    .099/.041    .130/.144    .237/.286    .209/.212    .408/.694
      D    .098/.041    .130/.150    .237/.288    .209/.214    .408/.661
      K    .071/.079    .208/.220    .362/.373    .184/.248    .115/.207
      KL   .065/.074    .236/.230    .262/.344    .197/.234    .100/.195
      NT   .071/.077    .218/.231    .369/.384    .263/.288    .115/.207
      BK   .085/.108    .303/.317    .472/.506    .285/.424    .100/.192
As    H    .092/.033    .119/-.043   .235/-.124   .103/-.078   .380/.955
      D    .091/.032    .119/-.038   .235/-.124   .102/-.077   .380/.955
      K    .065/.056    .111/-.076   .083/.100        —        .091/.177
      KL   .058/.047    .125/-.076   .091/.094    .065/.057    .118/.220
      NT   .060/.054    .141/-.123   .083/.100        —        .091/.226
      BK   .113/.116    .119/-.015   .168/.163    .127/.119    .062/.146
B     H    .275/.202        —            —            —        .251/.479
      D    .275/.202        —            —            —        .251/.479
      K    .296/.240        —            —            —        .000/.326
      KL   .292/.239        —            —            —        .000/.303
      NT   .296/.240        —            —            —        .000/.326
      BK   .351/.333    .629/.565    .756/.682    .697/.634    .000/.322
A+B   H    .101/.066        —            —            —        .237/.422
      D    .101/.066        —            —            —        .237/.422
      K    .090/.088        —            —            —        .018/.132
      KL   .072/.085        —            —            —        .006/.106
      NT   .090/.088        —            —            —        .018/.132
      BK   .105/.114    .379/.383    .458/.456    .407/.413    .010/.116
OF    H    .103/.047    .142/-.116   .216/.224    .177/.122    .362/.539
      D    .128/.078    .094/-.037   .129/.103    .117/.057    .296/.560
      K    .115/.115    .151/.106    .245/.255    .234/.305    .237/.353
      KL   .128/.202    .195/.196    .326/.347    .319/.405    .175/.271
      NT   .060/.034    .160/.117    .397/.440    .471/.470    .263/.365
      BK   .099/.171    .152/.179    .280/.315    .290/.337    .198/.203
OS    H    .091/.032    .121/-.103   .112/.080        —        .270/.328
      D    .088/.031    .121/-.103   .112/.080        —        .270/.328
      K    .084/.085    .203/-.197   .315/.317    .303/.094    .204/.230
      KL   .066/.091    .144/.015    .423/.446    .510/.439    .195/.256
      NT   .084/.085    .203/-.197   .315/.317    .303/.094    .204/.230
      BK   .127/.123    .170/-.170   .359/.357    .350/.185    .191/.240
Lf    H    .123/.093    .156/-.105   .319/.289    .194/.129    .377/1.526
      D    .127/.098    .148/-.021   .188/.155    .174/.111    .317/1.526
      K    .129/.113    .183/.109    .329/.368    .334/.376    .135/.425
      KL   .136/.105    .176/.166    .327/.386    .333/.418    .156/.418
      NT   .111/.107    .238/.218    .573/.601    .608/.604    .243/1.985
      BK   .145/.166    .256/.275    .478/.540    .477/.534    .139/1.042
```

Mean whole-episode probe ripple / pooled crossing force RMS, in N:

```text
Group       H           D           K           KL          NT          BK
A       .080/.245   .080/.245   .076/.223   .077/.212   .076/.225   .070/2.338
As      .078/4.421  .078/4.421  .090/.369   .098/.303   .092/.503   .078/.281
B       .071/—      .071/—      .067/—      .066/—      .067/—      .056/.827
A+B     .083/—      .083/—      .081/—      .082/—      .081/—      .075/.202
Lo      .082/.391   .083/.315   .069/.451   .062/.404   .074/.420   .063/.425
Lf      .077/.604   .076/.380   .067/.306   .068/.605   .067/6.147  .055/2.877
```

- **K1 failed:** original mean **0.305568 > 0.195**; fresh **0.222365 ≤ 0.265**.
- **K2 met:** case 7 **0.210832 < 0.5**, versus Huber 2.881194.
- **K3 failed:** case 8 ΔL **+0.029565**; the other 22 eligible cases stayed within ±0.02.
- **K4 met:** combined large-scene mean `kf_noT − kf` = **+0.485558**.
- **K5 failed:** worst cumulative durations and their run-specific maximum ratios follow. `kf` itself passes.
- **K6 failed:** case 7 improves to **0.201266**, but elsewhere group ΔL is **+0.393, −0.00477, +0.177, +0.0226, −0.0272, +0.423**, respectively.
- **K7 failed:** pooled original plateau median |e| is **0.124228 H → 0.107370 K / 0.101563 KL**: reductions **13.57% / 18.24%**, below 15% / 35%. Sustained-underestimate medians improve in every reporting group and both original large-scene splits.

```text
K5 worst       H        D        K        KL       NT       BK
Run            1        1       37        N1       N5       14
Seconds      .384     .384     .216      .560    1.880    2.816
Max ratio   5.418    5.411   17.893    15.991    7.309    3.867
```

The main failures are **loss of the original rib-crossing benefit, especially case 37; insufficient plateau-error reduction; and prolonged over-estimates in the longer, untransported and band-fed memories**. The fresh-scene average improvement does not resolve those failures.