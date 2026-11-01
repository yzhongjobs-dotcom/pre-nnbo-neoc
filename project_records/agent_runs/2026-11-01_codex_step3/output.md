已交付8个新增文件：[补丁](D:/pip_env/temp/claude/E--github-project-scan-impedance-update/03dd8b7f-853e-4b3c-8a7b-88d5bd165f74/scratchpad/step3_mpc.patch)、[压缩包](D:/pip_env/temp/claude/E--github-project-scan-impedance-update/03dd8b7f-853e-4b3c-8a7b-88d5bd165f74/scratchpad/step3_mpc_delivery.zip)、[完整报告](D:/pip_env/temp/claude/E--github-project-scan-impedance-update/03dd8b7f-853e-4b3c-8a7b-88d5bd165f74/scratchpad/step3_work/experiments/method_redesign/results/step3_mpc.md)。包含控制器、筛选与复现脚本、JSON/Markdown结果及逐拍耗时。未改原文件，未提交或推送，未使用final集。

61次调参后，选定H=10、Q/R=0.31859、终端权重=0.97764；**525/525网格点稳定，最大谱半径0.94955**。

诊断集每族50场景，RMSE中位数（N）：

|族|Proposed|MPC-U|MPC-C|MPC-oracle|Oracle|DOMPAC|
|---|---:|---:|---:|---:|---:|---:|
|S2|0.124|0.227|0.236|0.521|0.558|0.774|
|T3m|0.151|0.216|0.223|0.580|0.608|2.169|
|软组织推拉|0.077|0.080|0.081|0.057|0.019|2.121|
|硬组织推拉|0.283|0.519|2.227|0.306|0.058|6.300|
|快速推拉|0.597|2.024|3.245|0.763|0.075|7.780|

**不建议采用本次MPC-C替换现有控制层。** 五族平均RMSE均更差，配对95%区间均高于0；失控数22/250，对照2/250；长失联83/250，对照19/250。MPC-C单拍P99约1.54–1.56 ms，但15 N模型约束没有形成实际安全保证。

主要疑点是刚度上界与力偏置的可靠性、硬约束不可行后的回撤，以及线性接触预测的失配。真值输入明显改善硬组织结果，但仍远弱于现有Oracle，说明信息预测和本MPC控制设计都存在瓶颈，不能仅归因于预测。

补丁应用验证通过；固定30场景重跑的科学指标完全一致。