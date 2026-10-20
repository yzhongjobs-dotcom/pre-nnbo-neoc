# 参数与观测边界核查

保留两个独立工具用于核对候选共同前端及更新假设；没有完整MHE实现。

- `state_model_discriminator.py`及`state_model_discriminator/results.json`：区分已知空间状态移动、真实材料变化和估计/模型修订，说明物理变化慢不自动意味着描述符更新小。
- `contact_equivalence_audit.py`及`contact_equivalence_audit/`：在完整所用wrench模型中复核B0已有材料/表面等价退化，保留独立积分与相同历史数组。

这些是模型条件检查，不是新领域failure、主方法或闭环控制优势。参数及结果保留原样；按需运行对应脚本会覆盖同名结果，应在副本或正确Git版本中进行。主要依赖NumPy/SciPy，具体公式见源码与B0附录。
