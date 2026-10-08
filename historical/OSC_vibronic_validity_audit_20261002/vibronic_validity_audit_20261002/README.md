# 振动核物理适用性专题

先读 `REPORT.md`。这是对已交付理论的独立补充，不覆盖旧报告、原始实验、FF78 配置或任何器件拟合。

运行：

```bash
cd /workspace/shared/wxh/research/vibronic_validity_audit_20261002
export MPLBACKEND=Agg MPLCONFIGDIR=/tmp/vibronic-audit-mpl XDG_CACHE_HOME=/tmp/vibronic-audit-cache
python kernel_benchmark.py
python focused_checks.py
python plot_results.py
```

主要输出：
- `results/fully_quantum_continuum_comparison.csv`：同总重组能的全量子低频浴对照
- `results/matched_lambda_high_mode_comparison.csv`：高频一模/两模对照
- `results/classical_variance_thresholds.csv`：经典方差误差门槛
- `results/weak_coupling_scales.csv`：弱耦合和典型 LZ 尺度
- `results/prior_graph_thermalization_scales.csv`：旧图总离开率与热化诊断
- `results/representative_edge_kernel_sensitivity.csv`：旧图实际单边能量上的核敏感性
- `results/validation.json` 和 `focused_checks.json`：数值验证
- `sources_review/`：原文、MinerU OCR、原页核验与独立解析审查

全部新谱参数均是条件场景；不代表材料实测。未重跑完整产生图或器件 J–V，也没有修改55页原PDF。
