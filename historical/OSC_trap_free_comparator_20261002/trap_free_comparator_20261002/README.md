# 无缺陷跨隙对照

当前常温300 K入口：`room_temperature_300K/REPORT_300K.md`，图见该目录`figures/room_temperature_300K.png`。

原`REPORT.md`与数据保留为前一阶段记录；其中80 K例子不是当前主基准。共同物理参数仍见`parameters.json`。

常温补充运行：`python room_temperature_300K/audit_300K.py`、`python room_temperature_300K/plot_300K.py`。

运行：
```bash
cd /workspace/shared/wxh/research/trap_free_comparator_20261002
python run_audit.py
python derive_ledgers.py
python final_precision_checks.py
python plot_results.py
```

依赖Python、NumPy、SciPy、mpmath、Matplotlib。输入中保存了旧核和储库实现，代码可独立于原专题重现。

- `comparator_model.py`：两轨道四态与三轨道八态模型
- `results/fixed_budget_comparison.csv`：相同内部耦合预算
- `results/distance_law_comparison.csv`：同一距离律、不重归一预算
- `results/representative_edge_ledger.csv`：边、电荷、能量、电子/声子热账本
- `results/nonzero_field_equilibrium.csv`：非零场平衡
- `results/common_bath_construction.json`：共同浴位移结构
- `results/validation.json`：验证指标
- `sources/`：原文及核验记录

无缺陷模型的密度名为Npair。图率、条件电流、实测器件电流不得混同。
