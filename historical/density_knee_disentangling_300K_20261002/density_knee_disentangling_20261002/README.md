# 缺陷密度与反偏形状的常温分离检验

先读 `REPORT.md`，再看 `figures/density_amplitude_vs_shape.png`。所有计算均为300 K、100 nm、固定FF78公共参数的现象学模型，不是样品拟合。前向FF未强制调回78%。

## 内容

- `DESIGN_BEFORE_COMPUTING.md`：主扫描前固定的范围、判据和拒绝条件；哈希见provenance.json
- `baseline_ff78.json`：完整统一配置
- `code/run_density.py`：稳定QF/循环的暗亮自洽PF、匹配局域DD及前向光伏
- `code/refine_windows.py`：161节点的21–26 V局部0.125 V细扫
- `code/analyze_density.py`：冻结DD源、每中心归一、形状和电流交点、FF
- `code/validate_density.py`：电流/电荷/源、网格、偏压步长、逐事件独立重建
- `code/check_invariants.py`：有失败即退出的验收入口
- `code/plot_density.py`：4张中文图
- `reference/`：所需原物理代码与稳定算术依赖的原样快照；压缩包只收入运行必需文件
- `data/*.npz`：全部反偏节点状态和剖面；JSON含全部端流、内部量、收敛和前向曲线
- `data/selfconsistent_curves.csv`：统一主曲线、Nt归一、局域场/占据/电荷/源分量和导数
- `data/frozen_source_reference.csv`：冻结Nt=10^15局域DD n/p/φ后评估的PF净对源，**不是新端流**
- `data/shape_landmarks.csv`、`current_crossings.csv`、`forward_metrics.csv`：描述性指标
- `data/validation.json`、各checks.csv、`final_assertions.json`：验证范围和结果
- `logs/`：实际计算运行记录；构建阶段的修正另见IMPLEMENTATION_NOTES.md

## 重跑

已实际检查Python 3.12.14、NumPy 2.3.5、SciPy 1.17.0、Matplotlib 3.10.8。无需访问网络或原项目；在解压目录运行：

```bash
bash reproduce.sh
```

该入口会重写本目录的生成数据、图和日志，复现前可保留原包。不触碰原项目。中文图使用系统Noto Sans CJK字体；数值计算不依赖该字体。内存与计算量适合单机；如需加速，可自行并行独立run_density命令。

快速检查既有数据，不重新计算：

```bash
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python -B code/check_invariants.py
```

单组例子：

```bash
python -B code/run_density.py PF 1e15 0 81 .25 30
python -B code/run_density.py DD 1e15 1 81 .25 30
```

位置单位cm，电流A/cm²，能量eV，温度K，Nt为cm⁻³。CSV列名明确标注单位，绘图再换成nm、mA/cm²或MV/cm。NPZ电势phi为V，n/p为cm⁻³，field为V/cm，pair_rate_s为每中心每秒；Rn/Rp为cm⁻³ s⁻¹，bim同单位，Jn/Jp为A/cm²。

## 关键边界

反偏幅值U=−V；正提取电流I=−J。U25/U50/U100是人为电流交点，未跨越应为空。本轮都跨越。曲率最大点不自动是启动电压；严格唯一膝点判据本轮失败，较小局部肩部也有窗口依赖。高场平台来自既定势垒构造；等温大电流不表示真实器件能够持续承受。
