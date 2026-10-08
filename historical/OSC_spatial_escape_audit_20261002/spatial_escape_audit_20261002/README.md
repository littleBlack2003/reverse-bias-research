# 空间逃逸与相关截断：可复现入口

本目录是2026-10-02新专题，旧量子循环/储库/温度检查点未更改。先读REPORT.md，六页概要见output/中的PDF；完整先行能量定义见THEORY_BEFORE_CODE.md。

## 环境与复算

依赖Python、NumPy、SciPy、mpmath、Matplotlib；PDF额外用reportlab与Pillow。实际版本见results/provenance.json。以下脚本只覆写本专题结果目录，不修改输入副本和上游模型。

```bash
export OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1
export MPLCONFIGDIR=/tmp/mpl-spatial XDG_CACHE_HOME=/tmp/spatial-cache
python -B run_audit.py
python -B run_followup_controls.py
python -B plot_results.py
# 可选：重建六页概要
python -B create_brief_pdf.py
```

默认250位状态概率求解，10点用350位交叉；283个表列控制。原核仍双精度，因此状态账本极小误差并不提高核或材料参数的物理精度。代表运行约35秒，不需要大网格。

## 文件

- spatial_model.py：统一状态能、图、谱交换、物理账本；没有裁负概率
- run_audit.py：197主控制、规范/精度/积分、Schur与强可聚合测试
- run_followup_controls.py：86个多电荷/核/受限平衡控制，高场详细账本及瞬时/稳态比较
- results/：全部CSV、JSON、核验与来源；CSV科学计数法中的极小数是形式数学均值
- inputs/：本次使用的上游只读源码副本，原件哈希见SOURCE_LEDGER.json
- review/：独立物理审阅、冻结源码及独立小点核验
- figures/：4张静态图；连线不是实验拟合或已确定阈值
- output/：六页中文概要；完整方程/限定与文献见REPORT.md
- run_initial_130digit_failure.log：保留首次僵硬点概率失败，不裁剪后重跑

## 三条不可混淆的限定

1. 旧colocated是10nm参考；direct与explicit才共享12nm外部窗口
2. last-tag扫描改变相关截断/几何，不能叫数值网格收敛；distance衰减控制不保持相同实际总H²
3. single_pair参数是附加多粒子禁占诊断，改变物理状态空间；不是原32态模型的精确粗粒化

本阶段无Poisson/DD闭合、无器件J–V、FF、V50或击穿预测，无目标材料参数拟合。
