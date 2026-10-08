# 未采纳的密度依赖候选：仅供比较审阅

本目录与当前基线隔离。`review.py` 不加载它；它不是已经认可的新默认模型，也不是 100 K 自洽器件解。

原样摘录 3 个核心文件：
- `transport.py`：Pasveer 型密度因子、场因子，及 density_only / egdm / baseline 三种模式
- `density_device.py`：复用当前 FD/DD 核心的候选适配层，同一因子乘漂移与广义扩散一次
- `sclc.py`：用于检查归一化的单极 FD SCLC 求解器，理想 c=1/2 接触

另附已保存的 `sclc_trace_calibration.json`（各载流子/模式的全局幅值及拟合残差）和 `review_density_transport_summary.json`（历史数值检查摘要）。这些是原候选阶段资料，本次没有重新标定；全局幅值必须按载流子和模式选取，不可把 DensityDevice 默认 amplitude=1 当成历史标定值。

新增小检查的运行方式，在包根目录：

```sh
python -B optional_candidate_NOT_ADOPTED/check_candidate.py
```

只检查导入、公式极限、100 K 外推拒绝和有限边因子；不重算候选完整 J–V。未附完整 SCLC 数字化源图、标定驱动、候选全状态、有限场主方程网络，因而不能用本目录重现完整候选研究。

重要限制：
- 经验密度式的原拟合无序范围 2≤σ/(kBT)≤6；100 K 超出，默认拒绝
- c>0.1 的增强按显式约定饱和；保留 235000 V/cm 场截断，均非普适材料定律
- 有效 SCLC μ0 不自动等于 EGDM 的零密度幅值，需要前向测量模型处理
- 暖温空穴 SCLC 形状仍存在明显不匹配；数值一致性不构成物理采纳依据
- 候选减基线的 J–V 差异包含全局 SCLC 幅值归一化的影响，不能全归因于密度效应
- FD 扩散已由基线引入，不应再叠加第二个广义 Einstein 因子

公式来源标识：Pasveer 等，Physical Review Letters 94, 206601 (2005)，https://doi.org/10.1103/PhysRevLett.94.206601 。此处保留来源标识供专家查核，不表示此次打包重新完成文献审查。
