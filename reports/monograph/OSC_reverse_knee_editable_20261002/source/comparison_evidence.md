# 晶硅、有机与卤化物钙钛矿太阳能电池的反偏行为：限定条件下的比较

**读表约定：** 电流开始上翘（膝点）、强反向导电与永久损伤不是同一事件；表内文献实例不是普适阈值。温度比较使用反向电压绝对值 |V|，并须固定同一个判据、扫描程序和器件状态。

| 比较项 | 晶硅 | 有机（以无序体异质结为主） | 卤化物钙钛矿 |
|---|---|---|---|
| 主要机制 | 雪崩倍增；高掺杂窄结可有带间隧穿；缺陷、结边缘与局部高场可导致提前导通 [1,2] | 候选包括接触注入/隧穿、陷阱辅助或其他场增强产生，以及局部弱区导电；不能预先把本项目膝点归为 PF 或永久丝状短路。已有特定体系报告深陷阱关联的反向隧穿及可逆/不可逆阶段 [3] | 离子再分布可改变界面场并促进隧穿；注入载流子与卤素/电极的电化学反应、局部缺陷及热效应可共同参与 [5,6] |
| 输运特点 | 电子/空穴以能带输运为主，高场下可发生碰撞电离；不能只以“自由程长”概括雪崩条件 [1,2] | 常以局域态间跳跃描述，并受无序、占据及形貌连通性控制；跳跃距离不能直接等同于能带载流子的自由程 [4] | 电子/空穴输运与可移动离子并存；离子通常通过改变场、势垒和化学状态影响电子电流 [5,6] |
| 触发场和电压 | 随掺杂、耗尽层、结曲率和缺陷而变；局部提前导通与理想均匀结击穿必须区分，不给统一临界场或电池电压 [1,2] | 本项目若将 10 V 全部分配于 100 nm 活性层，名义平均场约 1 MV/cm；实际局部场需自洽求解。该换算不是材料普适阈值；特定优化 OSC 已报道不可逆击穿 |V|>35 V [3] | 早期特定器件报告 −1 至 −4 V [5]；优化阻挡层后已有 |Vbd|>20 V 的实例 [6]。端电压不能直接代表界面局部场 |
| 温度依赖 | 理想机制下，雪崩主导的 |Vbd| 通常随温度升高而增大，Zener 主导者可减小；混合机制及缺陷导电须另判 [1] | 热激活发射/注入与隧穿的温度响应不同；固定电流阈值与形状膝点也可有不同温度趋势。本项目尚无变温定论，不能预设统一正/负温度系数 | 离子运动、反应动力学与电子输运共同受温度影响；表观阈值还依赖测试时长和预处理，不能用单一系数概括 [5,6] |
| 时间和历史依赖 | 纯电子响应可较快，但自热、缺陷演化与应力损伤会引入时间/历史依赖；不能统称“弱” [2] | 本项目观察：在已测试范围，扫描快慢及是否采用脉冲对膝点影响不明显；具体脉宽等条件仍待补齐。这不是所有有机器件的通则 | 常见明显的偏压历史、扫描程序和保持时间影响，涉及离子再分布与反应；程度依结构与条件而变 [5,6] |
| 损伤与可逆性 | 局部电流集中可形成热点并导致永久损伤；进入反向导电本身不等于已经损坏 [2] | 本项目轻微跨过膝点后，正向 J–V 基本不变，再次反扫可重现；更强应力下才需另查永久短路、热损伤或电极/材料变化。不可把未验证的银迁移写成本器件既定机制 | 可有可恢复性能变化，也可有永久局部短路、金属丝、电极腐蚀及材料分解；不保证所有器件都依次经历“先软后硬” [5,6] |
| 常见缓解方向 | 模块旁路保护；控制缺陷与局部场集中，改善散热，并验证反偏耐受性 | 控制膜厚/形貌均匀性和缺陷，优化接触选择性与阻挡层，并采用限流/保护；具体方案须与主导机制匹配 [3] | 使用离子/金属阻挡层、更稳定的接触及传输层，降低注入与局部缺陷电流，并配合模块保护 [6] |

## 文献与证据边界

[1] M. Singh Tyagi, “Zener and avalanche breakdown in silicon alloyed p–n junctions—II: Effect of temperature on the reverse characteristics and criteria for distinguishing between the two breakdown mechanisms”, Solid-State Electronics 11, 117–128 (1968). DOI: https://doi.org/10.1016/0038-1101(68)90142-1 。本次核对出版社原始论文摘要，明确区分 Zener 与 avalanche 的电压温度系数。使用绝对值避开负电压符号歧义。不要把反向电流温度系数与击穿电压温度系数混用。

[2] “Influence of surface texture on the defect-induced breakdown behavior of multicrystalline silicon solar cells”, Progress in Photovoltaics (2013). DOI: https://doi.org/10.1002/pip.1226 。原始研究摘要确认纹理刻蚀/晶界与位错蚀坑影响局部击穿。补充原始研究：“Thermal and electrical investigation of the reverse bias degradation of silicon solar cells”, Microelectronics Reliability (2013), https://www.sciencedirect.com/science/article/abs/pii/S0026271413001868 ，确认反偏应力的热、电和永久损伤演变。

[3] J. Huang et al., “Perovskite–organic tandem solar cells with superior reverse-bias stability”, Nature Materials 25, 1419–1428 (2026). DOI: https://doi.org/10.1038/s41563-026-02541-6 。本次出版社正文抓取失败，已核对作者机构正式书目及原始摘要：https://research.polyu.edu.hk/en/publications/perovskiteorganic-tandem-solar-cells-with-superior-reverse-bias-s/ 。只引用摘要明确支持的深陷阱关联反向隧穿、可逆/不可逆阶段和优化 OSC 不可逆击穿绝对值超过 35 V；不将其泛化到本项目 D18:L8-BO。

[4] W. F. Pasveer et al., “Unified Description of Charge-Carrier Mobilities in Disordered Semiconducting Polymers”, Physical Review Letters 94, 206601 (2005). DOI: https://doi.org/10.1103/PhysRevLett.94.206601 。已核对 APS 原始摘要与作者机构公开全文 https://pure.tue.nl/ws/files/2107419/Metis188708.pdf ，支持无序局域态跳跃及温度/场/密度依赖，不支持“所有有机半导体自由程均为一个格点”。

[5] A. R. Bowring, L. Bertoluzzi, B. C. O'Regan, M. D. McGehee, “Reverse Bias Behavior of Halide Perovskite Solar Cells”, Advanced Energy Materials, 1702365 (2017). DOI: https://doi.org/10.1002/aenm.201702365 。已读作者公开全文：https://web.stanford.edu/group/mcgehee/publications/AEM2017.pdf 。第 2 页所研究多种接触器件为 −1 至 −4 V；论述部分恢复、离子致界面带弯曲/隧穿及局部短路差异。此范围不能扩展为所有钙钛矿。

[6] “Barrier reinforcement for enhanced perovskite solar cell stability under reverse bias”, Nature Energy 9, 1264–1274 (2024). DOI: https://doi.org/10.1038/s41560-024-01579-7 。已核对出版社原始研究摘要，支持碘/铜氧化还原、铜迁移与金属丝，以及 LiF/SnO2/ITO 阻挡设计后 |Vbd|>20 V、−1.6 V 下约 1000 h 的 T90。表格不需要纳入寿命数字。

## 修改说明（供编辑）

删除未严格限定的 0.3/1 MV/cm、硅“几 V 到二三十 V”、钙钛矿 0.02–0.1 MV/cm，避免混淆局部峰值场、平均场及材料临界场。保留 100 nm / 10 V 的纯单位换算。温度一行将有机的“击穿电压随温度降低”改成待判机制依赖，不能由热激活电流直接推导膝点温度系数。银迁移在本项目没有证据，不作为既定损伤类型。表内本项目实验描述来源于用户陈述，未核对原始曲线与脉冲波形。
