# 方程 → 实际代码

本文件为本次新增审阅说明，源代码为准。单位：cm、s、V、eV、A；kB=8.617333262145×10⁻⁵ eV/K，q=1.602176634×10⁻¹⁹ C。以 eV 计能量时，Vt=kB T 的数值对应伏特。

## 1. 状态、DOS 与 FD 占据

位置从左侧空穴选择接触到右侧电子选择接触。ψ=φ/Vt；u=EFn/Vt−(bp+V)/Vt；v=−EFp/Vt+bp/Vt（能量以相应 eV 数值表示）。Ec=Eg−φ，Ev=−φ。

ηn=u+ψ+(bp+V−Eg)/Vt；ηp=v−ψ−bp/Vt；sᵢ=σᵢ/(kBT)。

c(η,s)=∫ exp(−z²/2)/√(2π) · [1+exp(sz−η)]⁻¹ dz；n=N0 c(ηn,sn)，p=N0 c(ηp,sp)。

- `fd_eos.py:31`：direct_eos，801 点 Gauss–Legendre，z∈[−16,16]
- `fd_eos.py:52`：η∈[−240,0]、步长 0.02 的 C1 Hermite 插值；正 η 用粒子–空穴对称性
- `fd_eos.py:64`：GaussianDOS，导数、反函数与热力学因子
- `device.py:90`：Device.evaluate，状态转密度

数值检查范围 s≤8.7；这不是迁移率经验公式的适用范围。100 K 空穴 s≈8.5873，不能由 EOS 数值通过推论低温输运已验证。

## 2. 迁移率 μ(T,E)

μᵢ,0(T)=μᵢ,300 exp{−(4/9)(σᵢ/kB)²[1/T²−1/(300 K)²]}。

μn,300=8.4×10⁻⁴、μp,300=1.3×10⁻⁴ cm²/(V s)；σn=0.060 eV、σp=0.074 eV。100 K 分别为 4.0413529773×10⁻¹² 与 2.8954602151×10⁻¹⁷ cm²/(V s)。

- `device.py:36`：mobility；`Device.__init__` 设置 mn/mp
- `extend_device.py:21`：当前温度参数生成

Eedge=|Δψ|Vt/Δx；Ecapped=min(Eedge,235000 V/cm)，Er=1 V/cm。

Fᵢ=exp{γᵢ(T)[(Ecapped²+Er²)^(1/4)−√Er]}；μᵢ,edge=μᵢ,0 Fᵢ。

- `field_device.py:16–23`：同一正因子 Fᵢ 乘完整 SG 电流；不是只乘漂移项
- γn=0.001377 (cm/V)^(1/2)，不随温度变
- γp：225/250/275/300 K 源表点在 T 上线性插值；低于 225 K，γp=γp,225+240.33468456025972(1/T²−1/225²)
- γp,225=0.0036974093504264856；γp,100=0.02298352601266955；γp,300=0.001145
- `warm_device.py:24` 及 `reference/transport_at_target_temperatures.csv` 提供源表相对值

当前没有额外 μ(n,p) 密度因子。σ 固定；低温 μ、γ 是条件外推，未加迁移率下限。

## 3. 广义扩散与 SG 离散

局域 gᵢ=c/(dc/dη)=1/(d ln c/dη)，Dᵢ=μᵢ Vt gᵢ。

连续形式：Jn=q μn n E + q Dn ∂n/∂x；Jp=q μp p E − q Dp ∂p/∂x，E=−∂φ/∂x。

实际边上 ge=Δη/Δln(c)>0；极小 Δη 用中点导数倒数。SG 使用稳定 Bernoulli/亲和力形式，避免大正反通量相消，零准费米梯度给精确零通量。

- `device.py:40`：log_bernoulli
- `device.py:47–62`：sg_flux（实际小差阈值 |Δη|≤10⁻⁷）
- `fd_eos.py:130`：独立 secant_g 辅助接口；SG 实际在 device.py 中计算

## 4. Poisson 与连续性方程

∂²φ/∂x²=−q(p−n)/ε；∂Jn/∂x=q(R−L G)，∂Jp/∂x=−q(R−L G)。

L=0 为暗态，L=1 为此基线光强；G=1.5357919101745596×10²² cm⁻³ s⁻¹；ε=3.5 ε0。有限体积，内部控制体积 volᵢ=(Δxᵢ₋₁+Δxᵢ)/2；非均匀网格由 tanh 聚集两端，321 节点、强度 5。

- `device.py:105–120`：Device.residual，包括所有边界条件
- `device.py:121`：9 色中心差分 Jacobian
- `device.py:152`：阻尼 Newton、线搜索和数值门槛；`advance`：自适应电压/光强连续延拓
- `device.py:128`：ledger，电流守恒、Poisson 电荷、功率与闭合误差

## 5. 净复合与 β(T)

A=V/Vt+u+v=(EFn−EFp)/(kBT)，R=β(T)np[1−exp(−A)]。

`device.py:97–103` 使用保号对数实现；A=0 时 R=0；R A≥0。R 可以为负，表示这一热力学闭合下的净反向过程。只写 βnp 是大正 A 的近似；FD EOS 下也不能任意替换成 β(np−常数 ni²)。

β300=4.5633944521269105×10⁻¹² cm³/s。源表温度点按 β300×k2,author(T)/k2,author(300) 缩放；200–300 K 在 1/T 上插值 ln k2；T<200 K 用 β200 exp{−Ea/kB(1/T−1/200)}。

β200=2.3402022831420073×10⁻¹³；Ea=0.09117276441311352 eV；β100=1.1797795196788708×10⁻¹⁵ cm³/s。

- `extend_device.py:16–31`；`reference/author_fit_lines.csv`
- 这是独立经验 β 温度闭合，未自动令 β=q(μn+μp)/ε；无额外 SRH 或 CT 模块

## 6. 边界与输出

φ(0)=0，φ(d)=Eg−bn−bp−V；v(0)=0，u(d)=0；左端 Jn=0，右端 Jp=0。当前 bn=bp=0，因此相应多数载流子接触占据 c=1/2。这是理想 DOS 中心库，而非实验功函数标定。

J=平均(Jn+Jp)；光伏象限 J<0；Jsc=−1000 J(V=0)；Pout=−VJ；Voc 为带符号电流根；FF=Pmax/(Voc Jsc)（计算时保持一致单位）。

Rs=0；Jsh=0；Vterminal=Vintrinsic，Jterminal=Jintrinsic。

- `warm_device.py:45–63`：当前 observable，严格零外部分流
- `recompute_no_shunt.py`：原完整扫描、符号根和 MPP 驱动，依赖完整历史状态目录
- `review.py`：本次新增小包入口，加载最小原状态并调用同一 core

原数值门槛只约束残差、守恒与浮点相消，不覆盖网格、接触、材料参数或物理闭合误差。局域/体积分诊断采用原 finite-volume 定义；部分静态 BACE 代理使用含端点 trapezoid，不是瞬态 BACE 模拟。
