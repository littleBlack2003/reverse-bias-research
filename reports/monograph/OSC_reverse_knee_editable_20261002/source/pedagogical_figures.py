"""Eight authored conceptual TikZ reading aids. No numerical data or model runs."""
from pathlib import Path
import json
P=Path(__file__).resolve().parents[1]
D={};A=[]
START=r'''\begin{tikzpicture}[x=1cm,y=1cm,font=\small,
box/.style={draw=tolblue,line width=.75pt,fill=tolblue!3,align=center,inner sep=7pt},
arr/.style={-{Stealth},line width=.85pt,draw=tolblue},
rev/.style={-{Stealth},line width=.75pt,draw=tolorange},
both/.style={<->,>=Stealth,line width=.8pt,draw=tolblue}]
'''
def add(section,name,body,caption):
 D[name]=r'\resizebox{.90\linewidth}{!}{'+START+body+r'\end{tikzpicture}}'
 A.append({'section':section,'name':name,'caption':caption,'provenance':'Original conceptual schematic derived from the adjacent existing text and equations; no new model result or measured parameter.'})
add(10,'pair_contact_bookkeeping',r'''
\node[font=\bfseries]at(6,3.8){{\cjkfont 箭头表示粒子运动；先问电子从哪里来}};
\draw[fill=gray!12,draw=gray](0,0)rectangle(1.7,3.15);
\draw[fill=gray!12,draw=gray](10.3,0)rectangle(12,3.15);
\node[align=center]at(.85,1.7){{\cjkfont 左电极}\\{\cjkfont 空穴侧}};
\node[align=center]at(11.15,1.7){{\cjkfont 右电极}\\{\cjkfont 电子侧}};
\draw[dashed,gray](1.7,0)rectangle(10.3,3.15);
\node at(6,2.75){{\cjkfont 活性层：体内成对出现} $e^-+h^+$};
\node[draw=tolorange,circle,inner sep=3pt] (h) at(5.05,1.7){$h^+$};
\node[draw=tolblue,circle,inner sep=3pt] (e) at(6.95,1.7){$e^-$};
\draw[rev](h)--node[above]{{\cjkfont 空穴提取}}(1.9,1.7);
\draw[arr](e)--node[above]{{\cjkfont 电子提取}}(10.1,1.7);
\node at(6,.65){{\cjkfont 体源} $G-R$ {\cjkfont 同时进入两条连续性方程}};
\draw[arr](11.1,-.65)--node[below]{{\cjkfont 电极注入已有电子：边界输入}}(7.2,-.65);
\draw[rev](7.2,-1.4)--node[below]{{\cjkfont 同一边界也允许逆向提取}}(11.1,-1.4);
\draw[arr,black](1.1,-2.25)--node[below]{$x$ {\cjkfont 正向}}(10.9,-2.25);
\node[align=center]at(3.5,-1.10){{\cjkfont 产生与注入不可重复计数}\\{\cjkfont 粒子运动不等于常规电流方向}};
''','怎样读粒子账本：上半图是假设由体内产对并向两端收集的情形，箭头表示电子或空穴运动，不是常规电流方向。下半图单独画出右接触的电子注入与提取；它们是边界交换，不能再当成新增体产对。实际两类过程可共存，需用四个接触净通量和体源积分区分。电极标签只沿用本书坐标约定，不代表实际选择性已被测定。')
add(18,'pf_schottky_origins',r'''
\node[font=\bfseries]at(2.5,3){{\cjkfont 体内PF：真实带电中心}};
\node[font=\bfseries]at(9.0,3){{\cjkfont Schottky：平面像电荷}};
\node[draw=tolorange,circle,minimum size=.70cm]at(1.1,1.65){$+$};
\node[draw=tolblue,circle,minimum size=.52cm]at(3.65,1.65){$-$};
\draw[rev](3.25,1.35)--node[below]{{\cjkfont 库仑吸引}}(1.55,1.35);
\draw[arr](3.95,1.65)--node[above]{{\cjkfont 外场力}}(5.15,1.65);
\draw[both,gray](1.1,.15)--node[below]{$r$}(3.65,.15);
\node at(2.8,-.8){$U_{\rm C}(r)=-\dfrac{q^2}{4\pi\varepsilon r}$};
\fill[gray!15](5.8,.4)rectangle(8.0,2.5);
\draw[gray,line width=1pt](8,.35)--(8,2.55);
\node at(6.7,2.65){{\cjkfont 金属}};
\node[draw=tolorange,dashed,circle,minimum size=.65cm]at(6.7,1.65){$+$};
\node[draw=tolblue,circle,minimum size=.52cm]at(9.3,1.65){$-$};
\node[font=\footnotesize]at(6.7,1.1){{\cjkfont 等效像电荷}};
\draw[rev](9.0,.75)--(8.2,.75);
\node at(9.6,.42){{\cjkfont 像力吸引}};
\draw[arr](9.65,1.65)--node[above]{{\cjkfont 外场力}}(11.9,1.65);
\draw[both,gray](8,.05)--node[below]{$x$}(9.3,.05);
\node at(9.0,-.8){$U_{\rm im}(x)=-\dfrac{q^2}{16\pi\varepsilon x}$};
\node at(6,-1.9){{\cjkfont 相同理想介电常数与场幅值下：}\quad $\Delta_{\rm PF}=2\Delta_{\rm S}$};
''','两种降垒先从势的来源区分。左图的正中心吸引真实电子；右图用金属平面内的等效正像电荷表示像力，像电荷不是独立供给中心。蓝箭头表示使电子离开的外场力，橙箭头表示吸引力。两条势能式采用SI单位；像势含金属极化的自能因子，不能直接当作两独立电荷的相互作用能。二倍降垒关系只属于这些理想势及共同介电常数、场幅值，不直接标定真实有机界面。')
add(19,'balance_rates_cycles',r'''
\node[font=\bfseries]at(2.75,3.4){{\cjkfont 一条边：速率与通量分开}};
\node[box,minimum width=1.35cm] (i) at(.9,1.85){$i$\\$G_i,P_i$};
\node[box,minimum width=1.35cm] (j) at(4.6,1.85){$j$\\$G_j,P_j$};
\draw[arr](1.7,2.05)--node[above]{$k_{ij}$}(3.8,2.05);
\draw[rev](3.8,1.60)--node[below]{$k_{ji}$}(1.7,1.60);
\node at(2.75,.55){$k_{ij}/k_{ji}=e^{-\Delta G_{ij}/k_BT}$};
\node at(2.75,-.12){{\cjkfont 平衡：} $P_i k_{ij}=P_j k_{ji}$};
\node[align=center]at(2.75,-1.0){{\cjkfont 两个率共同乘以} $g>0$\\{\cjkfont 比值不变，等待时间改变}};
\node[font=\bfseries]at(9.2,3.4){{\cjkfont 闭合循环：逐边定义一致}};
\node[box,circle,minimum size=.65cm] (a)at(9.2,2.45){$i$};
\node[box,circle,minimum size=.65cm] (b)at(7.45,.8){$j$};
\node[box,circle,minimum size=.65cm] (c)at(10.95,.8){$k$};
\draw[both](a)--(b);\draw[both](b)--(c);\draw[both](c)--(a);
\node[align=center]at(9.2,-.2){{\cjkfont 同温、共同电化学势的平衡}\\$\mathcal A_{\rm cycle}=0$};
\node[align=center]at(9.2,-1.15){{\cjkfont 静电内建场可非零}\\{\cjkfont 仍不能驱动持续净循环}};
\node[draw=gray,dashed,inner sep=6pt]at(6,-2.2){{\cjkfont 持续循环需非平衡驱动，例如}\ $\mu_L\ne\mu_R$\quad {\cjkfont 驱动力和局部场功不能重复计入}};
''','左图箭头表示允许的跃迁，不表示已有净流；Pi、Pj为状态概率。详细平衡比较的是概率乘速率的通量，而不是要求两方向速率相等。右图是抽象的可逆状态图，节点不指定材料或电荷态；在同一热浴与共同电化学势下，闭合路径的状态能变化相消。不同储库的驱动必须进入同一能量与粒子定义，不能把非零静电场本身当成另一台泵。')
add(22,'finite_pf_endpoints',r'''
\node[font=\bfseries]at(6,3.8){{\cjkfont 只画最有利的一对端点；实线为电子转移}};
\node[box,minimum width=2.2cm] (v)at(1,2.1){{\cjkfont 价态端点}\\$x_h=x-\ell$};
\node[box,draw=tolorange,minimum width=2.2cm] (t)at(6,2.1){{\cjkfont 中心}\\$x$};
\node[box,draw=tolgreen,minimum width=2.2cm] (c)at(11,2.1){{\cjkfont 传导端点}\\$x_e=x+\ell$};
\draw[arr](2.15,2.4)--node[above]{{\cjkfont 空穴发射}}(4.85,2.4);
\draw[arr](7.15,2.4)--node[above]{{\cjkfont 电子发射}}(9.85,2.4);
\draw[rev,dashed](4.85,1.8)--node[below]{{\cjkfont 空穴捕获}}(2.15,1.8);
\draw[rev,dashed](9.85,1.8)--node[below]{{\cjkfont 电子捕获}}(7.15,1.8);
\node at(1,.85){{\cjkfont 最后留下} $h^+$};
\node at(6,.85){{\cjkfont 中心复位}};
\node at(11,.85){{\cjkfont 最后留下} $e^-$};
\draw[both,gray](1,.05)--node[above]{$\ell$}(6,.05);
\draw[both,gray](6,.05)--node[above]{$\ell$}(11,.05);
\draw[arr,black](1,-.6)--node[below]{$+x$ {\cjkfont 方向；图示有利情形下} $\varphi$ {\cjkfont 升高}}(11,-.6);
\draw[rev](9,-1.65)--node[above]{{\cjkfont 电场}\ $E=-\partial_x\varphi$}(3,-1.65);
\node[box,minimum width=5.6cm]at(2.9,-2.9){{\cjkfont 端点能量差（SI）}\\$E_g-q[\varphi(x_e)-\varphi(x_h)]$};
\node[box,draw=tolorange,minimum width=5.6cm]at(9.1,-2.9){{\cjkfont PF项改变过渡态势垒}\\{\cjkfont 每条边的正反率共享能量定义}};
''','从左到右读一次完整产对：价态电子进入中心后，中心电子到传导端点，最终在两侧留下同一对空穴和电子。虚线是必须保留的逆过程；图只抽出原模型最有利的非局域端点，未画其它方向和局域通道。图中横轴是真实模型位置而非反应坐标，左向电场有利于电子向右、空穴向左分离。下方把端点电势做功与过渡态降垒分开，场功只计一次；SI能量式显式含q。距离和电子可达性仍是模型假设，不由示意图证明。')
add(25,'power_heat_boundaries',r'''
\node[box,minimum width=2.7cm,minimum height=1.3cm](port)at(.9,1.5){{\cjkfont 外部端口}\\$P_{\rm ext}=VJ$};
\node[box,minimum width=4.8cm,minimum height=1.8cm](dev)at(6,1.5){{\cjkfont 器件内部场功}\\$P_{\rm int}=\int J_{\rm tot}E\,dx$\\{\cjkfont 非局域贡献仅是其中一部分}};
\node[box,draw=gray,minimum width=2.9cm,minimum height=1.3cm](bath)at(11.2,1.5){{\cjkfont 热浴与其它能流}\\{\cjkfont 需另行闭合}};
\draw[arr](port)--(dev);
\draw[both,dashed,gray](dev)--(bath);
\node[box,draw=gray,minimum width=6.1cm](contact)at(6,3.8){{\cjkfont 接触的化学与能量边界项}};
\draw[both,dashed,gray](contact)--(dev);
\node at(6,-.25){$P_{\rm int}=(V-V_{\rm bi})J\,,\qquad P_{\rm int}-P_{\rm ext}=-V_{\rm bi}J$};
\node[align=center]at(6,-1.3){{\cjkfont 内建势包含在内部静电势差中}\\{\cjkfont 只有完整能量账本闭合后，才可讨论热流与温升}};
''','这张图画能量账本的边界，而不是已计算出的热流分配。实线连接外部电功与器件内部描述；虚线提示还需要明确的接触、热浴和其它能流。两行功率关系使用本书电压、电流符号，均为单位面积的功率。内部场功与外部端口功的差不能直接当成热，非局域场功也不是全部器件功率；原模型没有据此求得温升。')
add(32,'wkb_path_count',r'''
\node[font=\bfseries]at(6,4.0){{\cjkfont 上：固定电子能量下的顺序微结（能量高度任意）}};
\draw[fill=tolblue!8,draw=tolblue](0,0)rectangle(1.8,2.8);
\draw[fill=tolgreen!8,draw=tolgreen](10.2,0)rectangle(12,2.8);
\node[align=center]at(.9,2.2){{\cjkfont 左有效库}\\$f_L(E)$};
\node[align=center]at(11.1,2.2){{\cjkfont 右有效库}\\$f_R(E)$};
\draw[draw=gray,line width=.9pt](2.1,.3)--(3.1,2.6)--(4.8,.3);
\draw[draw=gray,line width=.9pt](7.2,.3)--(8.2,2.6)--(9.9,.3);
\draw[dashed,gray](0,1.05)--(12,1.05);
\node[anchor=east]at(-.08,1.05){$E$};
\draw[tolorange,line width=1.3pt](5.45,1.05)--(6.55,1.05);
\node[align=center]at(6,.1){{\cjkfont 陷阱}\ $f_t$};
\draw[both](1.9,3.15)--node[above]{$\Gamma_L=\nu e^{-S_L}$}(5.35,3.15);
\draw[both](6.65,3.15)--node[above]{$\Gamma_R=\nu e^{-S_R}$}(10.1,3.15);
\node at(6,-.85){{\cjkfont 两库交换同一电子；不是额外的体产对}};
\draw[gray,line width=.85pt](.6,-2.1)--(.6,-3.6);
\draw[gray,line width=.85pt](11.4,-2.1)--(11.4,-3.6);
\foreach \yy in {-2.45,-3.15}{
 \draw[tolblue,line width=.8pt](.6,\yy)--(11.4,\yy);
 \foreach \xx in {3,6,9}{\node[circle,fill=white,draw=tolblue,minimum size=.23cm,inner sep=0pt]at(\xx,\yy){};}
}
\node at(6,-1.70){{\cjkfont 下：完整并联路径的几何计数示意}};
\node at(6,-4.05){{$N_{\rm path}$ {\cjkfont 数贯通通道的面密度，不数每个串联微结}}};
''','上半图仅说明固定能量E的电子通过两侧空间势垒，与一个陷阱顺序交换；Γ描述隧穿交换强度，实际进出还受储库费米占据与陷阱空满限制。折线势垒是概念草图，不是重新计算或实测势形，也不表示左右事件相干同时发生。下半图的两条线仅作计数例子：每条贯通路径可以含多个串联微结，不能把全部微结都当成独立并联通道。因此体陷阱数乘厚度不自动等于Npath。含多个串联微结的路径还需联合求解占据与交换，不能直接复用上方单微结的净率。')
add(49,'escape_collection',r'''
\node[font=\bfseries]at(6,3.7){{\cjkfont 局部离开成功，与到达电极，是两层问题}};
\node[box,minimum width=2.5cm](origin)at(1,2.0){{\cjkfont 返回母态}\\{\cjkfont 或原束缚态}};
\node[box,draw=tolorange,minimum width=2.5cm](near)at(5.7,2.0){{\cjkfont 刚形成的近态}\\{\cjkfont 从这里开始竞争}};
\node[box,draw=tolgreen,minimum width=2.5cm](away)at(10.7,2.0){{\cjkfont 继续分离}\\{\cjkfont 尚未到电极}};
\draw[rev](near)--node[above]{$b$ {\cjkfont 回跳}}(origin);
\draw[arr](near)--node[above]{$u$ {\cjkfont 离开}}(away);
\node at(5.7,.55){$P_{\rm escape}=\dfrac{u}{u+b}$};
\draw[gray,dashed](.1,-.35)rectangle(11.9,-2.35);
\node at(6,-.68){{\cjkfont 后续还要经过器件输运与收集}};
\node[draw=tolblue,circle,inner sep=3pt](e)at(1,-1.45){$e^-$};
\foreach \xx in {3,5,7,9}{\draw[tolgreen,line width=1pt](\xx-.25,-1.45)--(\xx+.25,-1.45);}
\draw[arr](1.5,-1.45)--(2.6,-1.45);\draw[arr](3.4,-1.45)--(4.6,-1.45);\draw[arr](5.4,-1.45)--(6.6,-1.45);\draw[arr](7.4,-1.45)--(8.6,-1.45);\draw[arr](9.4,-1.45)--(10.4,-1.45);
\node[box,draw=gray,minimum width=.9cm,minimum height=.8cm]at(11,-1.45){{\cjkfont 电极}};
\node at(6,-3.05){{\cjkfont 输运时间、再次复合、接触与空间电荷都可能改变收集}};
''','上排从已经形成的近态开始，只在两个恒定、互斥的指数等待出口u和b这一简化条件下，先成功离开的概率为u/(u+b)。箭头表示状态转移；离开这个近态不等于已被电极提取，也不排除更远处继续回跳。下排短线与箭头表示后续输运示意，没有假定固定跳数、能量阶梯或单向真实网络。最终端口电流仍需与占据、复合和自洽输运一起计算。')
add(57,'mechanism_diagnostics',r'''
\node[box,minimum width=2.8cm,minimum height=1.55cm](a)at(0,2.7){{\cjkfont 体产对}\\{\cjkfont 形成、回填、分离}};
\node[box,draw=tolgreen,minimum width=2.8cm,minimum height=1.55cm](b)at(3.35,2.7){{\cjkfont 接触供给}\\{\cjkfont 注入、传输、提取}};
\node[box,draw=tolorange,minimum width=2.8cm,minimum height=1.55cm](c)at(6.7,2.7){{\cjkfont 电热反馈}\\{\cjkfont 功率、温升、速率}};
\node[box,draw=violet,minimum width=2.8cm,minimum height=1.55cm](d)at(10.05,2.7){{\cjkfont 空间弱区}\\{\cjkfont 局部面积、局部场}};
\node[box,minimum width=12.9cm,minimum height=1.0cm](obs)at(5.025,.35){{\cjkfont 可以产生相似的电流上翘，且机制可以共存}\\{\cjkfont 单一阈值或一段直线不足以唯一鉴定}};
\draw[arr](a.south)--(0,.95);\draw[arr](b.south)--(3.35,.95);\draw[arr](c.south)--(6.7,.95);\draw[arr](d.south)--(10.05,.95);
\node[box,draw=gray,minimum width=12.9cm,minimum height=1.2cm](tests)at(5.025,-1.75){{\cjkfont 联合观察与共同约束}\\{\cjkfont 温度＋输运时间}\quad {\cjkfont 接触＋暗亮幅值}\quad {\cjkfont 时长＋散热}\quad {\cjkfont 面积＋空间分布}};
\draw[both](obs)--(tests);
\node at(5.025,-3.0){{\cjkfont 同像素、实际端压、相同历史；重复循环区分可逆变化与持久损伤}};
''','先把四类候选看成可共存的原因，再用下方多种观察共同约束。连线表示“需要检验的可能关系”，不是已识别的样品机制，也没有声称某一种观测是唯一指纹。例如低温电流变小可以来自产生或提取，时长变化也须配合散热与初态。图的作用是帮助安排理论问题和未来判别逻辑，不是新增实验结论或测量操作方案。')
(P/'source/pedagogical_figures.json').write_text(json.dumps(D,ensure_ascii=False,indent=2))
(P/'source/reading_aids.json').write_text(json.dumps(A,ensure_ascii=False,indent=2))
print('Authored',len(A),'conceptual reading aids; no numerical data generated')
