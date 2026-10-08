"""Native TikZ mechanisms and pgfplots curves from frozen, hashed CSV only."""
from pathlib import Path
import csv,json,hashlib,shutil
P=Path(__file__).resolve().parents[1]
D={}
start=r'''\begin{tikzpicture}[x=1cm,y=1cm,box/.style={draw=tolblue,fill=tolblue!4,line width=.8pt,align=center,inner sep=8pt,font=\small},bi/.style={<->,>=Stealth,line width=.85pt,draw=tolblue}]
'''
def keep(name,s):D[name]=r'\resizebox{.90\linewidth}{!}{'+start+s+r'\end{tikzpicture}}'
keep('bands_localized',r'''
% Conceptual electronic energy layout: no measured organic energy or gap.
\node[font=\small\bfseries]at(2.35,3.3){(a) {\cjkfont 周期无机晶体}};
\node[font=\small\bfseries]at(9.15,3.3){(b) {\cjkfont 典型无序有机半导体}};
\draw[-{Stealth},line width=.65pt](-.55,-.05)--(-.55,2.85);
\node[rotate=90,font=\small]at(-.95,1.4){{\cjkfont 电子能量}};
\draw[fill=tolblue!12,draw=tolblue,line width=.8pt](0,2.1)rectangle(4.6,2.75);
\draw[fill=tolblue!12,draw=tolblue,line width=.8pt](0,0)rectangle(4.6,.65);
\node[font=\small]at(2.3,2.43){{\cjkfont 导带：延展态}};
\node[font=\small]at(2.3,.33){{\cjkfont 价带：延展态}};
\node[anchor=west,font=\small]at(4.68,2.1){$E_C$};
\node[anchor=west,font=\small]at(4.68,.65){$E_V$};
\draw[<->,>=Stealth,line width=.75pt](2.0,.70)--node[right,font=\small]{$E_g=E_C-E_V$}(2.0,2.05);
\node[font=\small]at(2.3,-.6){{\cjkfont 周期势}\quad Bloch{\cjkfont 态与能带色散}};
\draw[-{Stealth},line width=.65pt](6.25,-.05)--(6.25,2.85);
% Shaded envelopes mark disorder spread, not a continuous Bloch band.
\fill[tolgreen!7](6.65,-.05)rectangle(11.65,.70);
\fill[tolgreen!7](6.65,2.02)rectangle(11.65,2.78);
\foreach \x/\y in {6.85/2.20,7.75/2.57,8.65/2.10,9.55/2.69,10.45/2.35}{
 \draw[draw=tolgreen!80!black,line width=1.2pt](\x,\y)--++(.72,0);
}
\foreach \x/\y in {6.85/.20,7.75/.52,8.65/.08,9.55/.39,10.45/.61}{
 \draw[draw=tolgreen!80!black,line width=1.2pt](\x,\y)--++(.72,0);
}
\node[font=\small]at(9.15,2.96){LUMO{\cjkfont 衍生态}};
\node[font=\small]at(9.15,-.25){HOMO{\cjkfont 衍生态}};
\node[font=\small,align=center]at(9.15,1.40){{\cjkfont 能量与空间分布依赖局部环境}\\{\cjkfont 短线：不同分子或有限片段的态}};
\node[font=\small]at(9.15,-.75){{\cjkfont 局域态与有限离域}\quad {\cjkfont 重组、跳跃与连通性}};
''')
keep('exciton_ct',r'''
\node[box,minimum width=3.4cm,minimum height=1.8cm] (a) at(0,0) {{\cjkfont 中性激子}\\{\cjkfont 同一局部环境}};
\node[box,draw=tolorange,fill=tolorange!5,minimum width=3.4cm,minimum height=1.8cm] (b) at(4.7,0) {{\cjkfont 界面}CT\\D$^{+}${\cjkfont 与}A$^{-}${\cjkfont 关联}};
\node[box,draw=tolgreen,fill=tolgreen!4,minimum width=3.4cm,minimum height=1.8cm] (c) at(9.4,0) {{\cjkfont 可输运载流子}\\{\cjkfont 分离与收集}};
\draw[bi] (a)--(b);\draw[bi] (b)--(c);
\node[font=\small] at(4.7,-1.6) {{\cjkfont 光学激发能} $\ne$ CT{\cjkfont 自由能} $\ne$ {\cjkfont 自由电荷端点差}};
''')
keep('stack',r'''
\node[box,minimum width=1.6cm,minimum height=2cm] (a) at(0,0) {ITO};
\node[box,draw=tolgreen,minimum width=2.4cm,minimum height=2cm] (b) at(2.02,0) {PEDOT:PSS};
\node[box,draw=tolorange,fill=tolorange!4,minimum width=3.6cm,minimum height=2cm] (c) at(5.05,0) {D18:L8-BO\\{\cjkfont 约}100 nm};
\node[box,draw=violet,minimum width=1.7cm,minimum height=2cm] (d) at(7.73,0) {PDINN};
\node[box,draw=red!60!black,minimum width=1.5cm,minimum height=2cm] (e) at(9.36,0) {Ag};
\node[font=\small] at(4.6,-1.6) {{\cjkfont 给体∶受体质量比} $1:1.2$\quad {\cjkfont 材料身份暂定}};
''')
keep('srh_cycle',r'''
% Top panels: every arrow denotes an electron transition, not hole motion.
\foreach \xx/\letter/\title in {0/a/电子捕获,3.5/b/电子发射,7/c/空穴捕获,10.5/d/空穴发射}{
 \begin{scope}[xshift=\xx cm]
  \node[font=\small\bfseries] at(1.05,2.9) {(\letter) {\cjkfont \title}};
  \draw[line width=.65pt] (0,2.2)--(2.15,2.2) node[right,font=\small] {$E_C$};
  \draw[dashed,line width=.65pt] (.2,1.1)--(1.95,1.1) node[right,font=\small] {$E_t$};
  \draw[line width=.65pt] (0,0)--(2.15,0) node[right,font=\small] {$E_V$};
 \end{scope}
}
\draw[-{Stealth},line width=.6pt](-.65,0)--(-.65,2.25);
\node[rotate=90,font=\small]at(-1,1.15){{\cjkfont 电子能量}};
\draw[-{Stealth},tolblue,line width=1.15pt](1.05,2.15)--(1.05,1.16);
\draw[-{Stealth},tolblue,line width=1.15pt](4.55,1.16)--(4.55,2.15);
\draw[-{Stealth},tolorange,line width=1.15pt](8.05,1.04)--(8.05,.06);
\draw[-{Stealth},tolorange,line width=1.15pt](11.55,.06)--(11.55,1.04);
\node[font=\small] at(1.05,-.50) {$0\to1$};
\node[font=\small] at(4.55,-.50) {$1\to0$};
\node[font=\small] at(8.05,-.50) {$1\to0$};
\node[font=\small] at(11.55,-.50) {$0\to1$};
\node[font=\small] at(1.05,-1.05) {$r_{nc}=N_t c_n n(1-f)$};
\node[font=\small] at(4.55,-1.05) {$r_{ne}=N_t e_n f$};
\node[font=\small] at(8.05,-1.05) {$r_{pc}=N_t c_p p f$};
\node[font=\small] at(11.55,-1.05) {$r_{pe}=N_t e_p(1-f)$};
\node[font=\small,align=center] at(1.05,-1.63) {{\cjkfont 移走一个}\ {\cjkfont 传导电子}};
\node[font=\small,align=center] at(4.55,-1.63) {{\cjkfont 增加一个}\ {\cjkfont 传导电子}};
\node[font=\small,align=center] at(8.05,-1.63) {{\cjkfont 填补已有价态空穴}};
\node[font=\small,align=center] at(11.55,-1.63) {{\cjkfont 留下一个价态空穴}};
% Bottom arrows denote the center's occupation changes, not energy or space.
\node[box,minimum width=2.8cm,minimum height=1.15cm] (empty) at(1.7,-3.55) {$\circ\quad 0$ {\cjkfont 空态}\\{\cjkfont 概率} $1-f$};
\node[box,draw=tolgreen,minimum width=2.8cm,minimum height=1.15cm] (filled) at(10.9,-3.55) {$\bullet\quad 1$ {\cjkfont 电子占据}\\{\cjkfont 概率} $f$};
\draw[-{Stealth},line width=.85pt,draw=tolblue] (3.2,-3.24)--node[above,font=\small] {{\cjkfont 电子捕获} $c_n n$ {\cjkfont ＋空穴发射} $e_p$}(9.4,-3.24);
\draw[-{Stealth},line width=.85pt,draw=tolorange] (9.4,-3.85)--node[below,font=\small] {{\cjkfont 电子发射} $e_n$ {\cjkfont ＋空穴捕获} $c_p p$}(3.2,-3.85);
''')
keep('charge_states',r'''
\node[box,minimum width=5cm,minimum height=2.1cm] at(0,0) {donor $+/0$\\{\cjkfont 电子发射后残余正核心}\\{\cjkfont 电子吸引尾势}};
\node[box,draw=tolorange,fill=tolorange!4,minimum width=5cm,minimum height=2.1cm] at(6.5,0) {acceptor $0/-$\\{\cjkfont 空穴发射后残余负核心}\\{\cjkfont 空穴吸引尾势}};
''')
keep('weak_regions',r'''
\draw[draw=tolblue,fill=tolblue!4,line width=.9pt] (0,0) rectangle(11,3.5);
\node[box,minimum width=7cm,minimum height=2.2cm] at(3.8,1.8) {{\cjkfont 良好区域}\\{\cjkfont 局部面积电流} $J_{\rm good}$};
\node[box,draw=tolorange,fill=tolorange!5,minimum width=1.8cm,minimum height=2.8cm] at(9.3,1.8) {{\cjkfont 弱区}\\$J_{\rm bad}$};
\node[font=\small] at(5.5,4.1) {{\cjkfont 总面积归一}\quad $J=(1-f)J_{\rm good}+fJ_{\rm bad}$};
''')
keep('eight_state',r'''
\node[box,minimum width=3.5cm,minimum height=1.25cm] (a) at(0,3) {100\\{\cjkfont 价态有电子}};
\node[box,minimum width=3.5cm,minimum height=1.25cm] (b) at(7,3) {010\\{\cjkfont 价空穴＋缺陷电子}};
\node[box,minimum width=3.5cm,minimum height=1.25cm] (c) at(7,0) {001\\{\cjkfont 空间分离束缚对}};
\node[box,minimum width=3.5cm,minimum height=1.25cm] (d) at(0,0) {000\\{\cjkfont 电子已导出}};
\draw[bi] (a)--node[above,font=\small] {{\cjkfont 价态} $\to$ {\cjkfont 缺陷}}(b);
\draw[bi] (b)--node[right,font=\small,align=left] {{\cjkfont 缺陷}\\$\to$ {\cjkfont 传导态}}(c);
\draw[bi] (c)--node[above,font=\small] {{\cjkfont 电子导出}}(d);
\draw[bi] (d)--node[left,font=\small] {{\cjkfont 回填}}(a);
''')
keep('model_closure',r'''
\node[box,minimum width=11cm,minimum height=.8cm] (a) at(0,3.9) {{\cjkfont 微观态与热核\quad 能级 电荷 耦合 振动谱}};
\node[box,draw=tolgreen,minimum width=11cm,minimum height=.8cm] (b) at(0,2.6) {{\cjkfont 完整可逆循环\quad 形成 返回 占据 复位}};
\node[box,draw=tolorange,minimum width=11cm,minimum height=.8cm] (c) at(0,1.3) {{\cjkfont 空间器件与接触\quad 净源 电荷 场 输运}};
\node[box,draw=violet,minimum width=11cm,minimum height=.8cm] (d) at(0,0) {{\cjkfont 可观测量与实验}\quad J--V {\cjkfont 温度 厚度 时间}};
\draw[bi](a)--(b);\draw[bi](b)--(c);\draw[bi](c)--(d);
''')
M=Path('/workspace/scratch/d9cc890e37cd/serif_delivery_tree/OSC_reverse_knee_editable_20261002/data/002高YZ/OSC_FF_reverse_knee_model/model_comparison/results_ff78')
IN=P/'data/latex_plot_inputs';IN.mkdir(parents=True,exist_ok=True)
provenance=[];series=[]
def read(name):
 original=M/name;dest=IN/original.name
 src=dest if dest.exists() else original
 if src!=dest:shutil.copy2(src,dest)
 provenance.append({'input':str(dest.relative_to(P)),'original':str(src),'sha256':hashlib.sha256(src.read_bytes()).hexdigest()})
 return list(csv.DictReader(src.open(encoding='utf-8-sig')))
def coords(name,rows,x,y):
 a=[(float(x(r)),float(y(r))) for r in rows]
 series.append({'series':name,'points':a,'count':len(a),'sha256':hashlib.sha256(json.dumps(a).encode()).hexdigest()})
 return ' '.join('('+format(x,'.17g')+','+format(y,'.17g')+')' for x,y in a)
rows=read('comparison_curves.csv')
defs={'local':(0,'局域SRH'),'pf':(1e15,'有限PF'),'hopping':(1e15,'MA跳跃'),'contact':(.6,'接触'),'thermal':(.0002,'电热'),'weak':(.01,'弱区'),'empirical':(.05,'经验')}
colors=['tolblue','tolgreen','tolorange','violet','red!70!black','gray','olive'];styles=['solid','solid','dashed','dashdotted','densely dashed','dotted','loosely dashed']
s=r'''\begin{tikzpicture}\begin{semilogyaxis}[width=.90\linewidth,height=.36\textheight,xmin=0,xmax=30,ymin=15,ymax=3000,xlabel={{\cjkfont 反偏电压幅值} / V},ylabel={{\cjkfont 光照电流幅值} / $\mathrm{mA\,cm^{-2}}$},axis lines=left,legend style={font=\scriptsize,at={(0.02,.98)},anchor=north west,legend columns=4},label style={font=\small},tick label style={font=\small}]
'''
for i,(model,(value,label)) in enumerate(defs.items()):
 a=[r for r in rows if r['model']==model and r['group']=='family' and float(r['d_nm'])==100 and float(r['bath_K'])==300 and float(r['parameter_value'])==value and float(r['light'])==1 and float(r['V'])<=0];a.sort(key=lambda r:abs(float(r['V'])))
 assert a
 s+=r'\addplot[color='+colors[i]+','+styles[i]+r',line width=.9pt,no marks] coordinates {'+coords('overview_'+model,a,lambda r:abs(float(r['V'])),lambda r:abs(float(r['J']))*1000)+r'};\addlegendentry{{\cjkfont '+label+'}}\n'
s+=r'\addplot[black,dashed,line width=.5pt,forget plot] coordinates {(0,50)(30,50)};\end{semilogyaxis}\end{tikzpicture}'
D['overview']=s
rows=read('pf_density_scan/sequence.csv');valid=[r for r in rows if r['V50_V']]
s=r'''\begin{tikzpicture}\begin{groupplot}[group style={group size=2 by 1,horizontal sep=1.3cm},width=.43\linewidth,height=.32\textheight,xmode=log,axis lines=left,label style={font=\small},tick label style={font=\small},xlabel={{\cjkfont 有效} $N_t$ / $\mathrm{cm^{-3}}$}]
\nextgroupplot[ylabel={FF / \%}]
\addplot[tolblue,mark=*,mark size=1.5pt,line width=.8pt] coordinates {'''+coords('density_FF',rows,lambda r:float(r['Nt_cm3']),lambda r:float(r['FF_percent']))+r'''};
\nextgroupplot[ylabel={$|V_{50}|$ / V}]
\addplot[tolorange,mark=square*,mark size=1.5pt,line width=.8pt] coordinates {'''+coords('density_V50',valid,lambda r:float(r['Nt_cm3']),lambda r:-float(r['V50_V']))+r'''};
\end{groupplot}\end{tikzpicture}'''
D['density']=s
(P/'source/native_figures.json').write_text(json.dumps(D,ensure_ascii=False,indent=2))
(P/'qa/native_figure_provenance.json').write_text(json.dumps({'schematics':list(D)[:8],'plots':['overview','density'],'inputs':provenance,'series':series,'scope':'same frozen CSV selections and unit transforms as the previous verified plotting script; no fitting, interpolation or new physics scans','line_styles':'additional dash and marker distinctions for accessibility; numeric points unchanged'},ensure_ascii=False,indent=2))
print('8 TikZ schematics +2 pgfplots charts; curve points',sum(x['count'] for x in series))
