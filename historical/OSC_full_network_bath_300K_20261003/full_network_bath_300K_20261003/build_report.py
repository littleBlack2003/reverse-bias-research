from pathlib import Path
import json,os
os.environ['MPLBACKEND']='Agg';os.environ['MPLCONFIGDIR']='/tmp/bath_mpl'
import pandas as pd,numpy as np,matplotlib.pyplot as plt
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,Image,PageBreak
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.cidfonts import UnicodeCIDFont
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
R=Path(__file__).resolve().parent;D=R/'results'
s=pd.read_csv(D/'summary.csv',dtype={'most_likely_state':str}).set_index('case');v=json.loads((D/'validation.json').read_text());b=s.loc['baseline_refined'];q=s.loc['full_quantum_refined']
bs=pd.read_csv(D/'baseline_refined_states.csv',dtype={'state':str}).set_index('state');qs=pd.read_csv(D/'full_quantum_refined_states.csv',dtype={'state':str}).set_index('state')
be=pd.read_csv(D/'baseline_refined_edges.csv',dtype={'source':str,'destination':str});qe=pd.read_csv(D/'full_quantum_refined_edges.csv',dtype={'source':str,'destination':str})
ratio=qe.forward_rate_s/be.forward_rate_s
fig,ax=plt.subplots(1,3,figsize=(12,3.5),layout='constrained')
col=['#547b9b','#d37a37'];labs=['Gaussian slow bath','Quantum slow bath']
ax[0].bar(labs,[b.R_s,q.R_s],color=col);ax[0].set_ylabel('Net turnover per graph (s$^{-1}$)');ax[0].set_title('A  Entire network');ax[0].tick_params(axis='x',labelsize=8)
for i,y in enumerate([b.R_s,q.R_s]):ax[0].text(i,y,f'{y:,.1f}',ha='center',va='bottom',fontsize=9)
ax[0].set_ylim(0,max(b.R_s,q.R_s)*1.2)
keys=bs.probability.nlargest(4).index;x=np.arange(len(keys));ax[1].bar(x-.18,bs.loc[keys].probability,.36,color=col[0]);ax[1].bar(x+.18,qs.loc[keys].probability,.36,color=col[1]);ax[1].set_xticks(x,keys);ax[1].set_ylabel('Steady probability');ax[1].set_title('B  Main residence states')
ax[2].scatter(be.energy_change_eV,ratio,c=['#547b9b' if not p.startswith('bath') else '#d37a37' for p in be.process],s=20);ax[2].axhline(1,color='gray',lw=1);ax[2].set_xlabel('Oriented state energy change (eV)');ax[2].set_ylabel('Quantum / Gaussian forward rate');ax[2].set_title('C  All 64 reversible edges');ax[2].set_ylim(.90,1.11)
fig.suptitle('300 K | 1.5 MV/cm | same 12 nm, 32-state graph | no material fit',fontsize=12)
fig.savefig(R/'network_comparison.png',dpi=180);plt.close(fig)
rows=[['量','经典慢浴基线','全量子慢浴'],['净周转率 R / s-1',f'{b.R_s:.6f}',f'{q.R_s:.6f}'],['最可能状态',str(b.most_likely_state),str(q.most_likely_state)],['该态概率',f'{b.p_most_likely:.8f}',f'{q.p_most_likely:.8f}'],['最大总离开率 / s-1',f'{b.maximum_exit_rate_s:.6g}',f'{q.maximum_exit_rate_s:.6g}'],['声子浴热 / eV s-1',f'{b.phonon_heat_eV_s:.6f}',f'{q.phonon_heat_eV_s:.6f}'],['电子库热 / eV s-1',f'{b.electron_bath_heat_eV_s:.6f}',f'{q.electron_bath_heat_eV_s:.6f}'],['总热=化学功 / eV s-1',f'{b.total_bath_heat_eV_s:.6f}',f'{q.total_bath_heat_eV_s:.6f}']]
lead=f'整网替换后，R 从 {b.R_s:.2f} 变为 {q.R_s:.2f} s⁻¹，改变 {(v["rate_ratio"]-1)*100:+.3f}%。单边速率比范围{ratio.min():.4f}–{ratio.max():.4f}，整网响应远小于单边最大变化，说明该点的稳态占据与串行网络显著减弱了净周转响应。这是固定模型下的谱敏感性，不能当成 D18:L8-BO 的材料预测。'
methods=[
('只比较一个物理点','300 K、F=1.5 MV/cm、显式12 nm五位点32态、η=0.5保持不变。两端电化学势差1.5 eV，静电降1.8 V，两者是原协议的独立输入；不能相加为两个电源。平衡检验只将两库设为共同μ=0.023 eV，保留同一个静电场。'),
('改变的唯一物理输入','经典 Gaussian 慢浴换为已有全量子 Ohmic 指数截止连续浴 Ec=20 meV；λs=50 meV、高频150 meV模式及S=1不变。20 meV是谱场景，不是测量值。全网64条可逆边均更新，包括32条内部边与32条储库边及其反向率和交换能量一阶矩。'),
('共同浴与预算','相同轨道键的ΣH²=4.095153730856492×10⁻¹⁰ eV²，两库各归一1 eV单侧DOS。每次添加/转移具有同一谱可由五位点共同位移 aᵢ=(e₀+eᵢ)/√2 构造，Gram特征值为0.5（四重）与3；这只证明可构造性。状态能沿用充分弛豫能量解释，不额外加入重组能。'),
('不匹配的谱量',f'保持重组能和核总积分，不保持有限温方差。全量子慢浴方差为 {v["refined"]["exact_variance_eV2"]:.8g} eV²，经典为 {2*.05*8.617333262145e-5*300:.8g} eV²。不能把相同λ称为完全相同谱预算。'),
('求解与独立诊断','原算法先独立算正反核/储库卷积，再用精确热力学比组装高精度矩阵。本轮另用未强制比率的双向积分直接组装矩阵，单独检查概率、电流和共同μ残余边流；高精度守恒不是核准确性的独立证据。'),
('一次数值收紧',f'主点：慢浴log P(E)网格2 meV、热侧带±30、Fourier容差2×10⁻¹⁰、分段Gauss-Legendre 24阶、250位稳态。复核：1 meV、±40、10⁻¹²、48阶、300位。两次均覆盖慢浴能量±9 eV，无外推、负率裁剪或截断后重归一。R相对变化 {v["full_quantum_R_refinement_relative"]:.3g}。这属于数值复核，不是新增物理场景。'),
('结论边界','无场扫描、温度扫描或参数拟合；无器件Poisson/漂移扩散、可用图密度、共享轨道、光照/FF闭合。单点不能确定上翘阈值。浴每跳前热化仍是条件；更精确求积不会验证材料热化时间。')]
checks=[f'原始双向详细平衡最大对数误差：{v["max_raw_DB_log_error"]:.3g}',f'全量子核离网格对数插值误差（121点诊断）：{v["refined"]["off_grid_max_log_error"]:.3g}',f'慢浴积分质量：{v["refined"]["slow_density_mass"]:.12g}；热侧带质量：{v["refined"]["sideband_mass"]:.12g}',f'共同μ直接组装最大逐边净流/双向总流：{s.loc["full_quantum_equilibrium_refined","direct_assembly_eq_max_edge_flux_over_gross"]:.3g}',f'直接组装与精确比率组装R差：{q.direct_assembly_R_difference_s:.3g} s⁻¹',f'最小稳态概率：{q.min_probability:.3g}；连续性相对误差：{q.max_node_continuity_relative:.3g}',f'总热/化学功闭合相对误差：{q.heat_chemical_relative:.3g}；完整12 nm位移检查通过']
md='# 300 K整网慢浴替换：固定图单点检验\n\n2026-10-03\n\n'+lead+'\n\n'+'\n'.join('| '+' | '.join(r)+' |' for r in [rows[0],['---']*3]+rows[1:])+'\n\n![全网比较](network_comparison.png)\n\n'
for h,t in methods:md+=f'## {h}\n\n{t}\n\n'
md+='## 验收结果\n\n'+'\n'.join('- '+x for x in checks)+'\n\n## 来源与复现\n\n来源是已授权恢复的空间逃逸审计与振动有效性审计。inputs为源代码原样副本，source_sha256.json为哈希；上游归档来源见SOURCE_PROVENANCE.json。run_comparison.py只计算本轮一个物理点及数值/平衡诊断；build_report.py生成此报告和图。依赖NumPy、SciPy、mpmath、Pandas、Matplotlib、ReportLab。运行 python -B run_comparison.py，然后 python -B build_report.py。results含7组完整状态、逐边净流、热账本和汇总，旧交付不变。\n'
(R/'REPORT_300K.md').write_text(md)
pdfmetrics.registerFont(UnicodeCIDFont('STSong-Light'));styles=getSampleStyleSheet();styles.add(ParagraphStyle(name='ZH',fontName='STSong-Light',fontSize=10,leading=15,spaceAfter=7));styles.add(ParagraphStyle(name='ZHHead',fontName='STSong-Light',fontSize=14,leading=19,spaceAfter=9));styles.add(ParagraphStyle(name='ZHTitle',fontName='STSong-Light',fontSize=21,leading=28,spaceAfter=15))
def clean(t):
 for a,b in {'⁻¹⁰':'^-10','⁻¹²':'^-12','⁻¹':'^-1','²':'^2','Σ':'sum ','ᵢ':'_i','₀':'_0','√2':'sqrt(2)','•':'-','–':'-'}.items():t=t.replace(a,b)
 return t
story=[];P=lambda t:Paragraph(clean(t),styles['ZH'])
story+=[Paragraph('300 K 整网慢浴替换',styles['ZHTitle']),P('固定图单点检验 | 2026-10-03'),P(lead),Spacer(1,8)]
t=Table([[P(c) for c in row] for row in rows],colWidths=[185,150,150]);t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e8eff5')),('LINEBELOW',(0,0),(-1,0),1,colors.HexColor('#547b9b')),('VALIGN',(0,0),(-1,-1),'TOP'),('BOTTOMPADDING',(0,0),(-1,-1),6)]));story+=[t,Spacer(1,13),Image(str(R/'network_comparison.png'),width=500,height=146),P('图A：整网净周转。图B：基线最常驻四态。图C：全部有向代表边的正向率比，蓝为内部跃迁、橙为储库交换；反向边同样更新。状态位序为 Vo、V、D、C、Co。'),PageBreak()]
for h,txt in methods[:5]:story+=[Paragraph(h,styles['ZHHead']),P(txt)]
story+=[PageBreak(),Paragraph('数值验收与结论边界',styles['ZHTitle'])]
for h,txt in methods[5:]:story+=[Paragraph(h,styles['ZHHead']),P(txt)]
story+=[Paragraph('独立核查',styles['ZHHead'])]+[P('• '+x) for x in checks]
story+=[Paragraph('来源与复现',styles['ZHHead']),P('使用已恢复的空间逃逸审计和振动有效性审计原样源码。原图基线复现为13094.610062875749 s-1。完整归档包含中文报告、脚本、源代码哈希、逐态/逐边CSV与验收JSON；运行 run_comparison.py 后运行 build_report.py 即可重建。未改动旧交付。')]
def footer(c,d):c.setFont('STSong-Light',9);c.setFillColor(colors.gray);c.drawString(45,26,'条件模型检验；非材料预测');c.drawRightString(550,26,str(d.page))
SimpleDocTemplate(str(R/'OSC_full_network_bath_300K_20261003.pdf'),pagesize=(595.28,841.89),leftMargin=45,rightMargin=45,topMargin=42,bottomMargin=45).build(story,onFirstPage=footer,onLaterPages=footer)
print(lead)
