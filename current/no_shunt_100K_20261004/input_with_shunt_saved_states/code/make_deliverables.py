import os
os.environ['MPLCONFIGDIR']='/tmp/pm6_latest_mpl'
from pathlib import Path
import json,csv,hashlib,sys
import numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
ROOT=Path(__file__).resolve().parents[1];TT=[100,125,150,175,200,225,250,275,300]
rr={T:json.loads((ROOT/f'data/latest_T{T}_N321.json').read_text()) for T in TT}
assert all(r['success'] for r in rr.values())
font=FontProperties(fname='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc')
plt.rcParams.update({'font.family':font.get_name(),'font.size':11,'axes.unicode_minus':False,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.17,'figure.facecolor':'white','axes.titleweight':'bold'})
# Explicit registration supports all Chinese captions.
from matplotlib import font_manager
font_manager.fontManager.addfont(font.get_file())
colors={T:plt.get_cmap('coolwarm')((T-100)/200) for T in TT}
rows=[];metrics=[];qa=[]
for T,r in rr.items():
 m=dict(r['metrics']);m['FF_percent']=m.pop('FF')*100;m['temperature_status']='conditional extrapolation' if T<225 else 'original warm branch reproduced';metrics.append(m)
 for p in r['curves']['light']:rows.append(dict(T_K=T,V_V=p['V'],J_terminal_mAcm2=p['J_terminal_Acm2']*1000,J_intrinsic_mAcm2=p['J_intrinsic_Acm2']*1000,J_shunt_mAcm2=p['J_shunt_Acm2']*1000,output_power_mWcm2=-p['external_terminal_Wcm2']*1000,nonlinear_residual=p['scaled_nonlinear_residual'],current_closure_uncertainty_Acm2=p['current_closure_uncertainty_Acm2'],current_spread_Acm2=p['current_spread_Acm2'],charge_relative=p['charge_relative'],max_density_fraction=p['max_density_fraction'],max_core_density_fraction=p['max_core_density_fraction'],field_capped_length_fraction=p['field_capped_length_fraction'],intrinsic_gate_passed=p['gate_passed'],temperature_status=m['temperature_status']))
 qa.append(dict(T_K=T,**r['audit_summary'],Voc_resolved_bracket=m['Voc_resolved_bracket'],Voc_bracket_half_width_V=m['Voc_bracket_half_width_V']))
def writecsv(path,rows):
 with path.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
writecsv(ROOT/'PM6Y6_latest_100-300K_JV.csv',rows);writecsv(ROOT/'PM6Y6_latest_100-300K_PV_parameters.csv',metrics);writecsv(ROOT/'validation/numerical_QA.csv',qa)
mesh=[]
for T in [100,150,200]:
 q=json.loads((ROOT/f'data/latest_T{T}_N641.json').read_text());assert q['success']
 for k in ['Jsc_mAcm2','Voc_V','FF','Pmax_mWcm2','Vmp_V','Jmp_mAcm2']:
  a=rr[T]['metrics'][k];b=q['metrics'][k];mesh.append(dict(T_K=T,quantity=k,N321=a,N641=b,difference=b-a,relative_difference=(b-a)/b))
writecsv(ROOT/'validation/mesh_321_641.csv',mesh)
c=json.loads((ROOT/'data/latest_T100_N321_cap1e+06.json').read_text());assert c['success']
cap=[]
for k in ['Jsc_mAcm2','Voc_V','FF','Pmax_mWcm2']:
 a=rr[100]['metrics'][k];b=c['metrics'][k];cap.append(dict(T_K=100,quantity=k,nominal_cap235000=a,check_cap1000000=b,relative_change=(b-a)/a))
writecsv(ROOT/'validation/100K_field_cap_sensitivity.csv',cap)
def curve(T):
 p=rr[T]['curves']['light'];return np.array([x['V'] for x in p]),np.array([x['J_terminal_Acm2']*1000 for x in p])
fig=plt.figure(figsize=(14.6,8.4),layout='constrained');gs=fig.add_gridspec(2,2,width_ratios=[1.72,1]);ax=fig.add_subplot(gs[:,0]);a125=fig.add_subplot(gs[0,1]);a100=fig.add_subplot(gs[1,1])
for T in TT[::-1]:
 v,j=curve(T);ax.plot(v,j,color=colors[T],lw=2.3,ls='--' if T<225 else '-',label=f'{T} K'+('  外推' if T<225 else ''))
ax.axhline(0,color='.35',lw=.8);ax.set(xlim=(0,1.06),ylim=(-28.5,2.5),xlabel='端电压 V (V)',ylabel='有符号端电流密度 J (mA/cm²)',title='100–300 K 全温区 | 全部线性坐标');ax.legend(ncol=2,loc='lower right',fontsize=10,framealpha=.93)
for T,aa,vs,js,xlab,ylab in [(125,a125,1,1,'端电压 V (V)','J (mA/cm²)'),(100,a100,1000,1000,'端电压 V (mV)','J (µA/cm²)')]:
 v,j=curve(T);aa.plot(v*vs,j*js,color=colors[T],lw=2.4,ls='--');aa.axhline(0,color='.35',lw=.8)
 m=rr[T]['metrics'];aa.scatter([m['Voc_V']*vs],[0],s=38,facecolor='white',edgecolor=colors[T],zorder=4)
 aa.set(xlim=(0,v.max()*vs),xlabel=xlab,ylabel=ylab,title=f'{T} K 放大 | 仍为线性坐标')
 aa.text(.035,.96,f"Jsc = {m['Jsc_mAcm2']*js:.6g} {'µA/cm²' if T==100 else 'mA/cm²'}\nVoc = {m['Voc_V']*vs:.6g} {'mV' if T==100 else 'V'}",transform=aa.transAxes,va='top',fontsize=10,bbox=dict(facecolor='white',alpha=.85,edgecolor='none'))
fig.suptitle('PM6:Y6 最新室温校准版：自洽 full-FD J–V',fontsize=19,fontweight='bold')
fig.supxlabel('固定同一组 300 K 校准；低于 225 K 为条件性外推。100 K 曲线为本轮新求解，非旧版结果。',fontsize=11)
fig.savefig(ROOT/'PM6Y6_latest_100-300K_LINEAR_JV.png',dpi=200);plt.close(fig)
fig,aa=plt.subplots(2,2,figsize=(12.2,8.3),layout='constrained')
for ax,(key,label,ylim) in zip(aa.ravel(),[('Jsc_mAcm2','Jsc (mA/cm²)',(0,29)),('Voc_V','Voc (V)',(0,1.06)),('FF','FF (%)',(0,75)),('Pmax_mWcm2','Pmax (mW/cm²)',(0,16))]):
 fac=100 if key=='FF' else 1;y=np.array([rr[T]['metrics'][key]*fac for T in TT]);ax.axvspan(225,300,color='#d6e5f4',alpha=.5)
 ax.plot(TT[:6],y[:6],'o--',color='#285d9a',lw=2);ax.plot(TT[5:],y[5:],'o-',color='#ba3d45',lw=2)
 ax.set(xlim=(95,305),ylim=ylim,xlabel='温度 T (K)',ylabel=label);ax.set_xticks(TT)
 if key=='Voc_V':ax.annotate('100 K: 1.383 mV',xy=(100,y[0]),xytext=(110,.15),arrowprops=dict(arrowstyle='->',color='.3'),fontsize=10)
fig.suptitle('最新校准版的光伏参数 | 全部线性坐标',fontsize=18,fontweight='bold')
fig.supxlabel('蓝色虚线：条件性低温外推；红色实线/阴影：原 225–300 K 分支原样复现。未逐温度重新拟合。',fontsize=10.5)
fig.savefig(ROOT/'PM6Y6_latest_100-300K_LINEAR_PV.png',dpi=200);plt.close(fig)
summary={'main_runs':len(rr),'JV_rows':len(rows),'all_main_intrinsic_gates_passed':all(q['all_intrinsic_gates'] for q in qa),'all_Voc_brackets_resolved':all(q['Voc_resolved_bracket'] for q in qa),'max_residual':max(q['max_scaled_nonlinear_residual'] for q in qa),'max_core_occupation':max(q['max_core_density_fraction'] for q in qa),'max_mesh_relative_Jsc':max(abs(r['relative_difference']) for r in mesh if r['quantity']=='Jsc_mAcm2'),'max_mesh_relative_Pmax':max(abs(r['relative_difference']) for r in mesh if r['quantity']=='Pmax_mWcm2'),'cap100_relative_Jsc':cap[0]['relative_change'],'cap100_relative_Pmax':cap[3]['relative_change']}
(ROOT/'validation/summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
text='''PM6:Y6 最新室温校准版：100–300 K 条件性外推\n2026-10-04\n\n本轮新求解，不采用较早版本的100 K结果。九个温度为100、125、150、175、200、225、250、275、300 K。所有图的横纵轴均为线性。绘制有符号端电流（光伏发电区J<0），并提供125 K和100 K放大图。\n\n1. 冻结的最新校准\n300 K：G=1.5357919101745596e22 cm^-3 s^-1，Eg=1.4210191871124291 eV，β_eff=4.5633944521269105e-12 cm^3/s；μn=8.4e-4，μp=1.3e-4 cm²/Vs；σn=60、σp=74 meV；N0=2.4e20 cm^-3；d=110 nm；εr=3.5。电子场系数γn=0.001377，空穴室温γp=0.001145 (cm/V)^1/2。Rs=0，Rsh=3575.0455923806257 Ω cm²。保持理想多数载流子库、少数阻挡边界，不逐温度改接触、带隙、光生率或并联电阻。\n225–300 K的原主分支rt_author_hole精确保留；本次重新求解与原始主要参数的差异在数值精度内。原暖温源包未覆盖。\n\n2. 新增且明确声明的低温假设\nμ(T)=μ300 exp[-(4/9)(σ/kB)²(T^-2-300^-2)]。迁移率来源温区约223–328 K；100–200 K明确开启外推，数值full-FD可算不等于材料迁移率已获验证。\nβ：200–300 K保持作者复合拟合线的相对形状、锚定最新β300；200 K以下按最后两个作者节点200/225 K的Arrhenius斜率延伸，Ea=0.0911727644 eV。β(T)=β(200) exp[-Ea/kB(1/T-1/200)]。这是一种可复现的条件性选择，来源没有100–175 K β测量。\nγp：225 K以下用GDM启发的γp=a+b/T²，连续锚定225 K，b由225/250 K两个已用节点决定，为240.334685 (cm/V)^1/2 K²；γn继续沿用原主分支的室温常数。该低温场系数并非原文测量或唯一可识别外推。\n局域场增强仍为μ(F)=μ0 exp[γ sqrt(|F|)]的原正则化实现，并保留235000 V/cm的原场上限。没有偷偷把γ或β钳为常数，也未任意加新密度迁移率/陷阱通道。\n\n3. 主要光伏参数\nT (K), Jsc (mA/cm²), Voc (V), FF (%), Pmax (mW/cm²)\n'''
for T in TT:
 m=rr[T]['metrics'];text+=f"{T}, {m['Jsc_mAcm2']:.9g}, {m['Voc_V']:.9g}, {100*m['FF']:.7g}, {m['Pmax_mWcm2']:.9g}\n"
text+='''\n100 K：输运显著冻结，μn0=4.04135e-12、μp0=2.89546e-17 cm²/Vs；Jsc约3.87436e-4 mA/cm²。维持室温Rsh使端Voc降至约1.383 mV，功率曲线近线性漏电负载，FF约25%。这不是100 K材料或器件的已验证预测。低温Voc回落不应解释为电子态统计数值失败；也不能据此断言实测器件会发生同样回落。\n\n4. 数值检查与不确定性\n自洽求解Poisson、双载流子连续性、Gaussian DOS full Fermi–Dirac积分及广义Einstein通量；保留带详细平衡的复合项。从300 K递降温度，以前一温度解作延拓初值，必要时支持自适应温度/电压桥接。每个温度覆盖0 V至略高于端Voc，根和最大功率点由PDE再次求值，未只做曲线插值。\n九个主结果全部通过原电流守恒、空间电荷和非线性残差门槛；全部正负端Voc夹逼可分辨。守恒门槛只衡量数值自洽，不表示模型获得实验验证。零端电流附近必须区分固有器件电流与并联漏电，两者已分别保存在CSV。\n'''
text+=f"主运行最大缩放非线性残差={summary['max_residual']:.4g}。核心区域最大占据率={summary['max_core_occupation']:.6g}；理想接触边界约0.5，仍按完整FD处理。\n"
text+=f"200/150/100 K用321→641节点核验：Jsc最大相对变化{summary['max_mesh_relative_Jsc']*100:.5g}%，Pmax最大相对变化{summary['max_mesh_relative_Pmax']*100:.5g}%。\n"
text+=f"100 K把场上限提高至1e6 V/cm的单项敏感性：Jsc变化{summary['cap100_relative_Jsc']*100:.4g}%，Pmax变化{summary['cap100_relative_Pmax']*100:.4g}%；主曲线仍采用原上限。\n"
text+='''100 K主分支短路态场上限在约3.18%厚度内激活，最大空穴场增强约6.74e4；这是γ外推的重要条件，并不声称其高场有效。该测试不能囊括β、γ、接触、温变Rsh/G等物理不确定性。\n\n5. 来源与边界\nPerdigón-Toro等，Advanced Energy Materials 12, 2103422 (2022)，https://doi.org/10.1002/aenm.202103422；原文及SI的SCLC/材料参数、作者复合拟合线审计表和最新室温冻结拟合均随计算参考表保留。作者β的200–300 K节点来自其学位论文Fig.8.3a。\nBässler, physica status solidi (b) 175, 15–56 (1993)，https://doi.org/10.1002/pssb.2221750102；GDM场温度形式仅作为延伸γ的物理动机，不是本器件低温定量验证。\n最新暖温分支原已有模型/独立文献温变Voc不一致，未因此重新逐温度拟合。请将此图视为固定最新校准加明确外推假设的模型结果，而非实验曲线或唯一低温预测。\n\n文件：两张LINEAR PNG、完整JV CSV、PV参数CSV、代码/状态/数值检查归档。未改变Library历史备份；旧模型从活动工作目录移入可恢复归档，恢复说明与清单另存。\n'''
(ROOT/'PM6Y6_latest_100-300K_说明.txt').write_text(text)
(ROOT/'reproduce.sh').write_text('#!/bin/sh\nset -eu\ncd "$(dirname "$0")"\nOPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python code/extend_device.py\nOPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python code/extend_device.py --N 641 --temperatures 200 150 100\nOPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python code/extend_device.py --N 321 --temperatures 100 --cap 1000000\npython code/make_deliverables.py\n')
