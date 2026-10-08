from analyze_density import *
import os
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'checkpoints/matplotlib'))
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
FONT='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'
plt.rcParams.update({'font.family':FontProperties(fname=FONT).get_name(),'axes.unicode_minus':False,'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':140,'savefig.dpi':180})
COL=['#2673a5','#df8c22','#8c4ea1'];nts=[1e14,1e15,1e16];Nlab=[r'$10^{14}$',r'$10^{15}$',r'$10^{16}$'];cache={}
for L in [0,1]:
 for N in nts:
  r=read('PF',N,L)['rows'];d=read('DD',N,L)['rows'];u=np.array([x['U'] for x in r]);I=-np.array([x['J_Acm2'] for x in r]);Id=-np.array([x['J_Acm2'] for x in d]);cache[N,L]=(u,I,I-Id,r)

def style(ax,x='反偏幅值 U = −V（V）',y=''):
 ax.set_xlabel(x);ax.set_ylabel(y);ax.grid(alpha=.18)
fig,axs=plt.subplots(2,2,figsize=(12.6,9.2),layout='constrained')
for i,N in enumerate(nts):
 for L,ls in [(0,'-'),(1,'--')]:
  u,I,ex,r=cache[N,L];axs[0,0].semilogy(u[1:],I[1:]*1000,ls,color=COL[i],label=Nlab[i]+(' 暗态' if not L else ' 光照'))
  valid=ex>1e-10;axs[0,1].semilogy(u[valid],ex[valid]*1000*1e15/N,ls,color=COL[i],label=Nlab[i]+(' 暗' if not L else ' 光'))
for a in axs.flat:style(a)
axs[0,0].set(title='A  密度改变幅值与电流交点',ylabel='提取电流 −J（mA/cm²）',xlim=(0,30),ylim=(1e-6,4e4));axs[0,0].axhline(50,color='#777777',lw=1,ls=':');axs[0,0].text(2,65,'50 mA/cm²人为交点',fontsize=9);axs[0,0].legend(ncol=2,fontsize=9)
axs[0,1].set(title='B  扣除同密度DD后，归一曲线接近重合',ylabel='(I_PF − I_DD) × 10^15/Nt（mA/cm²）',xlim=(4,30),ylim=(1e-3,4e3));axs[0,1].text(5,400,'10–20 V：每中心源偏差 ≤0.32%\n虚线：光照；实线：暗态',fontsize=10)
a=json.loads((ROOT/'data/analysis.json').read_text());x=np.arange(3)
for L,dx,c in [(0,-.16,'#425b6e'),(1,.16,'#e3b35b')]:
 vals=[next(t['U'] for t in a['thresholds'] if t['Nt_cm3']==N and t['light']==L and t['current_mAcm2']==50) for N in nts];axs[1,0].bar(x+dx,vals,.3,label='暗态' if L==0 else '光照',color=c)
 for k,v in enumerate(vals):axs[1,0].text(k+dx,v+.25,f'{v:.2f}',ha='center',fontsize=9)
axs[1,0].axhspan(23.5,23.75,color='#5a8268',alpha=.2,label='曲率标记范围：非启动阈值');axs[1,0].set_xticks(x,Nlab);axs[1,0].set(title='C  V50移动约8 V，形状标记未见对应移动',xlabel='Nt（cm^-3）',ylabel='反偏幅值 U（V）',ylim=(0,27));axs[1,0].legend(fontsize=8,loc='lower left')
for model,ls in [('PF','-o'),('DD','--s')]:
 vals=[next(t['FF'] for t in a['forward'] if t['Nt_cm3']==N and t['model']==model)*100 for N in nts];axs[1,1].plot(x,vals,ls,label=model,color='#8c4ea1' if model=='PF' else '#85939b');
 for k,v in enumerate(vals):
  if model=='PF':axs[1,1].text(k,v-.9,f'{v:.2f}%',ha='center',fontsize=10)
axs[1,1].set_xticks(x,Nlab);axs[1,1].set(title='D  同时计算正向FF，不重新调回78%',xlabel='Nt（cm^-3）',ylabel='实际FF（%）',ylim=(66,83));axs[1,1].legend(fontsize=9)
fig.suptitle('这组常温PF模型中，缺陷密度主要放大电流',fontsize=18)
fig.supxlabel('300 K · 100 nm · FF78固定输运参数 · 现象学有限距离PF + 自洽Poisson/DD · 81节点 · 高场恒温外推不代表承流能力',fontsize=10)
fig.savefig(ROOT/'figures/density_amplitude_vs_shape.png');plt.close(fig)

fig,axs=plt.subplots(2,2,figsize=(12.6,9.2),layout='constrained')
for i,N in enumerate(nts):
 for L,ls in [(0,'-'),(1,'--')]:
  u,I,ex,r=cache[N,L];pair=np.array([v['pair_source_Acm2'] for v in r]);axs[0,0].plot(u,pair*1000*1e15/N,ls,color=COL[i])
  ar=np.load(ROOT/'data'/f'PF_Nt{N:.0e}_L{L}_N81_h0.25.npz');k=np.argmin(abs(ar['V']+20));xx=np.linspace(0,100,81);axs[1,0].plot(xx,ar['field'][k]*-1/1e6,ls,color=COL[i]);axs[1,1].plot(xx,ar['f'][k],ls,color=COL[i],label=Nlab[i]+(' 暗' if not L else ' 光'))
 axs[0,1].plot(u,100*(np.array([r['pair_source_Acm2'] for r in cache[N,1][3]])/np.array([r['pair_source_Acm2'] for r in cache[N,0][3]])-1),color=COL[i],label=Nlab[i])
axs[0,0].set_yscale('symlog',linthresh=.001);axs[0,0].set(title='A  光照会改变占据和净对源',xlim=(0,30),ylim=(-.03,5e3));style(axs[0,0],y='净对源 × 10^15/Nt（mA/cm²）');axs[0,0].text(1,200,'低偏压亮态可为净复合\n不能把光照视为始终只加常数',fontsize=10)
axs[0,1].set(title='B  高反偏下暗亮净对源逐渐接近',xlim=(10,25),ylim=(-2.6,.1));style(axs[0,1],y='(光照净对源 / 暗态净对源 − 1)（%）');axs[0,1].legend(fontsize=9)
axs[1,0].set(title='C  −20 V局域电场：大体一致，仍有反馈');style(axs[1,0],x='位置 x（nm）',y='−F（MV/cm）')
axs[1,1].set(title='D  −20 V占据：均值相同不等于局域电荷为零',ylim=(-.05,1.05));style(axs[1,1],x='位置 x（nm）',y='陷阱占据 f');axs[1,1].legend(ncol=2,fontsize=8)
fig.suptitle('“暗亮形状接近”是有条件结果，并非严格不受光照影响',fontsize=17)
fig.supxlabel('300 K · 100 nm · 实线暗态 / 虚线光照 · 这些是模型内部量，未据此认定真实机制',fontsize=10);fig.savefig(ROOT/'figures/field_occupancy_and_light.png');plt.close(fig)

fig,axs=plt.subplots(1,2,figsize=(12.6,4.8),layout='constrained')
u,I,ex,r=cache[1e15,0]
for w,col in [(3,'#2673a5'),(5,'#df8c22'),(7,'#8c4ea1')]:
 ss,s1,s2=shape(u,ex,w);axs[0].plot(u,s2/max(abs(ex)),color=col,label=f'{w}点，跨度{(w-1)*.25:g} V')
axs[0].set(xlim=(18,27),title='最大正曲率临近模型势垒饱和',ylabel='归一二阶导数（V^-2）');style(axs[0],y='归一二阶导数（V^-2）');axs[0].legend(fontsize=9)
ss,s1,s2=shape(u,ex,5);axs[1].plot(u,s1*1e3,color='#425b6e');axs[1].axvline(23.5,color='#df8c22',ls=':');axs[1].set(xlim=(18,28),title='±2 V后的斜率已经回落，严格单膝判据不通过',ylabel='超额电流斜率（mA cm^-2 V^-1）');style(axs[1],y='超额电流斜率（mA cm^-2 V^-1）');axs[1].text(18.4,440,'不把23.5 V改名为物理启动电压\n也不把平台改名为击穿',fontsize=10)
fig.suptitle('形状指标的敏感性与拒绝条件也必须保留',fontsize=17);fig.supxlabel('Nt = 10^15 cm^-3 · 300 K · 暗态 · 曲率峰是描述性标记，不是材料常数',fontsize=10);fig.savefig(ROOT/'figures/shape_definition_sensitivity.png');plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(12.6,4.8),layout='constrained')
for i,N in enumerate(nts):
 for L,ls in [(0,'-'),(1,'--')]:
  u,I,ex,r=cache[N,L];ss,s1,s2=shape(u,ex,5);axs[0].plot(u,s2*1e3*1e15/N,ls,color=COL[i],label=Nlab[i]+(' 暗' if not L else ' 光'));valid=ex>1e-10;axs[1].plot(u[valid],s1[valid]/ex[valid],ls,color=COL[i])
style(axs[0],y='归一超额曲率（mA cm^-2 V^-2）');axs[0].set(xlim=(2,8),ylim=(0,.015),title='约5 V还存在较小的局部形状特征');axs[0].legend(ncol=2,fontsize=8)
style(axs[1],y='d ln(超额电流)/dU（V^-1）');axs[1].set(xlim=(2,24),ylim=(0,3),title='对数斜率并非只在一个电压突然开启')
fig.suptitle('局部特征与全局最大曲率不应混为同一膝点',fontsize=17);fig.supxlabel('300 K · 81节点 · 暗亮与窗口改变后，小特征约4.75–5.25 V；严格2 V背景显著性不通过，不能据此认证膝点',fontsize=10);fig.savefig(ROOT/'figures/early_shape_feature.png');plt.close(fig)
