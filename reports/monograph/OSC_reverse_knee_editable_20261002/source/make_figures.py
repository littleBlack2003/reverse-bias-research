from pathlib import Path
import os,json,shutil,hashlib,csv
P=Path(__file__).resolve().parents[1];os.environ['MPLCONFIGDIR']=str(P/'qa/mpl')
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,Circle,FancyArrowPatch
from matplotlib.font_manager import FontProperties
import numpy as np
font=FontProperties(fname=str(P/'assets/NotoSansSC-Regular.ttf'))
plt.rcParams.update({'font.family':font.get_name(),'font.size':11,'axes.titlesize':12,'axes.labelsize':11,'xtick.labelsize':10,'ytick.labelsize':10,'legend.fontsize':10,'svg.fonttype':'path','axes.spines.top':False,'axes.spines.right':False})
from matplotlib.font_manager import fontManager
fontManager.addfont(str(P/'assets/NotoSansSC-Regular.ttf'))
C=['#215B80','#087F8C','#BD6A27','#775899','#B33C41','#65717B']
def save(fig,name):
 fig.savefig(P/'figures'/f'{name}.png',dpi=250,bbox_inches='tight',facecolor='white');fig.savefig(P/'figures'/f'{name}.svg',bbox_inches='tight',facecolor='white');plt.close(fig)
def canvas(w=7,h=2.7):
 f,a=plt.subplots(figsize=(w,h));a.set(xlim=(0,10),ylim=(0,4));a.axis('off');return f,a
def box(a,x,y,w,h,label,c=0):
 a.add_patch(Rectangle((x,y),w,h,edgecolor=C[c],facecolor=C[c]+'10',lw=1.4));a.text(x+w/2,y+h/2,label,ha='center',va='center',fontproperties=font)
def arrow(a,x1,y1,x2,y2,label=None):
 a.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='<->',mutation_scale=12,color=C[0],lw=1.5));
 if label:a.text((x1+x2)/2,(y1+y2)/2+.25,label,ha='center',fontproperties=font,fontsize=10)
f,a=canvas();box(a,.2,.4,4.2,3.1,'周期势\nBloch态与能带\n散射  波包  有效质量',0);box(a,5.6,.4,4.2,3.1,'分子与非晶环境\n局域态与有限离域\n重组  跳跃  连通性',1);arrow(a,4.45,2,5.55,2);save(f,'bands_localized')
f,a=canvas();
for x,label,col in [(0.15,'中性激子\n同一局部环境',0),(3.6,'界面CT\nD+与A−关联',2),(7.05,'可输运载流子\n分离与收集',1)]:box(a,x,1.1,2.8,1.9,label,col)
arrow(a,3.0,2,3.55,2);arrow(a,6.45,2,7.0,2);a.text(5,.3,'光学激发能 ≠ CT自由能 ≠ 自由电荷端点差',ha='center',fontproperties=font);save(f,'exciton_ct')
f,(a,b)=plt.subplots(1,2,figsize=(7,2.7));e=np.linspace(-.3,.3,300);a.plot(e,np.exp(-e*e/(2*.08**2)),color=C[0]);a.set(xlabel='位点能偏移 / eV',ylabel='示意DOS / a.u.',title='能量分布');rng=np.random.default_rng(21);pts=rng.uniform(0,1,(18,2));
for i in range(18):
 for j in range(i):
  d=np.linalg.norm(pts[i]-pts[j]);
  if d<.35:b.plot([pts[i,0],pts[j,0]],[pts[i,1],pts[j,1]],color=C[1],lw=2*np.exp(-d/.2),alpha=.8)
b.scatter(pts[:,0],pts[:,1],c=np.arange(18),cmap='viridis',s=50);b.set_title('空间连通与耦合');b.axis('off');f.tight_layout();save(f,'disorder_network')
f,a=canvas();labels=['ITO','PEDOT:PSS','D18:L8-BO\n约100 nm','PDINN','Ag'];widths=[1.5,1.65,3.2,1.35,1.3];x=.3
for i,(lab,w) in enumerate(zip(labels,widths)):box(a,x,1,w-.06,2,lab,i%5);x+=w
a.text(5,.35,'给体∶受体质量比 1∶1.2   材料身份暂定',ha='center',fontproperties=font);save(f,'stack')
f,a=plt.subplots(figsize=(7,2.7));v=np.linspace(0,25,300);j=24+35*np.exp((v-15)/4)/(1+np.exp((v-23)/2));a.plot(v,j,color=C[0]);
for y,c in [(40,1),(50,2),(75,3)]:
 a.axhline(y,color=C[c],ls='--',lw=1);ix=np.argmin(abs(j-y));a.scatter(v[ix],j[ix],color=C[c]);a.text(v[ix]-.5,y+4,f'V{y}',color=C[c])
a.set(xlabel='反向电压幅值 / V',ylabel='电流幅值 / mA cm^-2',ylim=(0,165));f.tight_layout();save(f,'thresholds')
f,a=canvas();box(a,3.3,2.4,3.4,1,'中心空态  0',0);box(a,3.3,.35,3.4,1,'中心有电子  1',1);arrow(a,3.7,2.35,3.7,1.4);arrow(a,6.3,2.35,6.3,1.4);a.text(1.9,1.85,'电子捕获 ↓\n电子发射 ↑',ha='center',va='center',fontproperties=font);a.text(8.15,1.85,'空穴发射 ↓\n空穴捕获 ↑',ha='center',va='center',fontproperties=font);save(f,'srh_cycle')
f,a=plt.subplots(figsize=(7,2.7));g=np.logspace(-3,7,500);a.semilogx(g,2*g/(1+g),color=C[0],lw=2);a.axhline(2,ls='--',color=C[2]);a.set(xlabel='单支倍率 g',ylabel='R(g) / R(1)',ylim=(0,2.2));f.tight_layout();save(f,'supply_bottleneck')
f,a=plt.subplots(figsize=(7,2.9));r=np.linspace(.35,6,500);A=1.43996455/3.5
for F,c in [(0,0),(.05,1),(.2,2)]:
 a.plot(r,-A/r-F*r,color=C[c],label=f'F = {10*F:g} MV/cm')
 if F:rs=np.sqrt(A/F);a.scatter(rs,-2*np.sqrt(A*F),color=C[c])
a.axvspan(0,.6,color='#eeeeee');a.set(xlabel='离核心距离 r / nm',ylabel='相对势能 / eV',ylim=(-1.3,.08),xlim=(.3,6));a.legend(loc='upper right');f.tight_layout();save(f,'pf_potential')
f,a=canvas();box(a,.4,.8,4.1,2.5,'donor  + / 0\n电子发射后残余正核心\n电子吸引尾势',0);box(a,5.5,.8,4.1,2.5,'acceptor  0 / −\n空穴发射后残余负核心\n空穴吸引尾势',2);save(f,'charge_states')
f,a=plt.subplots(figsize=(7,2.7));dg=np.linspace(-.9,.4,500);T=.02585;lam=.2;km=np.exp(-(dg+lam)**2/(4*lam*T));kma=np.exp(-np.maximum(dg,0)/T);a.semilogy(dg,km,label='经典Marcus  归一化',color=C[0]);a.semilogy(dg,kma,label='MA  归一化',color=C[1]);a.axvline(-lam,color=C[2],ls='--');a.set(xlabel='自由能差 ΔG / eV',ylabel='相对率',ylim=(1e-7,1.5));a.legend();f.tight_layout();save(f,'ma_marcus_shapes')
f,a=plt.subplots(figsize=(7,2.7));x=np.linspace(-.3,1.35,400);a.plot(x,.6*x*x,label='初态自由能面',color=C[0]);a.plot(x,.6*(x-1)**2-.15,label='末态自由能面',color=C[1]);a.axvline(.375,color=C[2],ls='--');a.annotate('交叉构型',(.375,.084),(.62,.35),arrowprops={'arrowstyle':'->'},fontproperties=font);a.set(xlabel='集体核坐标 Q / arbitrary',ylabel='自由能 / eV',ylim=(-.2,.65));a.legend();f.tight_layout();save(f,'marcus_parabolas')
f,a=canvas();box(a,.8,.5,8.4,2.8,'',0);box(a,1,1.1,5.8,1.8,'良好区域\n局部面积电流 Jgood',0);box(a,7,.65,1.2,2.5,'弱区\nJbad',2);a.text(5,3.65,'总面积归一  J = (1−f) Jgood + f Jbad',ha='center',fontproperties=font);save(f,'weak_regions')
f,a=plt.subplots(figsize=(7,2.7));T=np.linspace(77,330,400);k=8.617333262e-5
for w,c in [(.01,0),(.02,1),(.15,2)]:
 x=w/(k*T);a.plot(T,x/2/np.tanh(x/2),label=f'{1000*w:g} meV',color=C[c])
a.set(xlabel='温度 / K',ylabel='量子涨落 / 经典涨落');a.legend();f.tight_layout();save(f,'vibration_quantization')
f,a=canvas();positions=[(1,2.4),(6.2,2.4),(6.2,.35),(1,.35)];labs=['100\n价态有电子','010\n价空穴＋缺陷电子','001\n空间分离束缚对','000\n电子已导出']
for (x,y),lab in zip(positions,labs):box(a,x,y,2.7,1.0,lab,0)
arrow(a,3.8,2.9,6.1,2.9,'价态 → 缺陷');arrow(a,7.55,2.3,7.55,1.45);arrow(a,6.1,.85,3.8,.85,'电子导出');arrow(a,2.35,1.45,2.35,2.3);a.text(.4,1.85,'回填',ha='center',fontproperties=font,fontsize=10);a.text(9.1,1.85,'缺陷 → 传导态',ha='center',fontproperties=font,fontsize=9);save(f,'eight_state')
f,a=canvas(h=3.0);layers=[('微观态与热核','能级 电荷 耦合 振动谱'),('完整可逆循环','形成 返回 占据 复位'),('空间器件与接触','净源 电荷 场 输运'),('可观测量与实验','J–V 温度 厚度 时间')]
for i,(head,sub) in enumerate(layers):
 y=3.1-i*.85;box(a,1.3,y,7.4,.68,head+'   '+sub,i%4)
 if i<3:arrow(a,5,y-.02,5,y-.14)
save(f,'model_closure')
# Add final frozen stages and copy data provenance.
manifest=json.loads((P/'data/source_manifest.json').read_text())
R=Path('/workspace/shared/wxh/research')
newfig={'reservoir_support':R/'quantum_reservoir_audit_20261002/figures/same_budget_spectral_support_control.png','reservoir_temperature':R/'quantum_reservoir_audit_20261002/figures/temperature_and_pair_accumulation.png','reservoir_heat':R/'quantum_reservoir_audit_20261002/figures/physical_heat_partition.png','lowT_cancellation':R/'qf_temperature_confidence_20261002/figures/cancellation_confidence.png'}
for name,src in newfig.items():
 dst=P/'figures'/f'{name}.png';shutil.copy2(src,dst);manifest.append({'id':name,'original':str(src),'copy':str(dst.relative_to(P)),'sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'kind':'new_frozen_stage'})
for base in [R/'quantum_reservoir_audit_20261002',R/'qf_temperature_confidence_20261002']:
 for src in base.rglob('*'):
  if src.is_file() and src.suffix in ['.md','.csv','.json'] and src.stat().st_size<10_000_000 and 'reference_solver' not in src.parts:
   dst=P/'data/research'/base.name/src.relative_to(base);dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst);manifest.append({'original':str(src),'copy':str(dst.relative_to(P)),'sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'kind':'new_frozen_stage_data'})
(P/'data/source_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
print('Original schematics and frozen additions ready')
