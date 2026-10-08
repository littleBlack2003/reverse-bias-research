from pathlib import Path
import csv,json,hashlib,os,shutil
P=Path(__file__).resolve().parents[1];os.environ['MPLCONFIGDIR']=str(P/'qa/mpl')
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import fontManager,FontProperties
import numpy as np
fontManager.addfont(str(P/'assets/NotoSerifSC-Regular.ttf'));fontManager.addfont(str(P/'assets/LiberationSerif-Regular.ttf'));fontManager.addfont(str(P/'assets/NotoSerif-SymbolFallback.ttf'));fp=FontProperties(fname=str(P/'assets/NotoSerifSC-Regular.ttf'))
plt.rcParams.update({'font.family':['Liberation Serif','Noto Serif SC','Noto Serif'],'mathtext.fontset':'stix','font.size':11.5,'axes.labelsize':11.5,'axes.titlesize':12,'legend.fontsize':10,'xtick.labelsize':10,'ytick.labelsize':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'path'})
M=P/'data/002高YZ/OSC_FF_reverse_knee_model/model_comparison/results_ff78';R=P/'data/research';C=['#215B80','#087F8C','#BD6A27','#775899','#B33C41','#65717B','#7C852B']
entries=[]
def read(path):
 entries.append({'original':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()});return list(csv.DictReader(path.open(encoding='utf-8-sig')))
def f(x,k):return float(x[k])
def save(fig,name):
 orig=P/'figures/original_archive';orig.mkdir(exist_ok=True)
 p=P/'figures'/f'{name}.png'
 if p.exists() and not (orig/p.name).exists():shutil.copy2(p,orig/p.name)
 fig.tight_layout(pad=.7);fig.savefig(p,dpi=250,bbox_inches='tight',facecolor='white');fig.savefig(p.with_suffix('.svg'),bbox_inches='tight',facecolor='white');plt.close(fig)
# One readable chart per page for the family overview, fixed default branch values.
rows=read(M/'comparison_curves.csv');defs={'local':(0,'局域SRH'),'pf':(1e15,'有限PF'),'hopping':(1e15,'MA跳跃'),'contact':(.6,'接触'),'thermal':(.0002,'电热'),'weak':(.01,'弱区'),'empirical':(.05,'经验')}
fig,ax=plt.subplots(figsize=(7,3.5))
for i,(m,(value,label)) in enumerate(defs.items()):
 a=[x for x in rows if x['model']==m and x['group']=='family' and f(x,'d_nm')==100 and f(x,'bath_K')==300 and f(x,'parameter_value')==value and f(x,'light')==1 and f(x,'V')<=0];a.sort(key=lambda x:abs(f(x,'V')))
 ax.semilogy([abs(f(x,'V')) for x in a],[abs(f(x,'J'))*1000 for x in a],label=label,color=C[i])
ax.axhline(50,color='#444',ls='--',lw=.8);ax.set(xlabel='反偏电压幅值 / V',ylabel='光照电流幅值 / mA cm⁻²',ylim=(15,3000),xlim=(0,30));ax.legend(ncol=4,loc='upper left');save(fig,'overview')
# Local A-D in photovoltaic window and reverse illuminated region; retain full data in source package.
rows=read(M/'four_rate_diagnostic/curves.csv')
fig,axs=plt.subplots(1,2,figsize=(7,3.3))
for i,m in enumerate('ABCD'):
 a=[x for x in rows if x['model']==m and f(x,'light')==1];a.sort(key=lambda x:f(x,'V'))
 pv=[x for x in a if 0<=f(x,'V')<=1.025];rv=[x for x in a if f(x,'V')<=0];rv.sort(key=lambda x:-f(x,'V'))
 axs[0].plot([f(x,'V') for x in pv],[-f(x,'J')*1000 for x in pv],color=C[i],label=m,lw=1.6,ls='--' if m=='D' else '-')
 axs[1].semilogy([-f(x,'V') for x in rv],[abs(f(x,'J'))*1000 for x in rv],color=C[i],label=m,lw=1.6,ls='--' if m=='D' else '-')
axs[0].set(xlabel='正向电压 / V',ylabel='输出电流 / mA cm⁻²',ylim=(-2,26),title='光伏工作区');axs[1].set(xlabel='反偏幅值 / V',ylabel='光照幅值 / mA cm⁻²',title='反向光照');axs[0].legend(ncol=2);save(fig,'abcd')
# A/B plus finite PF from exact stored table.
fig,axs=plt.subplots(1,2,figsize=(7,3.3))
for i,(m,lab) in enumerate([('A','局域SRH'),('B','局域全倍率')]):
 a=[x for x in rows if x['model']==m and f(x,'light')==1];pv=sorted([x for x in a if 0<=f(x,'V')<=1.025],key=lambda x:f(x,'V'));rv=sorted([x for x in a if f(x,'V')<=0],key=lambda x:-f(x,'V'))
 axs[0].plot([f(x,'V') for x in pv],[-f(x,'J')*1000 for x in pv],color=C[i],label=lab);axs[1].semilogy([-f(x,'V') for x in rv],[abs(f(x,'J'))*1000 for x in rv],color=C[i],label=lab)
allr=read(M/'comparison_curves.csv');a=[x for x in allr if x['model']=='pf' and x['group']=='family' and f(x,'d_nm')==100 and f(x,'parameter_value')==1e15 and f(x,'light')==1];pv=sorted([x for x in a if 0<=f(x,'V')<=1.025],key=lambda x:f(x,'V'));rv=sorted([x for x in a if f(x,'V')<=0],key=lambda x:-f(x,'V'))
axs[0].plot([f(x,'V') for x in pv],[-f(x,'J')*1000 for x in pv],color=C[3],label='有限PF');axs[1].semilogy([-f(x,'V') for x in rv],[abs(f(x,'J'))*1000 for x in rv],color=C[3],label='有限PF')
axs[0].set(xlabel='正向电压 / V',ylabel='输出电流 / mA cm⁻²',ylim=(-2,26));axs[1].set(xlabel='反偏幅值 / V',ylabel='光照幅值 / mA cm⁻²',ylim=(1e-4,1e4));axs[0].legend(fontsize=9);save(fig,'local_pf')
rows=read(M/'pf_density_scan/sequence.csv');fig,axs=plt.subplots(1,2,figsize=(7,3.2));axs[0].semilogx([f(x,'Nt_cm3') for x in rows],[f(x,'FF_percent') for x in rows],'o-',ms=3,color=C[0]);a=[x for x in rows if x['V50_V']];axs[1].semilogx([f(x,'Nt_cm3') for x in a],[-f(x,'V50_V') for x in a],'o-',ms=3,color=C[2]);axs[0].set(xlabel='有效 Nₜ / cm⁻³',ylabel='FF / %');axs[1].set(xlabel='有效 Nₜ / cm⁻³',ylabel='$|V_{50}|$ / V');save(fig,'density')
rows=read(M/'pf_thickness_scan/thickness_fits.csv');fig,ax=plt.subplots(figsize=(7,3.0))
for i,th in enumerate([40,50,75]):
 a=sorted([x for x in rows if f(x,'threshold_mAcm2')==th and x['status']=='full_range'],key=lambda x:f(x,'FF100_percent'));ax.plot([f(x,'FF100_percent') for x in a],[f(x,'K_V_per_nm') for x in a],'o-',ms=3,label=f'{th} mA/cm²阈值',color=C[i])
ax.set(xlabel='100 nm器件FF / %',ylabel='五厚度拟合 K / V nm⁻¹');ax.legend();save(fig,'thickness')
rows=read(M/'tat_results/curves.csv');fig,axs=plt.subplots(1,2,figsize=(7,3.3))
for light,c in [(0,0),(1,1)]:
 a=[x for x in rows if f(x,'d_nm')==100 and f(x,'T')==300 and f(x,'light')==light];rv=sorted([x for x in a if f(x,'V')<0],key=lambda x:-f(x,'V'));axs[1].plot([-f(x,'V') for x in rv],[abs(f(x,'J'))*1000 for x in rv],label='暗态' if not light else '光照',color=C[c])
 if light:
  pv=sorted([x for x in a if 0<=f(x,'V')<=1.025],key=lambda x:f(x,'V'));axs[0].plot([f(x,'V') for x in pv],[-f(x,'J')*1000 for x in pv],label='总电流',color=C[1]);axs[0].plot([f(x,'V') for x in pv],[-f(x,'J_baseline')*1000 for x in pv],label='DD基线',color=C[0],ls='--')
axs[0].set(xlabel='正向电压 / V',ylabel='输出电流 / mA cm⁻²',ylim=(-5,26));axs[1].set(xlabel='反偏幅值 / V',ylabel='总电流幅值 / mA cm⁻²');axs[0].legend();axs[1].legend();save(fig,'tat')
rows=read(M/'thermal_rate_curves.csv');fig,axs=plt.subplots(1,2,figsize=(7,3.2))
for i,rate in enumerate(sorted(set(f(x,'rate_Vs') for x in rows))):
 a=sorted([x for x in rows if f(x,'rate_Vs')==rate],key=lambda x:-f(x,'V'));axs[0].semilogy([-f(x,'V') for x in a],[-f(x,'J')*1000 for x in a],color=C[i%7],label=f'{rate:g} V/s');axs[1].plot([-f(x,'V') for x in a],[f(x,'T') for x in a],color=C[i%7])
axs[0].set(xlabel='反偏幅值 / V',ylabel='光照幅值 / mA cm⁻²',ylim=(20,1e4));axs[1].set(xlabel='反偏幅值 / V',ylabel='集总温度 / K');axs[0].legend(fontsize=9);save(fig,'thermal')
# Capacity vs localization length is the operative comparison in the associated chapter.
rows=read(R/'temperature_mechanism_audit_20261002/results/coupling_capacity_bounds.csv');fig,ax=plt.subplots(figsize=(7,3.0))
for i,(T,k,lab) in enumerate([(300,'classical','300 K经典'),(80,'classical','80 K经典')]):
 a=sorted([x for x in rows if f(x,'T_K')==T and x['kernel']==k],key=lambda x:f(x,'amplitude_decay_length_nm'));ax.semilogy([f(x,'amplitude_decay_length_nm') for x in a],[f(x,'symmetric_cycle_upper_J_Acm2')*1000 for x in a],'o-',color=C[i],label=lab)
ax.axhline(50,color=C[2],ls='--',label='50 mA/cm²');ax.set(xlabel='振幅衰减长度 a / nm',ylabel='理想两步上界 / mA cm⁻²');ax.legend();save(fig,'pair_capacity')
(P/'qa/redrawn_reproduction_provenance.json').write_text(json.dumps({'scope':'redrawn from unchanged numeric CSV, no rerun or fit','inputs':entries,'outputs':['overview','abcd','local_pf','density','thickness','tat','thermal','pair_capacity'],'display_notes':'Photovoltaic panels show0–1.025V and declared y windows; overview y15–3000mA/cm2; full curves retained in data. Source plots kept in figures/original_archive.'},ensure_ascii=False,indent=2))
print('8 result figures redrawn with readable labels')
