from pathlib import Path
import os,csv
os.environ.setdefault('MPLBACKEND','Agg');os.environ.setdefault('MPLCONFIGDIR','/tmp/trap-free-300K-mpl');os.environ.setdefault('XDG_CACHE_HOME','/tmp/trap-free-300K-cache')
import numpy as np
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parent
read=lambda n:list(csv.DictReader((root/'results'/n).open()))
curves=read('all_evaluated_300K.csv');metrics={x['scenario']:x for x in read('rise_metrics_300K.csv')}
colors={'direct_fixed':'#2878ad','trap_fixed':'#d17b29','direct_a05':'#6b50a1','direct_a10':'#45917b','direct_a20':'#a65769'}
labels={'direct_fixed':'No defect, equal budget','trap_fixed':'Midpoint defect','direct_a05':'No defect, a = 0.5 nm','direct_a10':'No defect, a = 1 nm','direct_a20':'No defect, a = 2 nm'}
plt.rcParams.update({'font.size':10.5,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(2,2,figsize=(12,8.2),constrained_layout=True)
for s in ['direct_fixed','trap_fixed']:
 rr=sorted([x for x in curves if x['scenario']==s and .5<=float(x['F_MVcm'])<=3],key=lambda x:float(x['F_MVcm']))
 F=np.array([float(x['F_MVcm']) for x in rr]);y=np.array([float(x['R_pair_s']) for x in rr]);j=y*1.602176634e-6
 ax[0,0].semilogy(F,j,'-',lw=1.8,c=colors[s],label=labels[s])
 ax[0,1].plot(F,y/float(metrics[s]['R_peak_s']),'-',lw=1.8,c=colors[s],label=labels[s])
 for f,p in [(.1,'10'),(.5,'50'),(.9,'90')]:ax[0,1].plot(float(metrics[s]['F_rise_'+p+'_MVcm']),f,'o',c=colors[s],ms=5)
ax[0,0].axhline(50,c='.35',ls='--',lw=1,label='50 mA/cm$^2$ (not reached)')
ax[0,0].set_title('Equal coupling budget: conditional current')
ax[0,0].set_xlabel('Local field (MV/cm)');ax[0,0].set_ylabel('J for stated pair density (mA/cm$^2$)');ax[0,0].legend(fontsize=9,loc='lower right')
for f in [.1,.5,.9]:ax[0,1].axhline(f,c='.8',lw=.6)
ax[0,1].set_xlim(1.05,1.65);ax[0,1].set_ylim(-.02,1.05);ax[0,1].set_xlabel('Local field (MV/cm)');ax[0,1].set_ylabel('R / own maximum in tested field window');ax[0,1].set_title('First rising-branch 10%, 50%, 90% crossings');ax[0,1].legend(fontsize=9,loc='lower right')
for s in ['trap_fixed','direct_a05','direct_a10','direct_a20']:
 rr=sorted([x for x in curves if x['scenario']==s and .5<=float(x['F_MVcm'])<=3],key=lambda x:float(x['F_MVcm']))
 ax[1,0].semilogy([float(x['F_MVcm']) for x in rr],[float(x['R_pair_s'])*1.602176634e-6 for x in rr],'-',c=colors[s],label=labels[s])
ax[1,0].set_title('Same distance law: no budget restoration');ax[1,0].set_xlabel('Local field (MV/cm)');ax[1,0].set_ylabel('J for stated pair density (mA/cm$^2$)');ax[1,0].legend(fontsize=9,loc='lower right')
ss=['direct_fixed','trap_fixed','direct_a05','direct_a10','direct_a20'];ys=[float(metrics[s]['J_peak_for_Npair_1e15_mAcm2']) for s in ss]
ax[1,1].bar(range(5),ys,color=[colors[s] for s in ss]);ax[1,1].set_yscale('log');ax[1,1].set_ylim(1e-9,100)
ax[1,1].axhline(50,c='.35',ls='--',lw=1);ax[1,1].set_xticks(range(5),['Direct\nbudget','Trap\nbudget','Direct\na=0.5 nm','Direct\na=1 nm','Direct\na=2 nm']);ax[1,1].set_title('Maximum in 0.5-3 MV/cm: amplitude matters');ax[1,1].set_ylabel('Conditional maximum J (mA/cm$^2$)')
for i,y in enumerate(ys):ax[1,1].text(i,y*1.5,f'{y:.3g}',ha='center',va='bottom',fontsize=9)
fig.suptitle('Room-temperature baseline: 300 K only\nNpair = 10$^{15}$ cm$^{-3}$, d = 100 nm; independent available graphs + full collection; not a device fit',fontsize=12)
fig.savefig(root/'figures'/'room_temperature_300K.png',dpi=170)
