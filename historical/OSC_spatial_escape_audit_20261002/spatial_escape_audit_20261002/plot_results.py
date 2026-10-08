from pathlib import Path
import csv,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).parent
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':150})
def load(fn):return list(csv.DictReader((R/'results'/fn).open()))
def vals(rows,k):return np.array([float(r[k]) for r in rows])
def select(rows,**kw):return [r for r in rows if all(str(r[k])==str(v) for k,v in kw.items())]
cs=load('controlled_cases.csv');cols={'colocated':'#777777','direct':'#dd8032','explicit':'#126f94'}
fig,ax=plt.subplots(2,2,figsize=(10.5,6.7),gridspec_kw={'height_ratios':[.72,1.]})
for panel,(v,title) in zip(ax[0],[('direct','Direct remote bath: no tagged outer charge'),('explicit','Explicit host escape: retain outer charges')]):
 panel.axhline(0,color='#cccccc',lw=2,zorder=0)
 xs=[-5,0,5] if v=='direct' else [-6,-5,0,5,6];labels=['V','D','C'] if v=='direct' else ['Vo','V','D','C','Co']
 panel.scatter(xs,[0]*len(xs),s=110,c=cols[v],zorder=3)
 for x,lab in zip(xs,labels):panel.text(x,.24,lab,ha='center',fontsize=11)
 for x in [-6,6]:panel.plot([x,x],[-.4,.5],color='#333333',lw=3);panel.text(x,-.57,f'{x} nm\nbath',ha='center',fontsize=9)
 panel.set(xlim=(-7,7),ylim=(-1,1),title=title);panel.axis('off')
for a,T in zip(ax[1],[80,300]):
 for v in cols:
  z=select(cs,T_K=T,variant=v);a.semilogy(vals(z,'F_MVcm'),vals(z,'R_s'),'-o',label=v+' (10 nm reference)' if v=='colocated' else v+' (12 nm)',color=cols[v],ms=4)
 a.set(xlabel='F (MV/cm)',ylabel='Net cycle rate (1/s per graph)',title=f'{T} K; fixed total H-squared budget');a.grid(alpha=.2);a.legend(fontsize=8)
fig.suptitle('Spatial escape audit | exploratory isothermal model',fontsize=14);fig.text(.5,.005,'Same electrochemical drive delta-mu = F eV; 12 nm controls have electrostatic drop 1.2 F V. Lines guide the eye.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.025,1,.96]);fig.savefig(R/'figures/geometry_and_controlled_rates.png',dpi=220);plt.close(fig)

fig,ax=plt.subplots(1,2,figsize=(10.5,3.6));tags=load('last_tag_placement.csv')
for T,c in [(80,'#126f94'),(300,'#dd8032')]:
 z=select(tags,T_K=T,F_MVcm='1.5');ax[0].semilogy(vals(z,'s_nm'),vals(z,'R_s'),'-o',label=f'{T} K',color=c)
ax[0].axvline(5.761993,ls=':',c='gray',label='Pair export gap = 0')
ax[0].set(xlabel='Last tagged-site position |x| = s (nm)',ylabel='Net rate (1/s per graph)',title='Baths fixed at +/-6 nm; F = 1.5 MV/cm');ax[0].legend(fontsize=8);ax[0].grid(alpha=.2)
occ=load('multiple_occupancy_controls.csv')
for T,c in [(80,'#126f94'),(300,'#dd8032')]:
 for cap,ls in [('False','-'),('True','--')]:
  z=select(occ,T_K=T,single_pair_constraint=cap);ax[1].semilogy(vals(z,'F_MVcm'),vals(z,'R_s'),ls,c=c,label=f'{T} K, '+('32 states' if cap=='False' else 'restricted 12 states'))
ax[1].set(xlabel='F (MV/cm)',ylabel='Net rate (1/s per graph)',title='High-field bypass uses multi-charge states');ax[1].legend(fontsize=8);ax[1].grid(alpha=.2)
fig.text(.5,.008,'Left: changes correlation cutoff/geometry, not numerical convergence. Right: restriction changes the physical occupancy model.',ha='center',fontsize=8)
fig.tight_layout(rect=[0,.04,1,1]);fig.savefig(R/'figures/cutoff_and_multicharge.png',dpi=220);plt.close(fig)

fig,ax=plt.subplots(1,2,figsize=(10.5,3.6));bud=load('fixed_budget_allocation.csv');de=load('distance_attenuation_controls.csv')
for T,c in [(80,'#126f94'),(300,'#dd8032')]:
 z=select(bud,T_K=T);ax[0].loglog(vals(z,'eta'),vals(z,'R_s'),'-o',color=c,label=f'{T} K')
ax[0].set(xlabel='Budget fraction eta for escape hop',ylabel='Net rate (1/s per graph)',title='Fixed H-squared: hop vs final exchange');ax[0].legend();ax[0].grid(alpha=.2)
for v,c in [('direct','#dd8032'),('explicit','#126f94')]:
 z=select(de,T_K=80,variant=v);ax[1].semilogy(vals(z,'decay_nm'),vals(z,'R_s'),'-o',color=c,label=v)
ax[1].set(xlabel='Assumed coupling decay length xi (nm)',ylabel='Net rate (1/s per graph)',title='80 K: distance attenuation, unequal actual budgets');ax[1].legend();ax[1].grid(alpha=.2)
fig.text(.5,.008,'F = 1.5 MV/cm. Attenuation H -> H exp(-distance/xi) is an uncalibrated sensitivity test; its prefactor is not renormalized.',ha='center',fontsize=8)
fig.tight_layout(rect=[0,.04,1,1]);fig.savefig(R/'figures/budget_and_distance.png',dpi=220);plt.close(fig)

fig,ax=plt.subplots(1,2,figsize=(10.5,3.5));z=select(cs,T_K=80,variant='explicit');F=vals(z,'F_MVcm');rate=vals(z,'R_s')
for k,c,l in [('phonon_heat_eV_s','#126f94','Phonon baths'),('electron_bath_heat_eV_s','#dd8032','Electron baths'),('total_bath_heat_eV_s','#333333','Total = delta-mu')]:ax[0].plot(F,vals(z,k)/rate,'-o',c=c,label=l)
ax[0].axhline(0,c='gray',lw=.8);ax[0].set(xlabel='F (MV/cm)',ylabel='Heat per net cycle (eV)',title='Explicit 32-state, 80 K: heat channels');ax[0].legend(fontsize=8);ax[0].grid(alpha=.2)
for v in ['direct','explicit']:
 z=select(cs,T_K=80,variant=v);w=select(cs,T_K=300,variant=v);ax[1].semilogy(vals(z,'F_MVcm'),vals(z,'R_s')/vals(w,'R_s'),'-o',c=cols[v],label=v)
ax[1].axhline(1,c='gray',ls=':');ax[1].set(xlabel='F (MV/cm)',ylabel='R(80 K) / R(300 K)',title='A local bottleneck is not a universal temperature law');ax[1].legend();ax[1].grid(alpha=.2)
fig.tight_layout();fig.savefig(R/'figures/heat_and_temperature.png',dpi=220);plt.close(fig)
