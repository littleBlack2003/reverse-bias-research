from pathlib import Path
import json,csv,numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1]
rows=[]
for f in sorted((R/'data').glob('*.json')):
 if not (f.name.startswith('gaussian') or f.name=='refine641.json'):continue
 rows+=json.loads(f.read_text())['rows']
keys=sorted(set().union(*(r.keys() for r in rows)))
with open(R/'data/all_device_points.csv','w') as f:
 w=csv.DictWriter(f,keys);w.writeheader();w.writerows(rows)
lookup={(r['mode'],r['barrier_eV'],r['T'],r['nodes'],r['V']):r for r in rows}
# At 300 K, freeze-at-300 control is exactly the same equations/state.
for V in [0,-5,-10,-15]:
 if ('gaussian',.4,300,321,V) in lookup:lookup[('gaussian_mu300',.4,300,321,V)]=lookup[('gaussian',.4,300,321,V)]
conv=[]
for b in [.4,.5]:
 for T in [300,340,380,420]:
  for V in [-5,-10,-15]:
   if all(('gaussian',b,T,n,V) in lookup for n in [81,161,321]):
    a=[lookup[('gaussian',b,T,n,V)]['J_Acm2'] for n in [81,161,321]]
    conv.append(dict(barrier=b,T=T,V=V,J81=a[0],J161=a[1],J321=a[2],relative_81_161=abs(a[1]/a[0]-1),relative_161_321=abs(a[2]/a[1]-1)))
with open(R/'data/grid_convergence.json','w') as f:json.dump(conv,f,indent=2)
S=R.parent/'S10_arrhenius_digitization_20261003/digitized_currents.csv'
s10=list(csv.DictReader(S.open()));kb=8.617333262145e-5
fits=[]
for mode,b in [('gaussian',.4),('gaussian',.5),('gaussian_constmu',.4),('gaussian_mu300',.4)]:
 for V in [-5,-10,-15]:
  subset=[lookup.get((mode,b,T,321,V)) for T in [300,340,380,420]]
  if any(r is None for r in subset):continue
  ts=np.array([r['T'] for r in subset]);js=np.array([abs(r['J_Acm2']) for r in subset]);coef=np.polyfit(1/ts,np.log(js),1)
  fits.append(dict(mode=mode,barrier=b,V=V,Ea_eV=-kb*coef[0],ratio_420_300=js[-1]/js[0],max_log_residual=float(max(abs(np.log(js)-np.polyval(coef,1/ts)))),J300_Acm2=js[0],J420_Acm2=js[-1]))
for V in [-5,-10,-15]:
 sub=[r for r in s10 if float(r['V_V'])==V];ts=np.array([float(r['T_K']) for r in sub]);js=np.array([float(r['J_mA_cm2']) for r in sub]);coef=np.polyfit(1/ts,np.log(js),1)
 fits.append(dict(mode='Huang_S10_digitized',barrier=None,V=V,Ea_eV=-kb*coef[0],ratio_420_300=js[-1]/js[0],max_log_residual=float(max(abs(np.log(js)-np.polyval(coef,1/ts)))),J300_Acm2=js[0]/1000,J420_Acm2=js[-1]/1000))
(R/'data/arrhenius_fits.json').write_text(json.dumps(fits,indent=2))
fig,ax=plt.subplots(1,3,figsize=(12,4),layout='constrained')
for i,V in enumerate([-5,-10,-15]):
 for mode,b,ls,label in [('gaussian',.4,'o-','EGDM b=0.4 eV'),('gaussian',.5,'s-','EGDM b=0.5 eV'),('gaussian_mu300',.4,'^-','Mobility T argument fixed at 300 K')]:
  sub=[lookup.get((mode,b,T,321,V)) for T in [300,340,380,420]]
  if any(r is None for r in sub):continue
  ax[i].semilogy([1000/r['T'] for r in sub],[abs(r['J_Acm2'])*1000 for r in sub],ls,label=label)
 sub=[r for r in s10 if float(r['V_V'])==V]
 ax[i].errorbar([1000/float(r['T_K']) for r in sub],[float(r['J_mA_cm2']) for r in sub],yerr=[[float(r['J_mA_cm2'])-float(r['J_lower']) for r in sub],[float(r['J_upper'])-float(r['J_mA_cm2']) for r in sub]],fmt='kD--',label='Huang S10 pixels')
 ax[i].set(xlabel='1000/T (1/K)',ylabel='|J| (mA/cm²)',title=f'{V} V terminal bias');ax[i].grid(alpha=.2)
ax[0].legend(fontsize=7);fig.suptitle('Conditional complete devices versus literature (different contact/material assumptions)')
fig.savefig(R/'temperature_device_comparison.png',dpi=180);fig.savefig(R/'temperature_device_comparison.pdf')
fig,ax=plt.subplots(1,3,figsize=(12,4),layout='constrained')
for i,V in enumerate([-5,-10,-15]):
 sub=[lookup.get(('gaussian',.4,T,321,V)) for T in [300,340,380,420]]
 if any(r is None for r in sub):continue
 for field,label in [('J_Acm2','Terminal |J|'),('pair_source_Acm2','Net Marcus pair source'),('bim_Acm2','Net bimolecular source')]:
  ax[i].semilogy([r['T'] for r in sub],[abs(r[field])*1000 for r in sub],'o-',label=label)
 ax[i].set(xlabel='T (K)',ylabel='Integrated current (mA/cm²)',title=f'{V} V; b=0.4 eV');ax[i].grid(alpha=.2)
ax[0].legend(fontsize=8);fig.suptitle('Current source ledger: terminal leakage dominates these contact scenarios')
fig.savefig(R/'source_current_ledger.png',dpi=180);fig.savefig(R/'source_current_ledger.pdf')
print(json.dumps(dict(completed_points=len(rows),convergence_points=len(conv),max_grid_error=max([r['relative_161_321'] for r in conv],default=0),fits=fits),indent=2))
