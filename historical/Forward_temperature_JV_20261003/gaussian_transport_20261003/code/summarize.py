import json,csv,numpy as np
from pathlib import Path
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];T=np.array([300,340,380,420]);modes=['legacy','mobility_only','gaussian_constmu','gaussian'];volts=np.array([-5,-10,-15])
allrows=[];curves={};summaries={}
for mode in modes:
 rr=[json.loads((ROOT/'data'/f'{mode}_T{t}_N81.json').read_text())['rows'] for t in T]
 arr=np.array([[abs(r['J_Acm2']) for r in rows[1:]] for rows in rr]);curves[mode]=arr
 summaries[mode]=dict(Ea_eV=(-8.617333262145e-5*np.polyfit(1/T,np.log(arr),1)[0]).tolist(),ratios420_300=(arr[-1]/arr[0]).tolist(),J_Acm2=arr.tolist())
 for rs in rr:allrows.extend(rs)
SI=ROOT/'reference/S10';raw=list(csv.DictReader((SI/'digitized_currents.csv').open()));meas=np.array([[float(next(r['J_mA_cm2'] for r in raw if float(r['T_K'])==t and float(r['V_V'])==v))/1000 for v in volts] for t in T]);fits=json.loads((SI/'analysis.json').read_text())['fits']
fig,axs=plt.subplots(2,2,figsize=(11,8))
colors=['#0072B2','#D55E00','#009E73']
for i,(v,c) in enumerate(zip(volts,colors)):
 axs[0,0].semilogy(T,meas[:,i]/meas[0,i],'o',color=c,label=f'S10 {v} V')
 axs[0,0].semilogy(T,curves['gaussian'][:,i]/curves['gaussian'][0,i],'-',color=c)
 axs[0,0].semilogy(T,curves['gaussian_constmu'][:,i]/curves['gaussian_constmu'][0,i],'--',color=c)
axs[0,0].set(xlabel='Temperature (K)',ylabel='|J(T)| / |J(300 K)|',title='Dots: digitized S10; solid: Gaussian; dashed: frozen μ')
axs[0,0].legend(fontsize=8)
for mode in modes:axs[0,1].plot(-volts,summaries[mode]['Ea_eV'],'o-',label=mode)
axs[0,1].errorbar(-volts,[r['Ea_eV'] for r in fits],yerr=[r['digitization_bound_eV'] for r in fits],fmt='ks',label='S10 pixel envelope')
axs[0,1].set(xlabel='Reverse voltage magnitude (V)',ylabel='Apparent activation energy (eV)',title='Four-temperature fit, not microscopic barrier')
axs[0,1].legend(fontsize=8)
for mode in ['legacy','gaussian']:axs[1,0].semilogy(-volts,curves[mode][0]*1000,'o-',label=mode)
axs[1,0].semilogy(-volts,meas[0]*1000,'ks--',label='S10 at 300 K')
axs[1,0].set(xlabel='Reverse voltage magnitude (V)',ylabel='|J| (mA/cm²)',title='Absolute-current mismatch remains')
axs[1,0].legend(fontsize=8)
for i,(v,c) in enumerate(zip(volts,colors)):
 axs[1,1].plot(T,curves['gaussian'][:,i]/curves['gaussian_constmu'][:,i],'-o',color=c,label=f'{v} V')
axs[1,1].set(xlabel='Temperature (K)',ylabel='J(full μ) / J(frozen μ)',title='Same Gaussian EOS, rates, and contacts')
axs[1,1].legend(fontsize=8)
for ax in axs.flat:ax.grid(alpha=.2)
fig.suptitle('UNCONVERGED Gaussian device curves: diagnostic only, NOT a literature fit',fontsize=12,color='red')
fig.tight_layout();fig.savefig(ROOT/'data/temperature_comparison.png',dpi=180)
mesh={}
for temp in [300,420]:
 r={n:json.loads((ROOT/'data'/f'gaussian_T{temp}_N{n}.json').read_text())['rows'] for n in [81,161,321] if (ROOT/'data'/f'gaussian_T{temp}_N{n}.json').exists()}
 if 321 in r:mesh[temp]={str(v):{str(n):rs[i+1]['J_Acm2'] for n,rs in r.items() if len(rs)>i+1} for i,v in enumerate(volts)}
(ROOT/'data/summary.json').write_text(json.dumps(dict(temperature=T.tolist(),V=volts.tolist(),models=summaries,mesh=mesh,S10_apparent_fits=fits,completed_state_count=len(allrows)),indent=2))
print(json.dumps(dict(models=summaries,mesh=mesh),indent=2))
