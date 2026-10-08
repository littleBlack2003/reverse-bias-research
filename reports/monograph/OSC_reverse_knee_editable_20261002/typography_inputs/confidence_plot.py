from temperature_probe import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'figure.dpi':160})
data=json.loads((ROOT/'data/confidence_map.json').read_text())['rows'];Ts=[80,150,220,300];Vs=[-5,-20]
a=np.array([[next(x['old_independent_pair_relative_error'] for x in data if x['T_K']==T and x['V']==V and x['electron_depth_eV']==.164) for V in Vs] for T in Ts]);b=np.array([[next(x['stable_physical_relative_error'] for x in data if x['T_K']==T and x['V']==V and x['electron_depth_eV']==.164) for V in Vs] for T in Ts])
fig,axes=plt.subplots(1,2,figsize=(10.8,4.7),constrained_layout=True)
for ax,m,title in zip(axes,[a,b],['Old source: independent cycle mismatch','Stable source: full-ledger error']):
 im=ax.imshow(np.log10(np.maximum(m,1e-16)),vmin=-16,vmax=22,cmap='coolwarm',aspect='auto')
 ax.set_xticks([0,1],['−5 V','−20 V']);ax.set_yticks(range(4),[f'{T} K' for T in Ts]);ax.set_title(title,fontsize=11)
 for i in range(4):
  for j in range(2):ax.text(j,i,f'{m[i,j]:.2e}',ha='center',va='center',color='white' if np.log10(max(m[i,j],1e-16))< -8 or np.log10(max(m[i,j],1e-16))>15 else 'black',fontsize=11)
fig.colorbar(im,ax=axes,label='log10(relative error)',shrink=.84)
fig.suptitle('Numerical stress test | Dn/Dp = 0.164/1.136 eV | classical rates held fixed\nNumerical precision is not physical validity at low temperature',fontsize=12)
fig.savefig(ROOT/'figures/cancellation_confidence.png');plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(10.5,4.1),constrained_layout=True)
for ax,D in zip(axes,[.65,.164]):
 for V,c in zip(Vs,['#2864a2','#cf722d']):
  vals=[-next(x['stable_current_Acm2'] for x in data if x['T_K']==T and x['V']==V and x['electron_depth_eV']==D) for T in Ts]
  ax.semilogy(Ts,vals,'o-',label=f'{V} V',color=c)
 ax.set(xlabel='Temperature (K)',ylabel='Formal dark |J| (A/cm²)',title=f'Dn/Dp = {D:.3f}/{1.3-D:.3f} eV');ax.grid(alpha=.2);ax.legend()
fig.suptitle('Fixed-parameter classical model, not a material temperature prediction',fontsize=12)
fig.savefig(ROOT/'figures/formal_classical_temperature.png');plt.close(fig)
