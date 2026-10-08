from rate_audit import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'figure.dpi':160,'axes.grid':True,'grid.alpha':.2})
fig,ax=plt.subplots(1,2,figsize=(10.8,4.2),constrained_layout=True)
for case,label in [('baseline_81','Both branches'),('n_only','Electron PF only'),('p_only','Hole PF only'),('neither','Neither branch')]:
 d=json.loads((ROOT/'data'/f'{case}.json').read_text());r=d['rows'][1:]
 ax[0].semilogy([-x['V'] for x in r],[-x['J_mAcm2'] for x in r],'o-',label=label)
ax[0].set(xlabel='Reverse bias magnitude (V)',ylabel='Dark |J| (mA/cm²)',title='Branch control with all inverse rates matched');ax[0].legend(fontsize=8)
for case in ['baseline_81','baseline_161','baseline_321']:
 d=json.loads((ROOT/'data'/f'{case}.json').read_text());r=d['rows'][1:]
 ax[1].plot([-x['V'] for x in r],[-x['J_mAcm2'] for x in r],'o-',label=case.replace('baseline_','N = '))
ax[1].set(xlabel='Reverse bias magnitude (V)',ylabel='Dark |J| (mA/cm²)',title='Mesh check: fixed physical 5 nm jumps');ax[1].legend()
fig.suptitle('FF78 constructed baseline | 100 nm | 300 K | Nt = 10¹⁵ cm⁻³ | hypothesis model')
fig.savefig(ROOT/'figures/branch_and_mesh.png');plt.close(fig)
d=device(161);a=np.load(ROOT/'data/baseline_161_-20.npz');x=d.x*100;xf=(x[1:]+x[:-1])/2
fig,ax=plt.subplots(1,2,figsize=(10.8,4.2),constrained_layout=True)
ax[0].plot(xf,(a['Jn']+a['Jp'])*1000,label='DD current');ax[0].plot(xf,a['Jnl_independent']*1000,label='Nonlocal displacement');ax[0].plot(xf,(a['Jn']+a['Jp']+a['Jnl_independent'])*1000,label='Total');ax[0].set(xlabel='Position (nm)',ylabel='J (mA/cm²)',title='Event displacement closes spatial current');ax[0].legend(fontsize=8)
ax[1].semilogy(x[1:-1],a['pair_gen'][1:-1],label='Pair generation per trap');ax[1].semilogy(x[1:-1],a['pair_rec'][1:-1],label='Pair recombination per trap');ax[1].semilogy(x[1:-1],a['En'][1:-1],label='Electron emission sum');ax[1].semilogy(x[1:-1],a['Ep'][1:-1],label='Hole emission sum');ax[1].set(xlabel='Position (nm)',ylabel='Rate (s⁻¹)',title='Generation and capture-emission sources');ax[1].legend(fontsize=8)
fig.suptitle('Fresh solution | −20 V dark | 161 nodes | exploratory PF model')
fig.savefig(ROOT/'figures/spatial_balance.png');plt.close(fig)
