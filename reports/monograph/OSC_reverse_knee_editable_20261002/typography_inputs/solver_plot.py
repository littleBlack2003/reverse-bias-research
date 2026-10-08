from solver_probe import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'figure.dpi':160,'axes.grid':True,'grid.alpha':.2})
a=json.loads((ROOT/'data/original_failed_newton.json').read_text())['history'];it=[x['iteration'] for x in a]
fig,ax=plt.subplots(1,2,figsize=(10,3.8),constrained_layout=True)
ax[0].semilogy(it,[x['step_max'] for x in a],'o-');ax[0].set(xlabel='Newton iteration',ylabel='Maximum raw log-density step',title='The depleted hole variable runs away')
ax[1].semilogy(it,[x['trials'][0][0] for x in a],'o-');ax[1].set(xlabel='Newton iteration',ylabel='Initial global damping factor',title='All variables become unable to advance')
fig.suptitle('Preserved failure at −8.4 V | Dn = 0.164 eV, Dp = 1.136 eV')
fig.savefig(ROOT/'figures/newton_failure.png');plt.close(fig)
fig,ax=plt.subplots(1,2,figsize=(10,3.8),constrained_layout=True)
for case in ['qf_asym','qf_asym_161']:
 d=json.loads((ROOT/'data'/f'{case}.json').read_text());r=d['rows'][1:];ax[0].semilogy([-x['V'] for x in r],[-x['J_Acm2'] for x in r],'o-',label=f"{d['nodes']} nodes")
ax[0].set(xlabel='Reverse bias magnitude (V)',ylabel='Dark |J| (A/cm²)',title='Recovered asymmetric branch');ax[0].legend()
v=json.loads((ROOT/'data/qf_validation.json').read_text())['physical_gate_on_old_false_pass'];labels=['Scaled residual only','With current budget'];ax[1].bar(labels,[v['no_gate']['error_Acm2'],v['gate']['error_Acm2']],color=['#be553b','#288b79']);ax[1].set_yscale('log');ax[1].axhline(v['gate']['budget_Acm2'],color='k',ls='--',label='Declared error budget');ax[1].set(ylabel='Maximum physical current error (A/cm²)',title='Same saved state at −5 V');ax[1].legend(fontsize=8)
fig.suptitle('Equivalent equations | contact-referenced quasi-Fermi variables | 300 K')
fig.savefig(ROOT/'figures/recovered_branch_and_gate.png');plt.close(fig)
