from pathlib import Path
import os,csv
os.environ.setdefault('MPLBACKEND','Agg');os.environ.setdefault('MPLCONFIGDIR','/tmp/trap-free-mpl');os.environ.setdefault('XDG_CACHE_HOME','/tmp/trap-free-cache')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parent
read=lambda n:list(csv.DictReader((root/'results'/n).open()))
a=read('fixed_budget_comparison.csv');d=read('distance_law_comparison.csv')
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,3,figsize=(14,4.6),constrained_layout=True)
for ax,T in zip(axes[:2],[80,300]):
 for topology,label,color in [('direct','No defect: one 10 nm transfer','#2878ad'),('trap','One midpoint defect: two 5 nm transfers','#d17b29')]:
  rr=[r for r in a if int(r['T_K'])==T and r['topology']==topology and float(r['F_MVcm'])>=1]
  ax.semilogy([float(r['F_MVcm']) for r in rr],[float(r['R_pair_s']) for r in rr],'o-',label=label,color=color)
 ax.set_xlabel('Local field (MV/cm)');ax.set_ylabel('Net pair turnover per available graph (s$^{-1}$)');ax.set_title(f'Equal internal coupling budget, {T} K');ax.legend(fontsize=8,loc='lower left' if T==300 else 'lower right')
for topology,label,color in [('direct','No defect, distance cost','#2878ad'),('trap','Midpoint defect','#d17b29')]:
 rr=[r for r in d if int(r['T_K'])==80 and float(r['F_MVcm'])==1.5 and r['topology']==topology]
 axes[2].semilogy([float(r['decay_nm']) for r in rr],[float(r['R_pair_s']) for r in rr],'o-',label=label,color=color)
axes[2].set_xlabel('Amplitude decay length a (nm)');axes[2].set_ylabel('Net pair turnover (s$^{-1}$)');axes[2].set_title('Same H(r) law, 80 K, 1.5 MV/cm');axes[2].legend(fontsize=9)
fig.suptitle('Trap-free versus trap-assisted: conditional dark, isothermal graphs\nSame 10 nm endpoints and finite reservoirs; exploratory parameters; lines connect tested points',fontsize=12)
fig.savefig(root/'figures'/'trap_free_comparison.png',dpi=180)
