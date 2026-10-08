from pathlib import Path
import csv,os
os.environ.setdefault('MPLBACKEND','Agg')
os.environ.setdefault('MPLCONFIGDIR','/tmp/vibronic-audit-mpl')
os.environ.setdefault('XDG_CACHE_HOME','/tmp/vibronic-audit-cache')
import numpy as np
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parent
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
read=lambda n:list(csv.DictReader((root/'results'/n).open()))
rows=read('fully_quantum_continuum_comparison.csv');multi=read('matched_lambda_high_mode_comparison.csv')
fig,ax=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
for j,T in enumerate([80.,300.]):
 for ec,c in zip([1,5,20],['#2878ad','#e68a2e','#bb3d3d']):
  r=[x for x in rows if float(x['T_K'])==T and float(x['Ec_meV'])==ec and float(x['delta_G_eV'])<=0]
  ax[0,j].plot([float(x['delta_G_eV']) for x in r],[float(x['quantum_over_existing']) for x in r],'o-',c=c,label=f'$E_c$ = {ec} meV')
 ax[0,j].axhline(1,c='.5',lw=.8);ax[0,j].set_title(f'Quantum slow bath / classical slow bath, {T:g} K')
 ax[0,j].set_xlabel(r'$\Delta G$ (eV)');ax[0,j].set_ylabel('Elementary rate ratio');ax[0,j].legend(fontsize=9)
 r=[x for x in multi if float(x['T_K'])==T and float(x['delta_G_eV'])<=0]
 ax[1,j].plot([float(x['delta_G_eV']) for x in r],[float(x['two_over_single']) for x in r],'o-',color='#7759a7')
 ax[1,j].axhline(1,c='.5',lw=.8);ax[1,j].set_yscale('log');ax[1,j].set_ylim(.18,7)
 ax[1,j].set_title(f'Two high modes / one high mode, {T:g} K\nMatched reorganization and zero-T variance',fontsize=11)
 ax[1,j].set_xlabel(r'$\Delta G$ (eV)');ax[1,j].set_ylabel('Elementary rate ratio (log scale)')
fig.suptitle('Conditional harmonic-bath benchmarks at fixed total reorganization 0.20 eV\nExploratory spectra; no fitted material parameters; lines only connect tested points',fontsize=12)
fig.savefig(root/'figures'/'spectral_sensitivity.png',dpi=170);plt.close(fig)
