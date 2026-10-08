"""Reproduce a bounded 300 K synthesis from archived CSVs; no model solves."""
from pathlib import Path
import os,json,hashlib
os.environ.setdefault('MPLCONFIGDIR','/tmp/mpl-300k-synthesis')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent.parent
OUT=Path(__file__).resolve().parent
sources={
 'metrics':'trap_free_comparator_20261002/room_temperature_300K/results/rise_metrics_300K.csv',
 'curves':'trap_free_comparator_20261002/room_temperature_300K/results/all_evaluated_300K.csv',
 'density':'density_knee_disentangling_20261002/data/selfconsistent_curves.csv',
 'crossings':'density_knee_disentangling_20261002/data/current_crossings.csv',
 'spatial':'spatial_escape_audit_20261002/results/controlled_cases.csv',
 'vibronic':'vibronic_validity_audit_20261002/results/matched_lambda_high_mode_comparison.csv',
 'thermal':'vibronic_validity_audit_20261002/results/prior_graph_thermalization_scales.csv'}
d={k:pd.read_csv(ROOT/v) for k,v in sources.items()}
m=d['metrics']; c=d['curves']; n=d['density']; x=d['crossings']
assert (m.T_K==300).all() and (c.T_K==300).all()
assert np.allclose(m.J_peak_for_Npair_1e15_mAcm2,1.602176634e-19*1e15*1e-5*1e3*m.R_peak_s,rtol=1e-12)
assert np.allclose(m.required_Npair_for_50mA_cm3,50/(1.602176634e-19*1e-5*1e3*m.R_peak_s),rtol=1e-12)
m.to_csv(OUT/'graph_metrics_300K.csv',index=False)
c[['scenario','T_K','F_MVcm','R_pair_s','conditional_J_for_Npair_1e15_and_d100nm_Acm2','internal_total_H2_eV2']].to_csv(OUT/'graph_curves_300K.csv',index=False)
nd=n[(n.model=='PF')&(n.light==0)&(n.nodes==81)].copy()
assert not nd.duplicated(['Nt_cm3','U']).any()
nd[['Nt_cm3','U','I_Acm2','pair_per_Nt','normalized_pair_vs_Nt15']].to_csv(OUT/'density_dark_300K.csv',index=False)
x.to_csv(OUT/'fixed_current_crossings_300K.csv',index=False)
for key in ['spatial','vibronic','thermal']:
 d[key][d[key].T_K==300].to_csv(OUT/(key+'_300K.csv'),index=False)
hi=n[(n.model=='PF')&(n.U>=10)&(n.U<=20)&(n.nodes==81)]
err=(hi.normalized_pair_vs_Nt15-1).abs().max()
u=x[(x.light==0)&(x.current_mAcm2==50)].sort_values('Nt_cm3')
checks={'temperature_K':300,'new_model_solves':0,'density_dark_U50_shift_V':float(u.U.iloc[0]-u.U.iloc[-1]),'max_density_normalized_source_deviation_10to20V':float(err),'max_graph_conditional_J_mAcm2':float(m.J_peak_for_Npair_1e15_mAcm2.max()),'direct_equal_vs_a05_peak_rate_ratio':float(m.loc[m.scenario=='direct_fixed','R_peak_s'].iloc[0]/m.loc[m.scenario=='direct_a05','R_peak_s'].iloc[0]),'all_graph_50mA_thresholds_unreached_at_Npair1e15':bool((m.J_peak_for_Npair_1e15_mAcm2<50).all())}
(OUT/'arithmetic_checks.json').write_text(json.dumps(checks,indent=2))
(OUT/'source_provenance.json').write_text(json.dumps({k:{'path':v,'sha256':hashlib.sha256((ROOT/v).read_bytes()).hexdigest()} for k,v in sources.items()},indent=2))
labels={'direct_fixed':'Direct, equal budget','trap_fixed':'Midpoint trap, equal budget','direct_a05':'Direct, a = 0.5 nm','direct_a10':'Direct, a = 1 nm','direct_a20':'Direct, a = 2 nm'}
colors=['#1565c0','#d35d14','#863eb0','#25865d','#60656d']
fig,axs=plt.subplots(2,2,figsize=(12,8.8),layout='constrained')
for scenario,col in zip(labels,colors):
 f=c[c.scenario==scenario].sort_values('F_MVcm'); peak=m[m.scenario==scenario].R_peak_s.iloc[0]
 axs[0,0].plot(f.F_MVcm,f.R_pair_s/peak,color=col,label=labels[scenario],lw=1.7)
 axs[0,1].semilogy(f.F_MVcm,1e3*f.conditional_J_for_Npair_1e15_and_d100nm_Acm2,color=col,lw=1.7)
axs[0,0].set(xlim=(1,1.8),ylim=(0,1.07),xlabel='Local field F (MV/cm)',ylabel='R / own window peak',title='A  Relative rise: graph turnover, not bare hop rate')
axs[0,0].legend(fontsize=8,loc='lower right')
axs[0,1].axhline(50,color='black',ls='--',lw=1);axs[0,1].text(.56,75,'50 mA/cm²: no crossing',fontsize=9)
axs[0,1].set(xlim=(.5,3),ylim=(1e-10,500),xlabel='Local field F (MV/cm)',ylabel='Conditional J (mA/cm²)',title='B  Same density does not guarantee useful amplitude')
axs[0,1].text(.55,1.5,'Npair = 10¹⁵ cm⁻³; d = 100 nm\nIndependent graphs; complete device collection',fontsize=8)
for nt,col in zip(sorted(nd.Nt_cm3.unique()),['#25865d','#1565c0','#d35d14']):
 f=nd[nd.Nt_cm3==nt].sort_values('U'); name=f'Nt = 10^{int(np.log10(nt))} cm⁻³'
 axs[1,0].semilogy(f.U,1e3*f.I_Acm2,color=col,label=name)
 h=f[(f.U>=10)&(f.U<=20)]
 axs[1,1].plot(h.U,100*(h.normalized_pair_vs_Nt15-1),color=col,label=name)
 cross=u[u.Nt_cm3==nt].U.iloc[0]; axs[1,0].plot(cross,50,'o',color=col);axs[1,0].annotate(f'{cross:.2f} V',(cross,50),xytext=(0,10),textcoords='offset points',ha='center',fontsize=8,color=col)
axs[1,0].axhline(50,color='black',ls='--',lw=1)
axs[1,0].set(xlim=(10,25),ylim=(.02,1e4),xlabel='Device reverse voltage magnitude U (V)',ylabel='Self-consistent dark current (mA/cm²)',title='C  Fixed-current crossings shift by 8.07 V')
axs[1,0].legend(fontsize=8,loc='lower right')
axs[1,1].set(xlim=(10,20),ylim=(-.35,.08),xlabel='Device reverse voltage magnitude U (V)',ylabel='Deviation of (pair source / Nt) from Nt = 10¹⁵ (%)',title='D  Density-normalized source changes by only ~0.31%')
axs[1,1].axhline(0,color='black',lw=.5)
for ax in axs.flat:ax.grid(alpha=.16)
fig.suptitle('300 K: distinguish relative rise, amplitude, and a fixed-current threshold\nTop: finite dark graphs. Bottom: separate self-consistent PF-inspired device model.',fontsize=13)
fig.savefig(OUT/'matched_evidence_300K.png',dpi=180);fig.savefig(OUT/'matched_evidence_300K.pdf');plt.close(fig)
print(json.dumps(checks,indent=2))
