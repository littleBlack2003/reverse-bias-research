from pathlib import Path
import json,csv,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
ROOT=Path(__file__).resolve().parents[1];SOURCE=ROOT/'input_with_shunt_saved_states'
if not SOURCE.exists():SOURCE=ROOT.parent/'latest_calibrated_100K_20261004'
Ts=list(range(100,301,25));data={T:json.loads((ROOT/'data'/f'no_shunt_T{T}_N321.json').read_text()) for T in Ts}
assert all(d['success'] and d['audit_summary']['all_intrinsic_gates'] and d['audit_summary']['external_shunt_identically_zero'] for d in data.values())
colors={T:plt.colormaps['viridis'](.07+.84*i/8) for i,T in enumerate(Ts)}
plt.rcParams.update({'font.size':11,'axes.titlesize':13,'axes.labelsize':11,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white','axes.grid':True,'grid.alpha':.18,'savefig.facecolor':'white'})
curves={}
for T,d in data.items():
 m=d['metrics'];curves[T]=[r for r in d['curves']['light'] if 0<=r['V']<=m['Voc_V']+.040000001]
 assert curves[T][0]['V']==0 and curves[T][-1]['V']>=m['Voc_V']+.03999999
fig=plt.figure(figsize=(14,10),layout='constrained');gs=fig.add_gridspec(2,3,height_ratios=[1.3,1]);ax=fig.add_subplot(gs[0,:])
for T in Ts:
 r=curves[T];ax.plot([q['V'] for q in r],[1000*q['J_terminal_Acm2'] for q in r],c=colors[T],lw=2,label=f'{T} K')
ax.axhline(0,color='#555',lw=.8);ax.set(xlabel='Voltage (V)',ylabel='Signed current density (mA cm$^{-2}$)',title='All nine temperatures | true open-circuit crossing included',xlim=(0,1.18));ax.legend(ncols=5,loc='lower left',frameon=True,fontsize=10)
for i,T in enumerate([100,125,150]):
 ax=fig.add_subplot(gs[1,i]);r=curves[T];m=data[T]['metrics'];factor=1e6 if T==100 else 1e3
 ax.plot([q['V'] for q in r],[factor*q['J_terminal_Acm2'] for q in r],c=colors[T],lw=2.2)
 ax.axhline(0,color='#555',lw=.8);ax.scatter([m['Voc_V']],[0],marker='o',facecolor='white',edgecolor=colors[T],zorder=4)
 ax.scatter([m['Vmp_V']],[-m['Jmp_mAcm2']*(1000 if T==100 else 1)],marker='D',c=[colors[T]],s=26,zorder=4)
 unit='$\\mu$A cm$^{-2}$' if T==100 else 'mA cm$^{-2}$'
 ax.set(xlabel='Voltage (V)',ylabel=f'Signed current density ({unit})',title=f'{T} K | linear low-current zoom',xlim=(0,m['Voc_V']+.04))
 ax.text(.035,.97,f'$V_{{OC}}$ = {m["Voc_V"]:.6f} V\nFF = {100*m["FF"]:.3f}%\nDiamond: maximum power',transform=ax.transAxes,va='top',fontsize=10)
fig.suptitle('PM6:Y6 | external shunt removed | 100–300 K',fontsize=18,fontweight='bold')
fig.supxlabel('Full Fermi–Dirac + self-consistent Poisson/drift–diffusion; all axes linear. Frozen parameters; no refit.\nLow-temperature transport/recombination extensions are conditional, not experimentally validated.',fontsize=10)
for ext in ['png','pdf']:fig.savefig(ROOT/'figures'/f'PM6Y6_no_shunt_100-300K_LINEAR_JV.{ext}',dpi=210)
plt.close(fig)
fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for ax,key,label in zip(axes.flat,['Jsc_mAcm2','Voc_V','FF','Pmax_mWcm2'],['Short-circuit current (mA cm$^{-2}$)','Open-circuit voltage (V)','Fill factor (%)','Maximum power (mW cm$^{-2}$)']):
 vals=[data[T]['metrics'][key]*(100 if key=='FF' else 1) for T in Ts]
 ax.plot(Ts,vals,'-',c='#456578',lw=1.8);ax.scatter(Ts,vals,c=[colors[T] for T in Ts],s=48,zorder=3)
 ax.set(xlabel='Temperature (K)',ylabel=label,xticks=Ts,xlim=(94,306));ax.tick_params(axis='x',labelsize=9)
 if key=='Jsc_mAcm2':ax.annotate('100 K: 0.000387436',xy=(100,vals[0]),xytext=(132,5),fontsize=10,arrowprops={'arrowstyle':'-','color':'#777'})
 if key=='Pmax_mWcm2':ax.annotate('100 K: 9.64582 × 10$^{-5}$',xy=(100,vals[0]),xytext=(132,3),fontsize=10,arrowprops={'arrowstyle':'-','color':'#777'})
 if key=='Voc_V':ax.set_ylim(.79,1.16)
fig.suptitle('Photovoltaic metrics | no external shunt | all axes linear',fontsize=17,fontweight='bold')
fig.supxlabel('N = 321 primary curves. Removed only Jsh = V/Rsh; series resistance remains zero.\n300 K calibration parameters retained without refitting; low-temperature predictions remain conditional.',fontsize=10)
for ext in ['png','pdf']:fig.savefig(ROOT/'figures'/f'PM6Y6_no_shunt_100-300K_LINEAR_PV.{ext}',dpi=210)
plt.close(fig)
rows=[]
for T in Ts:
 for q in curves[T]:rows.append(dict(T_K=T,V_V=q['V'],J_mAcm2=1000*q['J_terminal_Acm2'],Pout_mWcm2=-1000*q['external_terminal_Wcm2'],J_intrinsic_Acm2=q['J_intrinsic_Acm2'],J_shunt_Acm2=q['J_shunt_Acm2'],current_closure_uncertainty_Acm2=q['current_closure_uncertainty_Acm2'],sign_resolved=q['terminal_sign_resolved'],gate_passed=q['gate_passed'],scaled_nonlinear_residual=q['scaled_nonlinear_residual'],state_origin=q['category']))
def write_csv(path,rows):
 with path.open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
write_csv(ROOT/'data'/'PM6Y6_no_shunt_100-300K_JV.csv',rows)
metrics=[]
for T in Ts:
 m=data[T]['metrics'].copy();m['FF_percent']=100*m.pop('FF');metrics.append(m)
write_csv(ROOT/'data'/'PM6Y6_no_shunt_100-300K_PV.csv',metrics)
comparison=[]
for T in Ts:
 old=json.loads((SOURCE/'data'/f'latest_T{T}_N321.json').read_text())['metrics'];new=data[T]['metrics'];assert old['Jsc_mAcm2']==new['Jsc_mAcm2']
 comparison.append(dict(T_K=T,Jsc_mAcm2=new['Jsc_mAcm2'],old_with_shunt_Voc_V=old['Voc_V'],new_no_shunt_Voc_V=new['Voc_V'],old_with_shunt_FF_percent=100*old['FF'],new_no_shunt_FF_percent=100*new['FF'],old_with_shunt_Pmax_mWcm2=old['Pmax_mWcm2'],new_no_shunt_Pmax_mWcm2=new['Pmax_mWcm2']))
write_csv(ROOT/'data'/'before_after_external_shunt.csv',comparison)
summary={'all_primary_success':True,'temperatures_K':Ts,'new_primary_solved_states':sum(d['counts']['new_solved_states'] for d in data.values()),'same_grid_reused_states_including_equilibrium':sum(d['counts']['reused_same_grid'] for d in data.values()),'plotted_CSV_rows':len(rows),'metrics':metrics,'minimum_primary_root_bracket_signal_to_closure':min(d['metrics']['minimum_bracket_signal_to_closure'] for d in data.values()),'primary_Jsc_matches_prior_exactly':True}
(ROOT/'validation'/'primary_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
