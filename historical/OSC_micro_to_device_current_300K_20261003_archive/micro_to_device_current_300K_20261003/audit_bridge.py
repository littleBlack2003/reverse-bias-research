#!/usr/bin/env python3
"""Read-only upstream audit, dimensional capacities, and exact all-rate scaling.
No field scan, refit, new kinetic solver, or claim of material parameters.
"""
from pathlib import Path
import json,hashlib,os
os.environ['MPLBACKEND']='Agg';os.environ['MPLCONFIGDIR']='/tmp/micro_device_mpl'
import numpy as np,pandas as pd,matplotlib.pyplot as plt
P=Path(__file__).resolve().parent;S=P/'inputs' if (P/'inputs').exists() else P.parent;O=P/'results';O.mkdir(exist_ok=True)
src=S/'full_network_bath_300K_20261003/results'
s=pd.read_csv(src/'summary.csv').set_index('case');b=s.loc['baseline_refined'];qu=s.loc['full_quantum_refined']
q=1.602176634e-19;d=100e-7;ell=12e-7;eps=3.5*8.8541878128e-14;F=1.5e6
files=[src/'summary.csv',src/'baseline_refined_edges.csv',src/'full_quantum_refined_edges.csv',src/'baseline_refined_states.csv',S/'full_network_bath_300K_20261003/inputs/spatial_model.py',S/'matched_evidence_300K_20261002/spatial_300K.csv',S/'trap_free_comparator_20261002/room_temperature_300K/results/rise_metrics_300K.csv',S/'density_knee_disentangling_20261002/REPORT.md']
checks={};fluxrows=[];caps=[]
for tag,row in [('baseline',b),('full_quantum',qu)]:
 e=pd.read_csv(src/(tag+'_refined_edges.csv'));a=e.groupby('process')[['net_flux_s','gross_flux_s']].sum();a['forward_event_s']=(a.gross_flux_s+a.net_flux_s)/2;a['reverse_event_s']=(a.gross_flux_s-a.net_flux_s)/2;a['case']=tag;fluxrows.extend(a.reset_index().to_dict('records'))
 sign=np.where(a.index=='bathR',-1,1);err=max(abs(a.net_flux_s*sign/row.R_s-1));disp=e.electron_displacement_nm_s.sum();checks[tag+'_flux_relative_error']=err;checks[tag+'_displacement_relative_error']=abs(disp/(12*row.R_s)-1)
 for process,part in e.groupby('process'):
  kmax=part['reverse_rate_s' if process=='bathR' else 'forward_rate_s'].max()
  caps.append(dict(case=tag,process=process,forward_cut_cap_s=kmax,net_R_s=row.R_s))
 assert err<2e-14 and abs(disp/(12*row.R_s)-1)<2e-14
pd.DataFrame(fluxrows).to_csv(O/'flux_reconstruction.csv',index=False);pd.DataFrame(caps).to_csv(O/'fixed_rate_cut_bounds.csv',index=False)
rows=[]
for jm in [.1,1,10,25,50,100]:
 J=jm*1e-3;N=J/(q*d*b.R_s);Nq=J/(q*d*qu.R_s);Nlocal=J/(q*ell*b.R_s)
 rows.append(dict(target_J_mAcm2=jm,required_independent_N_cm3=N,quantum_required_N_cm3=Nq,required_5site_orbitals_cm3=5*N,mean_center_spacing_nm=N**(-1/3)*1e7,N_for_local_displacement_only_cm3=Nlocal,uncompensated_DeltaF_MVcm=q*N*b.mean_total_charge_e*d/eps/1e6,uncompensated_DeltaF_over_F=q*N*b.mean_total_charge_e*d/eps/F))
a=pd.DataFrame(rows);a.to_csv(O/'conditional_capacity.csv',index=False)
rows=[]
for N in [1e15,1e16,1e17,1e18,2e18]:
 scale=.05/(q*d*N*b.R_s);alpha=np.sqrt(scale)
 rows.append(dict(assumed_N_cm3=N,target_J_mAcm2=50,all_H_multiplier=alpha,all_rate_multiplier=scale,inner_H_microeV=10*alpha,total_graph_H2_eV2=b.total_coupling_H2_eV2*scale,scaled_max_exit_s=b.maximum_exit_rate_s*scale,relaxation_time_for_1pct_event_ns=-np.log(.99)/(b.maximum_exit_rate_s*scale)*1e9))
c=pd.DataFrame(rows);c.to_csv(O/'uniform_coupling_scaling.csv',index=False)
Ncap=1e19/5;N50=float(a.query('target_J_mAcm2==50').required_independent_N_cm3.iloc[0]);Nfield=.1*eps*F/(q*abs(b.mean_total_charge_e)*d)
checks.update(R_baseline_s=b.R_s,R_quantum_s=qu.R_s,relative_quantum_change=qu.R_s/b.R_s-1,J_at_N1e15_mAcm2=q*d*1e15*b.R_s*1e3,graph_current_A=q*b.R_s,required_N50_cm3=N50,example_orbital_pool_cm3=1e19,example_disjoint_graph_Ncap_cm3=Ncap,example_Jcap_mAcm2=q*d*Ncap*b.R_s*1e3,uncompensated_Ncap_for_10pct_field_cm3=Nfield,uncompensated_Jcap_for_10pct_field_mAcm2=q*d*Nfield*b.R_s*1e3,minimum_compensation_fraction_for_50mA_10pct_field=1-Nfield/N50,illustrative_12nm_cube_Ncap_cm3=ell**-3,illustrative_cube_Jcap_mAcm2=q*d*ell**-3*b.R_s*1e3,local_to_full_displacement_ratio=ell/d,global_uniform_scaling_identity='Q(alpha)=alpha^2 Q(1), p(alpha)=p(1), R(alpha)=alpha^2 R(1)',all_assertions_passed=True)
# Existing 300 K field points only, not a new scan.
z=pd.read_csv(S/'matched_evidence_300K_20261002/spatial_300K.csv');z=z.query("variant=='explicit'").copy();z['N_for_50mA_cm3']=.05/(q*d*z.R_s);z.to_csv(O/'existing_explicit_field_points.csv',index=False)
checks['R3_over_R1p5_existing']=float(z.loc[np.isclose(z.F_MVcm,3),'R_s'].iloc[0]/b.R_s)
(O/'checks.json').write_text(json.dumps(checks,indent=2))
(P/'SOURCE_PROVENANCE.json').write_text(json.dumps({str(f.relative_to(S)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files},indent=2))
fig,ax=plt.subplots(2,2,figsize=(11,8),layout='constrained');N=np.logspace(14,20,300)
ax[0,0].loglog(N,q*d*N*b.R_s*1e3,label='Full-device collection closure',color='#156082');ax[0,0].loglog(N,q*ell*N*b.R_s*1e3,'--',label='12 nm represented displacement only',color='#e97132');ax[0,0].axhline(50,color='0.4',ls=':');ax[0,0].scatter([N50],[50],color='#156082');ax[0,0].set(xlabel='Independent graph density (cm$^{-3}$)',ylabel='Conditional current (mA cm$^{-2}$)',title='A  The boundary closure changes the mapping');ax[0,0].legend(fontsize=8)
ax[0,1].loglog(a.required_independent_N_cm3,a.uncompensated_DeltaF_over_F,color='#a02b93',marker='o');ax[0,1].axhline(.1,color='0.4',ls=':');ax[0,1].set(xlabel='Independent graph density (cm$^{-3}$)',ylabel=r'Uncompensated $|\Delta F|/F$',title='B  Charge cannot be multiplied away');ax[0,1].text(.05,.85,'Frozen graph charge +0.10143 e\nNo reservoir/background compensation',transform=ax[0,1].transAxes,fontsize=9)
ax[1,0].loglog(c.assumed_N_cm3,c.all_H_multiplier,marker='o',color='#156082');ax[1,0].set(xlabel='Assumed independent graph density (cm$^{-3}$)',ylabel='Required uniform coupling multiplier',title='C  Exact fixed-generator rescaling to 50 mA');ax[1,0].text(.05,.82,'All physical bonds and contacts scaled\nNot a validated material change',transform=ax[1,0].transAxes,fontsize=9)
fl=pd.DataFrame(fluxrows).query("case=='baseline'").set_index('process').loc[['Vo_V','V_D','D_C','C_Co']];x=np.arange(4)
ax[1,1].bar(x-.18,fl.forward_event_s/1e3,.36,label='Forward event flux',color='#156082');ax[1,1].bar(x+.18,fl.reverse_event_s/1e3,.36,label='Reverse event flux',color='#e97132');ax[1,1].axhline(b.R_s/1e3,color='black',ls='--',label='Net R (already subtracted)');ax[1,1].set_xticks(x,fl.index);ax[1,1].set(ylabel='Events / 1000 s',title='D  Do not multiply a second recrossing factor');ax[1,1].legend(fontsize=8)
fig.suptitle('300 K | 1.5 MV/cm | fixed 12 nm graph | 50 mA/cm$^2$ is diagnostic only',fontsize=13)
fig.savefig(P/'micro_to_device_bridge.png',dpi=180);plt.close(fig)
print(json.dumps(checks,indent=2));print(a.to_string(index=False));print(c.to_string(index=False))
