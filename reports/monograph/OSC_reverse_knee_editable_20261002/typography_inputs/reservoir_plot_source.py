#!/usr/bin/env python3
import os
os.environ.setdefault('MPLBACKEND','Agg');os.environ.setdefault('MPLCONFIGDIR','/tmp/quantum_reservoir_mpl')
from pathlib import Path
import json,csv,sys,time,hashlib,platform
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import expit
import matplotlib;import matplotlib.pyplot as plt
import mpmath as mm
import reservoir as r
ROOT=Path(__file__).resolve().parent;R=ROOT/'results';FIG=ROOT/'figures'
def write(name,rows):
 with (R/name).open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def save(name):plt.tight_layout();plt.savefig(FIG/name,dpi=170);plt.close()
def independent_quad(T,da,mur,alpha,kind,W):
 lo,hi=(-W,0) if alpha=='V' else (0,W)
 # Independent adaptive scalar quadrature, rescaled to avoid rare-rate underflow.
 outputs=[]
 for add in [True,False]:
  def logf(u):
   lf=-np.logaddexp(0,(u-mur)/(r.KB*T)) if add else -np.logaddexp(0,-(u-mur)/(r.KB*T))
   return float(r.core.logkernel(da-u if add else u-da,T,kind,H=r.HRES))+lf-np.log(W)
  mesh=np.linspace(lo,hi,1001);anchor=max(logf(u) for u in mesh)
  fun=lambda u:np.exp(logf(u)-anchor)
  mass,err=quad(fun,lo,hi,epsabs=1e-13,epsrel=2e-12,limit=500)
  moment,err2=quad(lambda u:u*fun(u),lo,hi,epsabs=1e-13,epsrel=2e-12,limit=500)
  outputs.append((np.log(mass)+anchor,moment/mass,err/mass))
 return outputs

def elastic_limit_tests():
 rows=[];T=300;W=1.;mu=-.2;lo=-1.;hi=0.
 for a in [-.25,0]:
  target=expit((mu-a)/(r.KB*T))*(.5 if a==0 else 1)
  for lam in [.01,.001,.0001,.00001,.000001,1e-8]:
   sd=np.sqrt(2*lam*r.KB*T);mean=a+lam
   zlo=max((lo-mean)/sd,-14);zhi=min((hi-mean)/sd,14)
   integ=quad(lambda z:np.exp(-z*z/2)/np.sqrt(2*np.pi)*expit((mu-mean-sd*z)/(r.KB*T)),zlo,zhi,epsabs=1e-14,epsrel=1e-12)[0]
   rows.append(dict(T_K=T,a_from_band_edge_eV=a,mu_from_band_edge_eV=mu,lambda_eV=lam,rate_over_Gamma0=integ,expected_elastic_limit=target,relative_error=abs(integ/target-1),note='at hard edge the smooth Gaussian limit has half spectral mass'))
 write('elastic_limit_matching.csv',rows)
 return rows[-1]['relative_error']

def main():
 start=time.time();tests={};Ts=r.core.P['temperatures_K'];Fs=np.round(np.arange(0,3.00001,.005),8)
 params={'scope':'new finite normalized-DOS boundary control; not material model or parameter fit','internal_model':'inputs/model.py and inputs/parameters.json, unchanged copied checkpoint','temperature_grid_K':Ts,'field_range_MVcm':[0,3],'field_step_MVcm':.005,'W_baseline_eV':1,'DOS_V':'uniform in [EV-W,EV], normalized to 1','DOS_C':'uniform in [EC,EC+W], normalized to 1','H_reservoir_eV':r.HRES,'each_reservoir_H_squared_eV2':r.HRES**2,'two_reservoirs_H_squared_budget_eV2':2*r.HRES**2,'internal_H_squared_budget_eV2':2e-10,'matching':'Gamma0=2*pi*Hres^2/(hbar*W0)=1e6/s at W0=1eV; Hres fixed for width controls','main_drive':'delta_mu_eV = numerical F_MVcm over 10 nm span; independent controls break this tie','quadrature':'24 Gauss-Legendre nodes per interval split at Gaussian centers and Fermi edge/shoulders','precision_dps':130,'bath_nuclear_kernel':'same classical/thermal vibronic kernel as internal transfer, exploratory assumption','conditional_current_conversion':'q*1e15/cm^3*100nm*R, assuming independent uniformly populated graphs and full collection'}
 (ROOT/'parameters.json').write_text(json.dumps(params,indent=2))
 # Quadrature and exact reverse spectral relation before the graph scans.
 conv=[];ind=[]
 for T in [80,150,300,330]:
  for alpha in ['V','C']:
   for da,mur in [(0,.65 if alpha=='V' else -.65),(.041141844222857,.1),(-.041141844222857,-.1),(.123425532668571,0.)]:
    for kind in ['classical','quantum']:
     vals=[r.exchange(T,da,mur,alpha,kind,1.,n) for n in [12,24,48]]
     conv.append(dict(T_K=T,alpha=alpha,da_eV=da,mu_relative_eV=mur,kernel=kind,order12_vs24_max_logerror=max(abs(vals[0][z]-vals[1][z]) for z in ['log_add','log_out']),order24_vs48_max_logerror=max(abs(vals[1][z]-vals[2][z]) for z in ['log_add','log_out']),DB_logerror=vals[1]['direct_DB_log_error'],mean_energy_error_eV=vals[1]['conditional_mean_energy_error_eV'],DOS_mass_error=abs(vals[1]['DOS_mass']-1),nodes_order24=vals[1]['node_count']))
 for T,da,mur,alpha,kind in [(80,.041141844222857,.65,'V','quantum'),(80,-.041141844222857,-.65,'C','quantum'),(80,.123425532668571,0.,'V','quantum'),(80,0,.1,'C','classical'),(300,.041141844222857,.65,'V','quantum'),(330,-.041141844222857,-.1,'C','classical')]:
  gl=r.exchange(T,da,mur,alpha,kind);aq=independent_quad(T,da,mur,alpha,kind,1.)
  ind.append(dict(T_K=T,da_eV=da,mu_relative_eV=mur,alpha=alpha,kernel=kind,adaptive_vs_segmented_add_logerror=aq[0][0]-gl['log_add'],adaptive_vs_segmented_out_logerror=aq[1][0]-gl['log_out'],adaptive_vs_segmented_mean_add_eV=aq[0][1]-gl['mean_relative_energy_add_eV'],adaptive_vs_segmented_mean_out_eV=aq[1][1]-gl['mean_relative_energy_out_eV'],adaptive_relative_error_estimate=max(a[2] for a in aq)))
 write('quadrature_convergence.csv',conv);write('independent_adaptive_quadrature.csv',ind)
 tests['quadrature24_vs48_max_logerror']=max(x['order24_vs48_max_logerror'] for x in conv);tests['spectral_DB_max_logerror']=max(abs(x['DB_logerror']) for x in conv);tests['conditional_energy_mean_max_error_eV']=max(abs(x['mean_energy_error_eV']) for x in conv);tests['normalized_DOS_mass_max_error']=max(x['DOS_mass_error'] for x in conv);tests['independent_quad_max_logerror']=max(abs(x[z]) for x in ind for z in ['adaptive_vs_segmented_add_logerror','adaptive_vs_segmented_out_logerror'])
 tests['elastic_hard_edge_limit_final_relative_error']=elastic_limit_tests()
 # Main same-field/same-drive comparison with the old boundary, no retuning.
 rows=[]
 for T in Ts:
  for kind in ['classical','quantum']:
   for bath in ['legacy','finite_dos']:
    for F in Fs:rows.append(r.solve(float(F),T,kind,bath_model=bath))
    print('main',T,kind,bath,'elapsed',round(time.time()-start,1),flush=True)
 write('matched_boundary_full_scan.csv',rows)
 tests['main_scan_rows']=len(rows);tests['minimum_probability']=min(x['min_probability'] for x in rows);tests['continuity_max_relative_error_nonzero_drive']=max(x['continuity_relative_error'] for x in rows if x['delta_mu_eV']!=0)
 tests['heat_vs_chemical_max_relative_error_nonzero_drive']=max(x['heat_vs_chemical_relative_error'] for x in rows if x['delta_mu_eV']!=0);tests['heat_partition_max_relative_error_nonzero_drive']=max(x['heat_partition_relative_error'] for x in rows if x['delta_mu_eV']!=0);tests['entropy_vs_heat_max_relative_error_nonzero_drive']=max(x['entropy_vs_heat_relative_error'] for x in rows if x['delta_mu_eV']!=0)
 # Independent static field versus electrochemical driving.
 controls=[];eq=[];edgerows=[];states=[];cycles=[]
 basis=json.loads((ROOT/'inputs/cycle_basis.json').read_text())
 edge_index={(i,j):(ee,1) for ee,(i,j,n) in enumerate(r.core.EDGES)}
 edge_index.update({(j,i):(ee,-1) for ee,(i,j,n) in enumerate(r.core.EDGES)})
 for T in [80,300]:
  for F in [0,.5,1,1.5,2,3]:
   for dmu in [-1,-.5,0,.5,1,1.5,2,3]:
    for kind in ['classical','quantum']:
     row,p,rat,Q,ed,E,mu=r.solve(F,T,kind,delta_mu=dmu,detail=True);controls.append(row)
     for cy in basis:
      aff=mm.mpf(0);nv=0
      for e in cy['oriented_edges']:
       i=r.core.INDEX[tuple(map(int,e['from']))];j=r.core.INDEX[tuple(map(int,e['to']))];ee,sgn=edge_index[i,j]
       aff+=sgn*mm.log(rat[ee][0]/rat[ee][1])
       if e['process']=='V':nv+=sgn
      pred=mm.mpf(str(dmu))*nv/(mm.mpf(str(r.KB))*T)
      cycles.append(dict(T_K=T,F_MVcm=F,delta_mu_eV=dmu,kernel=kind,cycle=cy['cycle'],measured_affinity=float(aff),expected_affinity=float(pred),absolute_error=float(abs(aff-pred))))
     if dmu==0:
      weights=[mm.exp(-ee/(mm.mpf(str(r.KB))*T)) for ee in E];Z=sum(weights);pi=[z/Z for z in weights]
      ep=max(abs(p[i]/pi[i]-1) for i in range(8));edgeerr=max(abs(p[i]*kf-p[j]*kr)/max(p[i]*kf,p[j]*kr) for (i,j,n),(kf,kr,*_) in zip(r.core.EDGES,rat))
      eq.append(dict(T_K=T,F_MVcm=F,kernel=kind,common_mu_eV=0.,stationary_relative_error=float(ep),edge_relative_error=float(edgeerr),net_R_s=row['R_s']))
  for F in [.5,1.5,2,3]:
   for kind in ['classical','quantum']:
    for bath in ['legacy','finite_dos']:
     row,p,rat,Q,ed,E,mu=r.solve(F,T,kind,bath_model=bath,detail=True)
     for e in ed:edgerows.append(dict(T_K=T,F_MVcm=F,kernel=kind,bath=bath,**e))
     for s,pp in zip(r.core.STATES,p):states.append(dict(T_K=T,F_MVcm=F,kernel=kind,bath=bath,state=''.join(map(str,s)),probability=float(pp)))
 for T in Ts:
  if T in [80,300]:continue
  for F in [1.5,2,3]:
   row,p,*_=r.solve(F,T,'quantum',detail=True)
   for state,pp in zip(r.core.STATES,p):states.append(dict(T_K=T,F_MVcm=F,kernel='quantum',bath='finite_dos',state=''.join(map(str,state)),probability=float(pp)))
 write('fundamental_cycle_affinities.csv',cycles);tests['fundamental_cycle_max_affinity_error']=max(x['absolute_error'] for x in cycles)
 write('independent_field_and_mu_controls.csv',controls);write('common_mu_nonzero_field.csv',eq);write('representative_edge_rates_fluxes.csv',edgerows);write('representative_state_probabilities.csv',states)
 tests['equilibrium_max_relative_state_error']=max(x['stationary_relative_error'] for x in eq);tests['equilibrium_max_relative_edge_error']=max(x['edge_relative_error'] for x in eq);tests['sign_follows_delta_mu']=all(x['R_s']*x['delta_mu_eV']>0 for x in controls if x['delta_mu_eV']!=0)
 # Explicit nonzero common mu controls, as in the preceding phase.
 eq23=[]
 for T in [80,300]:
  for F in [1,2,3]:
   row,p,rat,Q,ed,E,mu=r.solve(F,T,'quantum',delta_mu=0.,mean_mu=.023,detail=True)
   lp=[-(E[i]-mm.mpf('.023')*sum(s))/(mm.mpf(str(r.KB))*T) for i,s in enumerate(r.core.STATES)];Z=sum(mm.exp(z) for z in lp)
   eq23.append(dict(T_K=T,F_MVcm=F,common_mu_eV=.023,max_relative_probability_error=float(max(abs(p[i]/(mm.exp(lp[i])/Z)-1) for i in range(8))),net_R_s=row['R_s']))
 write('common_mu_0023_control.csv',eq23)
 # Width changes redistribute a FIXED reservoir coupling budget.
 width=[]
 for T in [80,300]:
  for F in [1,1.5,2,3]:
   for W in [.25,.5,1.,2.]:
    for kind in ['classical','quantum']:width.append(r.solve(F,T,kind,W=W))
 write('fixed_coupling_bandwidth_controls.csv',width)
 # Same W/H/spectral weight: centered broad contact allows states on gap side.
 supportrows=[]
 for T in [80,300]:
  for F in [.5,1,1.5,2,3]:
   for kind in ['classical','quantum']:
    for support in ['one_sided','centered']:
     supportrows.append(r.solve(F,T,kind,support=support))
 write('one_sided_vs_broad_contact_support.csv',supportrows)
 supportconv=[]
 for T in [80,300]:
  for alpha in ['V','C']:
   for da,mur in [(0,.65),(.123425532668571,0),(-.041141844222857,-.65)]:
    a=r.exchange(T,da,mur,alpha,order=24,support='centered');b=r.exchange(T,da,mur,alpha,order=48,support='centered')
    supportconv.append(dict(T_K=T,alpha=alpha,da_eV=da,mu_relative_eV=mur,max_lograte_error=max(abs(a[z]-b[z]) for z in ['log_add','log_out']),DB_error=a['direct_DB_log_error'],DOS_mass_error=abs(a['DOS_mass']-1)))
 write('centered_support_convergence.csv',supportconv)
 tests['centered_support_quadrature_max_logerror']=max(x['max_lograte_error'] for x in supportconv)
 # Separate internal versus boundary quantum corrections (four kernels, fixed H).
 factorial=[]
 for T in [80,300]:
  for F in [1,1.5,2,3]:
   for ki in ['classical','quantum']:
    for kb in ['classical','quantum']:factorial.append(r.solve(F,T,ki,bath_kind=kb))
 write('internal_vs_bath_kernel_factorial.csv',factorial)
 # Wrong reverse Fermi occupancy produces spurious common-mu currents.
 bad=[]
 for T in [80,300]:
  for F in [0,.5,1,1.5,2,3]:
   for wrong in [False,True]:
    row=r.solve(F,T,'quantum',delta_mu=0.,mean_mu=.023,wrong_fermi=wrong)
    bad.append(dict(T_K=T,F_MVcm=F,wrong_reverse_fermi=wrong,common_mu_eV=.023,R_s=row['R_s'],spectral_DB_error=row['direct_integral_DB_max_log_error'],rate_entropy_kB_s=row['entropy_production_kB_s'],physical_heat_over_kBT_s=row['total_bath_heat_eV_s']/(r.KB*T)))
 write('wrong_reverse_fermi_counterexample.csv',bad)
 # Independent tree and high-precision check of the new integrated-rate graph.
 trees=r.core.spanning_trees();cross=[]
 for T in [80,300]:
  for F,dmu in [(0,0),(.5,.5),(1.5,1.5),(2,2),(3,3),(3,-1),(1.5,.01)]:
   a,p,rat,Q,ed,E,mu=r.solve(F,T,delta_mu=dmu,detail=True)
   tp=r.core.tree_stationary(rat,trees);pe=np.array([float(x) for x in p]);treeerr=max(abs(tp/pe-1));j1=sum(p[i]*rr[0]-p[j]*rr[1] for (i,j,n),rr in zip(r.core.EDGES,rat) if n=='DC')
   b,p2,rat2,*_=r.solve(F,T,delta_mu=dmu,dps=170,detail=True);j2=sum(p2[i]*rr[0]-p2[j]*rr[1] for (i,j,n),rr in zip(r.core.EDGES,rat2) if n=='DC')
   cross.append(dict(T_K=T,F_MVcm=F,delta_mu_eV=dmu,tree_probability_relative_error=treeerr,precision130_vs170_relative_R_error=float(abs(j2/j1-1)) if dmu else 0))
 write('independent_graph_solution_checks.csv',cross)
 tests['tree_max_probability_relative_error']=max(x['tree_probability_relative_error'] for x in cross);tests['precision130_vs170_max_relative_current_error']=max(x['precision130_vs170_relative_R_error'] for x in cross)
 # Peak and main representative tables.
 peaks=[];reps=[]
 for T in Ts:
  for kind in ['classical','quantum']:
   for bath in ['legacy','finite_dos']:
    cr=[x for x in rows if x['T_K']==T and x['internal_kernel']==kind and x['bath_model']==bath];pk=max(cr,key=lambda z:z['R_s']);peaks.append({key:pk[key] for key in ['T_K','internal_kernel','bath_model','F_MVcm','R_s','conditional_J_Acm2']})
    reps.extend(x for x in cr if x['F_MVcm'] in [0,.5,1,1.5,2,3])
 write('peak_capacities.csv',peaks);write('representative_main_points.csv',reps)
 # Figures: all curves explicitly exploratory and local-field based.
 fig,axes=plt.subplots(1,2,figsize=(11,4.4),sharey=True)
 for kind,ax in zip(['classical','quantum'],axes):
  for T,c in [(80,'tab:blue'),(300,'tab:red')]:
   for bath,ls in [('legacy','--'),('finite_dos','-')]:
    cr=[x for x in rows if x['T_K']==T and x['internal_kernel']==kind and x['bath_model']==bath and x['F_MVcm']>0]
    ax.semilogy([x['F_MVcm'] for x in cr],[x['R_s'] for x in cr],ls,color=c,label=f'{T} K, '+('constant Gamma' if bath=='legacy' else 'finite band'))
  ax.set(xlabel='Local favorable field (MV/cm)',ylabel='Net cycles per graph (s$^{-1}$)',title=kind,ylim=(1e-40,2e6));ax.grid(alpha=.2);ax.legend(fontsize=8)
 fig.suptitle('Fixed coupling budgets: local-band spectral exchange changes the bottleneck');save('finite_bands_vs_constant_gamma.png')
 fig,ax=plt.subplots(1,2,figsize=(11,4.4))
 for F,c in [(1.5,'tab:blue'),(2,'tab:orange'),(3,'tab:green')]:
  for bath,ls in [('legacy','--'),('finite_dos','-')]:
   cr=[x for x in rows if x['internal_kernel']=='quantum' and x['bath_model']==bath and x['F_MVcm']==F]
   ax[0].semilogy(Ts,[x['R_s'] for x in cr],ls,marker='.',color=c,label=f'{F:g} MV/cm, '+('Gamma' if bath=='legacy' else 'band'))
  cr=sorted([x for x in states if x['kernel']=='quantum' and x['bath']=='finite_dos' and x['F_MVcm']==F and x['state']=='001'],key=lambda z:z['T_K']);ax[1].plot([x['T_K'] for x in cr],[x['probability'] for x in cr],marker='o',label=f'{F:g} MV/cm')
 ax[0].set(xlabel='Temperature (K)',ylabel='Quantum net cycles (s$^{-1}$)',title='Local band resets restore activation');ax[1].set(xlabel='Temperature (K)',ylabel='Probability of bound-pair state 001',title='Charge-conditioned pair export bottleneck')
 for a in ax:a.grid(alpha=.2);a.legend(fontsize=8)
 save('temperature_and_pair_accumulation.png')
 fig,axes=plt.subplots(1,2,figsize=(11,4.4))
 for T,ax in zip([80,300],axes):
  for F in [0,1,2,3]:
   cr=[x for x in controls if x['T_K']==T and x['internal_kernel']=='quantum' and x['F_MVcm']==F];ax.plot([x['delta_mu_eV'] for x in cr],[x['R_s'] for x in cr],marker='.',label=f'F={F:g} MV/cm')
  ax.set_yscale('symlog',linthresh=1e-5);ax.set(xlabel='Independent reservoir muV - muC (eV)',ylabel='Signed quantum net cycles (s$^{-1}$)',title=f'{T} K; independent field and reservoir drive');ax.grid(alpha=.2);ax.legend(fontsize=8)
 save('independent_field_and_drive.png')
 fig,ax=plt.subplots(1,2,figsize=(11,4.4))
 for T,a in zip([80,300],ax):
  cr=[x for x in rows if x['T_K']==T and x['internal_kernel']=='quantum' and x['bath_model']=='finite_dos' and x['F_MVcm']>.2]
  for key,label in [('phonon_bath_heat_eV_s','phonon bath'),('electron_bath_heat_eV_s','electron reservoirs'),('total_bath_heat_eV_s','total')]:a.plot([x['F_MVcm'] for x in cr],[x[key]/x['R_s'] for x in cr],label=label)
  a.set(xlabel='Local favorable field (MV/cm)',ylabel='Heat per net cycle (eV)',title=f'{T} K; signed bath-resolved heat');a.axhline(0,color='k',lw=.6);a.grid(alpha=.2);a.legend(fontsize=8)
 save('physical_heat_partition.png')
 fig,ax=plt.subplots(figsize=(7,4.4))
 for T in [80,300]:
  cr=[x for x in bad if x['T_K']==T and x['wrong_reverse_fermi']];ax.plot([x['F_MVcm'] for x in cr],[x['R_s'] for x in cr],marker='o',label=f'{T} K')
 ax.set_yscale('symlog',linthresh=1e-20);ax.set(xlabel='Static field with common mu (MV/cm)',ylabel='False equilibrium cycle current (s$^{-1}$)',title='Wrong reverse Fermi factor breaks thermal balance');ax.legend();ax.grid(alpha=.2);save('wrong_fermi_counterexample.png')
 fig,axes=plt.subplots(1,2,figsize=(11,4.4))
 for T,ax in zip([80,300],axes):
  for support,ls in [('one_sided','-'),('centered','--')]:
   cr=[x for x in supportrows if x['T_K']==T and x['internal_kernel']=='quantum' and x['DOS_support']==support]
   ax.semilogy([x['F_MVcm'] for x in cr],[x['R_s'] for x in cr],ls,marker='o',label=support+' DOS')
  cr=[x for x in reps if x['T_K']==T and x['internal_kernel']=='quantum' and x['bath_model']=='legacy' and x['F_MVcm']>0]
  ax.semilogy([x['F_MVcm'] for x in cr],[x['R_s'] for x in cr],':',marker='s',label='constant Gamma')
  ax.set(xlabel='Local favorable field (MV/cm)',ylabel='Quantum net cycles (s$^{-1}$)',title=f'{T} K; same reservoir H squared and width');ax.grid(alpha=.2);ax.legend(fontsize=8)
 save('same_budget_spectral_support_control.png')
 assert tests['quadrature24_vs48_max_logerror']<1e-9 and tests['independent_quad_max_logerror']<1e-9
 assert tests['spectral_DB_max_logerror']<1e-10 and tests['conditional_energy_mean_max_error_eV']<1e-12
 assert tests['normalized_DOS_mass_max_error']<1e-13 and tests['minimum_probability']>0
 assert tests['equilibrium_max_relative_state_error']<1e-65 and tests['equilibrium_max_relative_edge_error']<1e-65
 assert tests['continuity_max_relative_error_nonzero_drive']<1e-65 and tests['heat_vs_chemical_max_relative_error_nonzero_drive']<1e-65
 assert tests['heat_partition_max_relative_error_nonzero_drive']<1e-65 and tests['entropy_vs_heat_max_relative_error_nonzero_drive']<1e-65
 assert tests['sign_follows_delta_mu'] and tests['tree_max_probability_relative_error']<1e-10
 assert tests['fundamental_cycle_max_affinity_error']<1e-70 and tests['centered_support_quadrature_max_logerror']<1e-9
 tests['all_assertions_passed']=True;tests['elapsed_s']=time.time()-start
 (R/'validation_summary.json').write_text(json.dumps(tests,indent=2));print(json.dumps(tests,indent=2),flush=True)
 prov={'date':'2026-10-02','python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'mpmath':mm.__version__,'matplotlib':matplotlib.__version__,'platform':platform.platform(),'command':'python -B run_reservoir_audit.py','sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'reservoir.py',ROOT/'run_reservoir_audit.py',ROOT/'parameters.json',ROOT/'THEORY_BEFORE_CODE.md',ROOT/'inputs/model.py',ROOT/'inputs/parameters.json']}}
 (R/'provenance.json').write_text(json.dumps(prov,indent=2))
if __name__=='__main__':main()
