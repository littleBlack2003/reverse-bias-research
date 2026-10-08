#!/usr/bin/env python3
import os
os.environ.setdefault('MPLBACKEND','Agg');os.environ.setdefault('MPLCONFIGDIR','/tmp/quantum_cycle_mpl')
from pathlib import Path
import json,csv,hashlib,time,sys,platform,itertools
import numpy as np
import scipy
from scipy.special import logsumexp
from scipy.linalg import expm
import matplotlib;import matplotlib.pyplot as plt
import mpmath as mm
import model as m
ROOT=Path(__file__).resolve().parent; R=ROOT/'results';FIG=ROOT/'figures'
R.mkdir(exist_ok=True);FIG.mkdir(exist_ok=True)
def write(name,rows):
 with (R/name).open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def savefig(name):
 plt.tight_layout();plt.savefig(FIG/name,dpi=170);plt.close()
def cycle_basis():
 # Deterministic spanning-tree fundamental cycles, signed orientations.
 adj=[[] for _ in range(8)];visited={0};queue=[0];tree=[]
 for i in queue:
  for e,(a,b,_) in enumerate(m.EDGES):
   j=b if a==i else a if b==i else None
   if j is not None and j not in visited:
    visited.add(j);queue.append(j);tree.append(e);adj[a].append((b,e));adj[b].append((a,e))
 result=[]
 for e,(a,b,_) in enumerate(m.EDGES):
  if e in tree:continue
  queue=[(b,[])];seen={b};path=None
  for v,pp in queue:
   if v==a:path=pp;break
   for w,ee in adj[v]:
    if w not in seen:seen.add(w);queue.append((w,pp+[(ee,1 if m.EDGES[ee][0]==v else -1)]))
  result.append([(e,1)]+path)
 return result

def finite_closed_reservoir():
 # Five finite orbitals RV,V,D,C,RC, two electrons, no external particle bath.
 # Coulomb includes only V,D,C; the two finite bath levels have stated energies.
 # High-precision exponential avoids probability drift at long time.
 F=2.;T=300;mm.mp.dps=100;beta=1/(mm.mpf(str(m.KB))*T)
 sts=[s for s in itertools.product([0,1],repeat=5) if sum(s)==2];ix={s:i for i,s in enumerate(sts)}
 inner=m.energy(F,True);es=[inner[m.INDEX[s[1:4]]]+mm.mpf(F)/2*(s[0]-s[4]) for s in sts];n=len(sts);Q=mm.zeros(n);edges=[]
 for i,s in enumerate(sts):
  for a,b in [(0,1),(1,2),(2,3),(3,4)]:
   if s[a] and not s[b]:
    t=list(s);t[a]=0;t[b]=1;j=ix[tuple(t)];dg=es[j]-es[i]
    if a in [1,2]:kf=mm.exp(mm.mpf(str(float(m.logkernel(float(dg),T)))));kr=kf*mm.exp(beta*dg)
    else:kf=mm.mpf('1e6')/(1+mm.exp(beta*dg));kr=mm.mpf('1e6')/(1+mm.exp(-beta*dg))
    Q[j,i]+=kf;Q[i,j]+=kr;Q[i,i]-=kf;Q[j,j]-=kr;edges.append((i,j,a,kf,kr))
 pp=[mm.exp(-e*beta) for e in es];ZZ=sum(pp);pe=[p/ZZ for p in pp];pi=mm.zeros(n,1);pi[ix[(1,1,0,0,0)]]=1
 rows=[]
 for t in np.r_[0,np.logspace(-9,-1,81)]:
  p=mm.expm(Q*mm.mpf(str(t)))*pi
  current=sum(p[i]*kf-p[j]*kr for i,j,a,kf,kr in edges if a==3)
  assert min(p)>=0 and max(p)<=1 and abs(sum(p)-1)<mm.mpf('1e-85')
  rows.append(dict(T_K=T,F_MVcm=F,time_s=t,finite_drain_occupation=float(sum(p[i]*s[4] for i,s in enumerate(sts))),finite_source_occupation=float(sum(p[i]*s[0] for i,s in enumerate(sts))),instantaneous_drain_flux_s=float(current),min_probability=float(min(p)),norm_error=float(abs(sum(p)-1))))
 write('closed_finite_bath_relaxation.csv',rows)
 pairerr=max(abs(pe[i]*kf-pe[j]*kr)/max(pe[i]*kf,pe[j]*kr) for i,j,a,kf,kr in edges)
 return {'state_count':n,'total_electrons':2,'canonical_edge_relative_error':float(pairerr),'final_drain_occupation':rows[-1]['finite_drain_occupation'],'final_drain_flux_s':rows[-1]['instantaneous_drain_flux_s'],'max_transient_normalization_error':max(r['norm_error'] for r in rows),'all_transient_probabilities_in_0_1':True,'precision_dps':100}

def main():
 start=time.time();tests={};Ts=m.P['temperatures_K'];Fs=np.round(np.arange(0,3.000001,.005),8)
 base=json.loads((ROOT/'inputs/baseline_ff78.json').read_text())['full_local_parameters']
 assert np.allclose([base['Eg'],base['eps_r'],base['Nt'],base['d']],[1.3,3.5,1e15,1e-5],rtol=1e-14,atol=0)
 assert np.isclose(sum(v*v for v in m.P['couplings_eV'].values()),m.P['sum_squared_couplings_eV2'],rtol=1e-14,atol=0)
 tests['baseline_load_check']=True;tests['fixed_budget_eV2']=2e-10
 # Units checked independently in SI, including the one conversion of eV^-1 DOS.
 HJ=1e-5*m.QE;hbarJ=m.HBAR*m.QE;lambdaJ=.2*m.QE;kTJ=m.KB*300*m.QE
 ksi=2*np.pi*HJ**2/hbarJ/np.sqrt(4*np.pi*lambdaJ*kTJ)
 kev=np.exp(m.logkernel(-.2,300,'classical'));tests['hbar_eV_vs_SI_relative_error']=abs(ksi/kev-1)
 kern=[]
 for T in Ts:
  l,lp=m.sidebands(T);mass=np.exp(logsumexp(lp))
  for dg in np.linspace(-1.3,1.3,105):
   a=float(m.logkernel(dg,T));b=float(m.logkernel(-dg,T));c=float(m.logkernel(dg,T,L=80))
   kern.append(dict(T_K=T,dg_eV=dg,sideband_mass_error=abs(mass-1),DB_log_error=a-b+dg/(m.KB*T),cutoff60_vs80_log_error=a-c))
 tests['kernel_DB_log_error']=max(abs(r['DB_log_error']) for r in kern);tests['sideband_mass_error']=max(r['sideband_mass_error'] for r in kern);tests['cutoff60_vs80_error']=max(abs(r['cutoff60_vs80_log_error']) for r in kern)
 write('kernel_tests.csv',kern)
 fc=[]
 for T in [80,150,300,330]:
  for dg in [-1.2,-.6,-.35,0,.35,.6,1.2]:
   fc.append(dict(T_K=T,dg_eV=dg,double_sum_N64_lograte=m.franck_condon_double_lograte(dg,T),sideband_lograte=float(m.logkernel(dg,T))))
 tests['independent_FC_max_log_error']=max(abs(z['double_sum_N64_lograte']-z['sideband_lograte']) for z in fc);write('independent_franck_condon.csv',fc)
 # All 8 states, charges and energies at representative fields.
 states=[]
 for F in [0,1,2,3]:
  for i,s in enumerate(m.STATES):
   states.append(dict(F_MVcm=F,state=''.join(map(str,s)),electron_number=sum(s),qV_e=1-s[0],qD_e=-s[1],qC_e=-s[2],total_charge_e=1-sum(s),energy_eV=m.energy(F)[i]))
 write('states_and_charges.csv',states)
 eq=[];crows=[];basis=cycle_basis();maxeqedge=0
 for T in Ts:
  for F in [0,.2,1,2,3]:
   for kind in ['classical','quantum']:
    row,p,rates,Q,edges,E,mu=m.solve(F,T,kind,equilibrium=True,detail=True)
    lp=[-(E[i]-mm.mpf('.023')*sum(s))/(mm.mpf(str(m.KB))*T) for i,s in enumerate(m.STATES)]
    Z=sum(mm.exp(x) for x in lp);boltz=[mm.exp(x)/Z for x in lp]
    pr=max(abs(p[i]/boltz[i]-1) for i in range(8));ee=max(abs(f-r)/max(f,r) for i,j,name,f,r,*_ in edges)
    eq.append(dict(T_K=T,F_MVcm=F,kernel=kind,stationary_vs_grandcanonical_relative_error=float(pr),max_edge_relative_flux_error=float(ee),net_cycle_s=row['R_s']))
    for ci,cyc in enumerate(basis):
     aff=sum(sign*mm.log(rates[e][0]/rates[e][1]) for e,sign in cyc)
     crows.append(dict(T_K=T,F_MVcm=F,kernel=kind,equilibrium=True,cycle=ci,computed_affinity=float(aff),predicted_affinity=0.,affinity_error=float(abs(aff))))
    # Driven cycle affinity depends only on bath transfer, never Coulomb residue.
    _,_,rr,_,_,_,mus=m.solve(F,T,kind,detail=True)
    for ci,cyc in enumerate(basis):
     aff=sum(sign*mm.log(rr[e][0]/rr[e][1]) for e,sign in cyc)
     dnV=sum(sign for e,sign in cyc if m.EDGES[e][2]=='V');pred=(mus['V']-mus['C'])*dnV/(mm.mpf(str(m.KB))*T)
     crows.append(dict(T_K=T,F_MVcm=F,kernel=kind,equilibrium=False,cycle=ci,computed_affinity=float(aff),predicted_affinity=float(pred),affinity_error=float(abs(aff-pred))))
 write('equilibrium_nonzero_field.csv',eq);write('all_fundamental_cycle_affinities.csv',crows)
 tests['equilibrium_probability_max_relative_error']=max(r['stationary_vs_grandcanonical_relative_error'] for r in eq);tests['equilibrium_edge_max_relative_error']=max(r['max_edge_relative_flux_error'] for r in eq);tests['cycle_affinity_max_error']=max(r['affinity_error'] for r in crows)
 (R/'cycle_basis.json').write_text(json.dumps([{'cycle':ci,'oriented_edges':[{'from':''.join(map(str,m.STATES[m.EDGES[e][0 if s==1 else 1]])),'to':''.join(map(str,m.STATES[m.EDGES[e][1 if s==1 else 0]])),'process':m.EDGES[e][2]} for e,s in c]} for ci,c in enumerate(basis)],indent=2))
 # Main scan, same all parameters except the classical versus quantum kernel.
 rows=[]
 for T in Ts:
  for kind in ['classical','quantum']:
   for F in Fs:rows.append(m.solve(float(F),T,kind))
   print('scan',T,kind,'elapsed',round(time.time()-start,1),flush=True)
 write('full_reversible_cycle_scan.csv',rows)
 tests['scan_rows']=len(rows);tests['scan_min_probability']=min(r['minimum_probability'] for r in rows);tests['scan_max_continuity_relative_error_nonzero_field']=max(r['continuity_relative_error'] for r in rows if r['F_MVcm']>0)
 tests['scan_max_generator_residual']=max(r['generator_scaled_residual'] for r in rows);tests['scan_max_heat_vs_chemical_error_nonzero_field']=max(r['heat_vs_chemical_relative_error'] for r in rows if r['F_MVcm']>0);tests['scan_min_entropy_production']=min(r['entropy_production_kB_s'] for r in rows)
 # Extraction limitation is separately diagnosed; no baseline parameters retuned.
 ext=[]
 for T in [80,300]:
  for F in [.5,1,1.5,2,3]:
   for kind in ['classical','quantum']:
    for gamma in [1,1e2,1e4,1e6,1e8,1e10,1e14]:ext.append(m.solve(F,T,kind,gamma=gamma))
 write('extraction_rate_diagnostic.csv',ext)
 # Asymmetric finite collection verifies both supplies can block independently.
 asym=[]
 for gv,gc in [(1e2,1e6),(1e6,1e2),(1e2,1e2),(1e6,1e6)]:asym.append(m.solve(1.5,80,'quantum',gamma={'V':gv,'C':gc}))
 write('asymmetric_extraction.csv',asym)
 # Failed implementations, retained as counterexamples.
 bad=[]
 for T in [80,300]:
  for F in [0,.5,1,1.5,2,3]:
   for kind,irr in [('quantum',False),('naive_ground_mlj',False),('quantum',True)]:
    rr=m.solve(F,T,kind,equilibrium=True,naive_irreversible_contacts=irr)
    bad.append(dict(T_K=T,F_MVcm=F,implementation=kind+('_irreversible_contacts' if irr else '_reversible_contacts'),common_mu_eV=.023,spurious_equilibrium_R_s=rr['R_s'],heat_bath_eV_s=rr['heat_bath_eV_s'],chemical_power_eV_s=rr['chemical_power_eV_s'],rate_entropy_production_kB_s=rr['entropy_production_kB_s'] if not irr else '',physical_heat_over_kBT_s=rr['heat_bath_eV_s']/(m.KB*T)))
 write('counterexamples.csv',bad)
 # Independent tree theorem and higher precision; expose ordinary-double errors.
 trees=m.spanning_trees();checks=[]
 for T in [80,300]:
  for F in [0,.05,.5,1,1.5,2,3]:
   for kind in ['classical','quantum']:
    rr,p,rat,Q,*_=m.solve(F,T,kind,detail=True);tp=m.tree_stationary(rat,trees)
    pp=np.array([float(z) for z in p]);treeerr=max(abs(tp/pp-1))
    high,hp,hr,hQ,he,*_=m.solve(F,T,kind,dps=170,detail=True)
    r110=sum(p[i]*r[0]-p[j]*r[1] for (i,j,name),r in zip(m.EDGES,rat) if name=='DC')
    r170=sum(e[5] for e in he if e[2]=='DC');herror=float(abs(r170/r110-1)) if F else 0
    try:
     dp,cond=m.numpy_stationary(Q);j=sum(dp[i]*float(kf)-dp[j]*float(kr) for (i,j,name),(kf,kr,*_) in zip(m.EDGES,rat) if name=='DC')
     derr=abs(j/rr['R_s']-1) if F else abs(j);dmin=float(min(dp));status='solved'
    except np.linalg.LinAlgError:cond=np.inf;derr=np.inf;dmin=np.nan;j=np.nan;status='singular'
    checks.append(dict(T_K=T,F_MVcm=F,kernel=kind,tree_relative_probability_error=treeerr,precision110_vs170_relative_R_error=herror,double_matrix_condition=cond,double_min_probability=dmin,double_R_s=j,validated_R_s=rr['R_s'],double_current_relative_error_if_Fpositive=derr,double_status=status))
 write('precision_and_tree_tests.csv',checks)
 tests['independent_spanning_tree_count']=len(trees);tests['tree_max_probability_relative_error']=max(r['tree_relative_probability_error'] for r in checks);tests['precision110_vs170_max_relative_R_error']=max(r['precision110_vs170_relative_R_error'] for r in checks)
 tests['ordinary_double_negative_probability_cases']=sum(r['double_min_probability']<0 for r in checks);tests['ordinary_double_bad_current_cases']=sum(r['double_current_relative_error_if_Fpositive']>1e-5 for r in checks if r['F_MVcm']>0)
 tests['closed_finite_bath']=finite_closed_reservoir()
 # Frozen-parameter threshold indicators, expressly not sample knee predictions.
 crossings=[];peaks=[]
 for T in Ts:
  for kind in ['classical','quantum']:
   curve=[r for r in rows if r['T_K']==T and r['kernel']==kind];best=max(curve,key=lambda r:r['R_s'])
   peaks.append({k:best[k] for k in ['T_K','kernel','F_MVcm','R_s','conditional_J_Acm2']})
   for threshold in [.0001,.001,.05]:
    found=[]
    for a,b in zip(curve[:-1],curve[1:]):
     aa=a['conditional_J_Acm2']-threshold;bb=b['conditional_J_Acm2']-threshold
     if aa*bb<0:
      ff=a['F_MVcm']+(b['F_MVcm']-a['F_MVcm'])*(-aa)/(bb-aa);found.append((ff,'rising' if bb>aa else 'falling'))
    if not found:crossings.append(dict(T_K=T,kernel=kind,threshold_Acm2=threshold,F_MVcm='',direction='',reason='no crossing within 0 to 3 MV/cm'))
    else:
     for ff,direction in found:crossings.append(dict(T_K=T,kernel=kind,threshold_Acm2=threshold,F_MVcm=ff,direction=direction,reason='linear interpolation; local-field source-capacity diagnostic only'))
 write('conditional_threshold_crossings.csv',crossings);write('scan_peaks.csv',peaks)
 # Field step convergence of peak / threshold tested by decimation, not new tuning.
 conv=[]
 for T in [80,300]:
  for kind in ['classical','quantum']:
   curve=[r for r in rows if r['T_K']==T and r['kernel']==kind]
   for stride in [1,2,4]:
    best=max(curve[::stride],key=lambda r:r['R_s']);conv.append(dict(T_K=T,kernel=kind,field_step_MVcm=.005*stride,peak_F_MVcm=best['F_MVcm'],peak_R_s=best['R_s']))
 write('field_step_convergence.csv',conv)
 # Plots preserve full data, log plot excludes exact F=0 and has disclosed floor.
 fig,ax=plt.subplots(1,2,figsize=(11,4.4),sharey=True)
 for k,a in zip(['classical','quantum'],ax):
  for T in [80,150,220,300,330]:
   cr=[r for r in rows if r['T_K']==T and r['kernel']==k and r['F_MVcm']>0]
   a.semilogy([r['F_MVcm'] for r in cr],[r['R_s'] for r in cr],label=f'{T} K')
  a.set(xlabel='Local favorable field (MV/cm)',ylabel='Net completed cycles per graph (s$^{-1}$)',title=k+'; finite reversible endpoints',ylim=(1e-40,2e6));a.grid(alpha=.2);a.legend(fontsize=8)
 fig.suptitle('Fixed exploratory geometry and coupling budget; dark isothermal benchmark')
 savefig('cycle_vs_field.png')
 fig,ax=plt.subplots(1,2,figsize=(11,4.3))
 for F in [.5,1,1.5,2,3]:
  q=[r for r in rows if r['kernel']=='quantum' and r['F_MVcm']==F];c=[r for r in rows if r['kernel']=='classical' and r['F_MVcm']==F]
  ax[0].semilogy(Ts,[r['R_s'] for r in q],marker='o',label=f'{F:g} MV/cm');ax[1].semilogy(Ts,[a['R_s']/b['R_s'] for a,b in zip(q,c)],marker='o',label=f'{F:g} MV/cm')
 ax[0].set(xlabel='Temperature (K)',ylabel='Quantum net cycles (s$^{-1}$)',title='Cooling response depends on field')
 ax[1].set(xlabel='Temperature (K)',ylabel='Quantum / classical net cycle rate',title='Same parameters and total reorganization')
 for a in ax:a.grid(alpha=.2);a.legend(fontsize=8)
 savefig('temperature_and_quantum_ratio.png')
 fig,ax=plt.subplots(1,2,figsize=(11,4.3))
 for T in [80,300]:
  for F in [1.5,2,3]:
   cr=[r for r in ext if r['kernel']=='quantum' and r['T_K']==T and r['F_MVcm']==F]
   ax[0].loglog([r['GammaC_s'] for r in cr],[r['R_s'] for r in cr],marker='.',label=f'{T} K, {F:g} MV/cm')
 for F in [1.5,2]:
  cr=[r for r in ext if r['kernel']=='quantum' and r['T_K']==80 and r['F_MVcm']==F]
  ax[1].semilogx([r['GammaC_s'] for r in cr],[r['nC'] for r in cr],marker='o',label=f'nC, {F:g} MV/cm')
  ax[1].semilogx([r['GammaC_s'] for r in cr],[1-r['nV'] for r in cr],marker='s',linestyle='--',label=f'hV, {F:g} MV/cm')
 ax[0].set(xlabel='Equal finite exchange rates Gamma (s$^{-1}$)',ylabel='Quantum net cycles (s$^{-1}$)',title='Extraction is part of the complete cycle')
 ax[1].set(xlabel='Equal finite exchange rates Gamma (s$^{-1}$)',ylabel='Tagged endpoint occupancy',title='Blocking at 80 K')
 for a in ax:a.grid(alpha=.2);a.legend(fontsize=8)
 savefig('extraction_blocking.png')
 fig,ax=plt.subplots(figsize=(7,4.3))
 for T in [80,300]:
  cr=[r for r in bad if r['T_K']==T and r['implementation']=='naive_ground_mlj_reversible_contacts'];ax.plot([r['F_MVcm'] for r in cr],[r['spurious_equilibrium_R_s'] for r in cr],marker='o',label=f'ground-only sign flip, {T} K')
 ax.set_yscale('symlog',linthresh=1e-20);ax.set(xlabel='Static field with common mu (MV/cm)',ylabel='False equilibrium cycle current (s$^{-1}$)',title='Counterexample: equilibrium current must be zero');ax.legend(fontsize=8);ax.grid(alpha=.2);savefig('failed_naive_reverse.png')
 # Hard checks cover mathematics, not material applicability.
 assert tests['kernel_DB_log_error']<1e-11 and tests['sideband_mass_error']<1e-12
 assert tests['cutoff60_vs80_error']<1e-12 and tests['independent_FC_max_log_error']<1e-11
 assert tests['equilibrium_probability_max_relative_error']<1e-70 and tests['equilibrium_edge_max_relative_error']<1e-70
 assert tests['scan_min_probability']>0 and tests['scan_max_continuity_relative_error_nonzero_field']<1e-70
 assert tests['scan_max_heat_vs_chemical_error_nonzero_field']<1e-70 and tests['precision110_vs170_max_relative_R_error']<1e-12
 assert tests['tree_max_probability_relative_error']<1e-10
 # Explicit complete representative cycle, including the Coulomb reset steps.
 ledger=[];seq=[(1,0,0),(0,1,0),(0,0,1),(0,0,0),(1,0,0)]
 for F in [0,1.5,3]:
  E=m.energy(F);mu={'V':F/2,'C':-F/2};names=['VD','DC','C out','V in']
  for s1,s2,name in zip(seq[:-1],seq[1:],names):
   i=m.INDEX[s1];j=m.INDEX[s2];dn=sum(s2)-sum(s1);chem=(mu['V'] if name=='V in' else mu['C'] if name=='C out' else 0)*dn
   dfield=sum(-(s2[z]-s1[z])*.1*F*m.P['positions_nm'][z] for z in range(3));dbare=sum((s2[z]-s1[z])*m.P['bare_electron_energies_eV'][z] for z in range(3));de=E[j]-E[i]
   ledger.append(dict(F_MVcm=F,from_state=''.join(map(str,s1)),to_state=''.join(map(str,s2)),process=name,delta_energy_eV=de,delta_bare_eV=dbare,delta_field_eV=dfield,delta_coulomb_eV=de-dfield-dbare,chemical_work_eV=chem,delta_thermodynamic_eV=de-chem,heat_to_bath_eV=chem-de))
 write('representative_cycle_energy_ledger.csv',ledger)
 tests['every_internal_edge_conserves_electrons']=all(sum(m.STATES[i])==sum(m.STATES[j]) for i,j,n in m.EDGES if n in ['VD','DC'])
 tests['every_reservoir_addition_has_one_electron']=all(sum(m.STATES[j])-sum(m.STATES[i])==1 for i,j,n in m.EDGES if n in ['V','C'])
 tests['charge_balance_identity']=all((1-sum(s))+sum(s)==1 for s in m.STATES)
 tests['elapsed_s']=time.time()-start;tests['all_assertions_passed']=True
 (R/'validation_summary.json').write_text(json.dumps(tests,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x)))
 provenance={'created_UTC':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'scipy':scipy.__version__,'mpmath':mm.__version__,'matplotlib':matplotlib.__version__,'command':'python -B run_audit.py','scope':'new reconstruction, not recovered historical graph; isothermal exploratory local source graph','sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'parameters.json',ROOT/'model.py',ROOT/'run_audit.py',ROOT/'inputs/baseline_ff78.json',ROOT/'inputs/temperature_audit_source.py']}}
 (R/'provenance.json').write_text(json.dumps(provenance,indent=2));print(json.dumps(tests,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x)),flush=True)
if __name__=='__main__':main()
