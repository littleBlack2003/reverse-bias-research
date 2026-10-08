from pathlib import Path
import json,csv,sys,hashlib,platform,itertools
import numpy as np
import mpmath as mp
import comparator_model as m
ROOT=Path(__file__).resolve().parent;R=ROOT/'results'
def csvout(name,rows):
 with (R/name).open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)

def checks_for(rows):
 return {k:max(abs(x[k] or 0) for x in rows) for k in ['normalization_error','generator_residual','continuity_relative','heat_vs_power_relative','heat_partition_relative','entropy_balance_relative','displacement_relative','direct_DB_max_log_error','conditional_spectrum_mean_error_eV']}

def main():
 rows=[];edge=[];tests={}
 Fs=[0,.5,1.,1.2,1.3,1.5,2.,3.]
 for T in [80,300]:
  for topology in ['direct','trap']:
   for F in Fs:
    row,p,Q,rates,ee,E,states,edges=m.solve(F,T,topology,detail=True)
    rows.append(row)
    if F in [1.,1.5]:edge.extend(dict(T_K=T,F_MVcm=F,topology=topology,**x) for x in ee)
 csvout('fixed_budget_comparison.csv',rows);csvout('representative_edge_ledger.csv',edge)
 tests['main']=checks_for(rows)
 # Fixed physical prefactor law without restoring unequal budgets.
 dist=[]
 for T in [80,300]:
  for F in [1.,1.5,2.]:
   for a in [.5,1.,2.]:
    for topology in ['direct','trap']:dist.append(m.solve(F,T,topology,mode='distance',decay=a))
 csvout('distance_law_comparison.csv',dist);tests['distance']=checks_for(dist)
 # Common mu at nonzero local field; gauge change; independent cycle drive.
 eq=[];ctrl=[];treeerrors=[];referrs=[]
 for T in [80,300]:
  for F in [.5,1.5,3.]:
   for topology in ['direct','trap']:
    eq.append(m.solve(F,T,topology,delta_mu=0,mean_mu=.023))
 for T in [80,300]:
  for topology in ['direct','trap']:
   base=m.solve(1.5,T,topology);gauge=m.solve(1.5,T,topology,gauge=.371)
   ctrl.append(dict(test='gauge_shift_0.371eV',T_K=T,topology=topology,relative_rate_error=abs(gauge['R_pair_s']/base['R_pair_s']-1),baseline_rate_s=base['R_pair_s'],control_rate_s=gauge['R_pair_s']))
   off=m.solve(1.5,T,topology,bath_scale=(1.,0.))
   ctrl.append(dict(test='remove_right_extraction_reservoir',T_K=T,topology=topology,relative_rate_error=abs(off['R_pair_s'])/base['R_pair_s'],baseline_rate_s=base['R_pair_s'],control_rate_s=off['R_pair_s']))
   for order,dps in [(48,220),(24,320)]:
    a=m.solve(1.5,T,topology,order=order,dps=dps)
    ctrl.append(dict(test=f'order_{order}_dps_{dps}',T_K=T,topology=topology,relative_rate_error=abs(a['R_pair_s']/base['R_pair_s']-1),baseline_rate_s=base['R_pair_s'],control_rate_s=a['R_pair_s']))
  for F in [0,1.,1.5,3.]:
   _,p,Q,*_=m.solve(F,T,'direct',detail=True)
   pt=m.independent_tree_probabilities(Q);treeerrors.append(float(max(abs(a/b-1) for a,b in zip(pt,p))))
  ref=m.old.solve(1.5,T,dps=220)['R_s'];new=m.solve(1.5,T,'trap')['R_pair_s'];referrs.append(abs(new/ref-1))
 csvout('nonzero_field_equilibrium.csv',eq);csvout('focused_controls.csv',ctrl)
 tests['equilibrium_max_probability_relative']=max(x['equilibrium_probability_relative'] for x in eq)
 tests['equilibrium_max_edge_net_over_gross']=max(x['equilibrium_edge_net_over_gross'] for x in eq)
 tests['independent_four_state_tree_max_relative']=max(treeerrors)
 tests['prior_three_site_model_rate_max_relative_difference']=max(referrs)
 tests['focused_control_max_relative']=max(x['relative_rate_error'] for x in ctrl)
 # Binding and finite reservoir supply controls.
 extra=[]
 for T in [80,300]:
  for topology in ['direct','trap']:
   extra.append(m.solve(1.5,T,topology,epsr=None))
   for ratefactor in [.01,1e-6]:extra.append(m.solve(1.5,T,topology,bath_scale=(ratefactor,ratefactor)))
 csvout('binding_and_supply_controls.csv',extra);tests['extra']=checks_for(extra)
 # State energy and common-mode displacement Gram construction.
 stateledger=[]
 for topology in ['direct','trap']:
  row,p,Q,rates,ee,E,states,edges=m.solve(1.5,80,topology,detail=True)
  for i,s in enumerate(states):stateledger.append(dict(topology=topology,state=''.join(map(str,s)),electron_count=sum(s),charge_e=1-sum(s),energy_eV=float(E[i]),probability=float(p[i])))
 csvout('state_charge_energy_ledger.csv',stateledger)
 grams={}
 for n in [2,3]:
  # One common plus independent local displacement per site, each carrying half
  # each sector's one-electron reorganization. Surface minima energies explicit.
  A=np.zeros((n,n+1));A[:,0]=np.sqrt(.5)
  for i in range(n):A[i,i+1]=np.sqrt(.5)
  G=A@A.T
  grams[str(n)]={'site_displacement_gram_normalized':G.tolist(),'site_addition_norm2':np.diag(G).tolist(),'all_transfer_norm2':[float(np.dot(A[i]-A[j],A[i]-A[j])) for i in range(n) for j in range(i+1,n)],'sector_reorganization_eV':{'slow':.05,'high':.15},'scope':'a possible shared harmonic bath; not a measured mode structure'}
  assert np.allclose(np.diag(G),1) and all(abs(v-1)<1e-14 for v in grams[str(n)]['all_transfer_norm2'])
 (R/'common_bath_construction.json').write_text(json.dumps(grams,indent=2))
 for batch in ['main','distance','extra']:
  for k,v in tests[batch].items():assert v<(1e-9 if k in ['direct_DB_max_log_error','conditional_spectrum_mean_error_eV'] else 1e-45),(batch,k,v)
 assert tests['equilibrium_max_probability_relative']<1e-45
 assert tests['equilibrium_max_edge_net_over_gross']<1e-45
 assert tests['independent_four_state_tree_max_relative']<1e-100
 assert tests['prior_three_site_model_rate_max_relative_difference']<1e-10
 assert tests['focused_control_max_relative']<1e-8
 tests['all_assertions_passed']=True
 tests['python']=platform.python_version();tests['numpy']=np.__version__;tests['mpmath']=mp.__version__
 (R/'validation.json').write_text(json.dumps(tests,indent=2));print(json.dumps(tests,indent=2))
if __name__=='__main__':main()
