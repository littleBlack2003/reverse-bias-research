"""Room-temperature continuation; reuse existing 300 K results before new solves.
Only a one-dimensional field sweep with prior fixed physical parameters.
Normalized F10/F50/F90 characterize a model rising branch, never device V50.
"""
from pathlib import Path
import sys,csv,json,hashlib
import numpy as np
from scipy.optimize import minimize_scalar,brentq
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent
sys.path.insert(0,str(BASE))
import comparator_model as m
SCENARIOS={'direct_fixed':('direct','fixed_budget',None),'trap_fixed':('trap','fixed_budget',None),'direct_a05':('direct','distance',.5),'direct_a10':('direct','distance',1.),'direct_a20':('direct','distance',2.)}
CACHE={};OLDKEYS=set();NEWS=[]
def key(s,F):return s,round(float(F),11)
def save(name,rows):
 with (ROOT/'results'/name).open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)

def loadold():
 for name in ['fixed_budget_comparison.csv','distance_law_comparison.csv']:
  for x in csv.DictReader((BASE/'results'/name).open()):
   if x['T_K']!='300':continue
   if x['mode']=='fixed_budget':s='direct_fixed' if x['topology']=='direct' else 'trap_fixed'
   elif x['topology']=='direct':s={.5:'direct_a05',1.:'direct_a10',2.:'direct_a20'}[float(x['decay_nm'])]
   else:continue
   F=float(x['F_MVcm']);r={}
   for k,v in x.items():
    if v=='':r[k]=None
    else:
     try:r[k]=float(v)
     except ValueError:r[k]=v
   r['reuse_source']=name;r['scenario']=s
   CACHE[key(s,F)]=r;OLDKEYS.add(key(s,F))

def get(s,F):
 k=key(s,F)
 if k not in CACHE:
  top,mode,decay=SCENARIOS[s];r=m.solve(k[1],300,top,mode=mode,decay=decay)
  r['reuse_source']='new_300K_solve';r['scenario']=s;CACHE[k]=r;NEWS.append(r)
 return CACHE[k]

def main():
 loadold();rows=[];summ=[];crossing=[];checks=[]
 fieldgrid=np.round(np.arange(.5,3.00001,.1),10)
 for s in SCENARIOS:
  vals=[get(s,F) for F in fieldgrid]
  for r in vals:rows.append(r)
  ys=np.array([r['R_pair_s'] for r in vals]);i=int(np.argmax(ys))
  # Refine every local-maximum bracket, rather than assuming the first is global.
  candidates=[(float(fieldgrid[0]),float(ys[0])),(float(fieldgrid[-1]),float(ys[-1]))]
  for j in range(1,len(ys)-1):
   if ys[j]>=ys[j-1] and ys[j]>=ys[j+1]:
    opt=minimize_scalar(lambda F:-get(s,F)['R_pair_s'],bounds=(fieldgrid[j-1],fieldgrid[j+1]),method='bounded',options={'xatol':2e-7})
    candidates.append((float(opt.x),-float(opt.fun)))
  Fpk,Rpk=max(candidates,key=lambda p:p[1]);row=dict(scenario=s,T_K=300,F_peak_within_window_MVcm=Fpk,R_peak_s=Rpk,J_peak_for_Npair_1e15_mAcm2=Rpk*m.QELEC*1e10*1000,required_Npair_for_50mA_cm3=.05/(m.QELEC*1e-5*Rpk),reference_window_MVcm='0.5..3.0',device_V50_V=None,V50_status='not reached at stated Npair; dark graph is not illuminated device JV')
  for fraction,label in [(.1,'10'),(.5,'50'),(.9,'90')]:
   target=fraction*Rpk
   # First ascending crossing in fixed window, not a later crossing after inversion.
   bracket=None
   for a,b in zip(fieldgrid[:-1],fieldgrid[1:]):
    if get(s,a)['R_pair_s']<target<=get(s,b)['R_pair_s']:
     bracket=(a,b);break
   assert bracket is not None,(s,fraction)
   Fq=brentq(lambda F:get(s,F)['R_pair_s']/Rpk-fraction,*bracket,xtol=1e-8,rtol=1e-12)
   rq=get(s,Fq)['R_pair_s'];row['F_rise_'+label+'_MVcm']=Fq
   crossing.append(dict(scenario=s,fraction=fraction,F_MVcm=Fq,R_s=rq,target_s=target,relative_target_error=abs(rq/target-1),metric_warning='fraction of each model own bounded-window maximum; not a common current threshold'))
  row['rise_width_F90_minus_F10_MVcm']=row['F_rise_90_MVcm']-row['F_rise_10_MVcm']
  summ.append(row)
  # Check peak and middle crossing against higher quadrature/precision.
  top,mode,decay=SCENARIOS[s]
  for F in [Fpk,row['F_rise_50_MVcm']]:
   base=get(s,F)
   for order,dps in [(48,220),(24,320)]:
    test=m.solve(F,300,top,mode=mode,decay=decay,order=order,dps=dps)
    checks.append(dict(scenario=s,F_MVcm=F,order=order,dps=dps,rate_relative_error=abs(test['R_pair_s']/base['R_pair_s']-1)))
 save('field_curves_300K.csv',rows);save('rise_metrics_300K.csv',summ);save('crossing_checks.csv',crossing);save('precision_checks.csv',checks)
 # Preserve every adaptive field evaluation, including non-grid evaluations.
 save('all_evaluated_300K.csv',list(CACHE.values()))
 maxima={k:max(abs(r[k] or 0) for r in CACHE.values()) for k in ['normalization_error','generator_residual','continuity_relative','heat_vs_power_relative','heat_partition_relative','entropy_balance_relative','displacement_relative','direct_DB_max_log_error']}
 assert max(x['relative_target_error'] for x in crossing)<1e-6
 assert max(x['rate_relative_error'] for x in checks)<1e-7
 assert maxima['direct_DB_max_log_error']<1e-9
 assert all(maxima[k]<1e-45 for k in maxima if k!='direct_DB_max_log_error')
 result=dict(temperature_K=300,reused_unique_points=len(OLDKEYS),new_unique_field_evaluations=len(NEWS),grid_step_MVcm=.1,field_window_MVcm=[.5,3.],peak_search='local maximum brackets from fixed grid, bounded refinement, largest refined result in window; not a proof of global maximum for arbitrary untested spectra',max_crossing_relative_error=max(x['relative_target_error'] for x in crossing),max_refined_point_precision_error=max(x['rate_relative_error'] for x in checks),conservation_and_balance=maxima,all_assertions_passed=True,all_parameters_unchanged=True)
 (ROOT/'results/validation_300K.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
 for x in summ:print(x)
if __name__=='__main__':main()
