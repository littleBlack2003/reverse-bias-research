from pathlib import Path
import json,csv,sys,platform,hashlib,time
import numpy as np,scipy,mpmath as mp
import spatial_model as sm
ROOT=Path(__file__).parent;OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
def writecsv(name,rows):
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with (OUT/name).open('w',newline='',encoding='utf8') as f:
  w=csv.DictWriter(f,keys);w.writeheader();w.writerows(rows)
def compare(a,b):return abs(a-b)/max(abs(a),abs(b),1e-250)
def schur(F,T,s=6.):
 full=sm.solve(F,T,s=s,detail=True);row,p,Q,rates,details,E,mu,g,states,edges=full
 keep=[i for i,st in enumerate(states) if st[0]==1 and st[-1]==0];elim=[i for i in range(len(states)) if i not in keep]
 A=mp.matrix([[Q[i,j] for j in keep] for i in keep]);BB=mp.matrix([[Q[i,j] for j in elim] for i in elim]);AB=mp.matrix([[Q[i,j] for j in elim] for i in keep]);BA=mp.matrix([[Q[i,j] for j in keep] for i in elim])
 X=mp.matrix(len(elim),len(keep))
 for k in range(len(keep)):
  sol=mp.lu_solve(BB,-BA[:,k]);X[:,k]=sol
 eff=A+AB*X
 trace=sm.core.stationary(eff);weights=[1+sum(X[:,k]) for k in range(len(keep))];norm=sum(weights[k]*trace[k] for k in range(len(keep)));pA=trace/norm;pB=X*pA
 rec=mp.matrix([0]*len(states))
 for k,i in enumerate(keep):rec[i]=pA[k]
 for k,i in enumerate(elim):rec[i]=pB[k]
 R=mp.mpf(0);Rwrong=mp.mpf(0)
 for (i,j,ty,k),r in zip(edges,rates):
  if r['name']=='V_D':
   R+=rec[i]*r['kf']-rec[j]*r['kr'];Rwrong+=(rec[i]*r['kf']-rec[j]*r['kr'])*norm
 # Strong lumpability: outgoing aggregate rates to every other core-occupancy block must agree.
 groups={}
 for i,st in enumerate(states):groups.setdefault(st[1:4],[]).append(i)
 violation=mp.mpf(0);example=None
 for src,inds in groups.items():
  for dst,targets in groups.items():
   if dst==src:continue
   vals=[sum(Q[t,i] for t in targets) for i in inds];diff=max(vals)-min(vals)
   if diff>violation:violation=diff;example={'source_core':src,'destination_core':dst,'aggregate_rates_s':[float(x) for x in vals],'hidden_states':[''.join(map(str,states[i])) for i in inds]}
 return dict(F_MVcm=F,T_K=T,s_nm=s,R_full_s=row['R_s'],R_reconstructed_s=float(R),probability_max_relative=float(max(abs(rec[i]/p[i]-1) for i in range(len(p)))),current_relative=float(abs(R-mp.mpf(str(row['R_s'])))/max(abs(R),mp.mpf('1e-200'))),retained_probability_mass=float(sum(p[i] for i in keep)),trace_clock_factor=float(norm),naive_trace_clock_R_s=float(Rwrong),effective_generator_colsum=float(max(abs(sum(eff[:,j])) for j in range(eff.cols))),lumpability_violation_s=float(violation),lumpability_example=example)

def main():
 start=time.time();allrows=[];main=[]
 for T in [80,150,300]:
  for F in [.5,1.,1.5,2.,3.]:
   for v in ['colocated','direct','left_only','right_only','explicit']:
    r=sm.solve(F,T,v);main.append(r);allrows.append(r)
 writecsv('controlled_cases.csv',main);print('main',len(main),time.time()-start,flush=True)
 tags=[]
 for T in [80,300]:
  for F in [1.,1.5,2.]:
   for s in [5.25,5.5,5.75,5.9,6.]:
    r=sm.solve(F,T,s=s);tags.append(r);allrows.append(r)
 writecsv('last_tag_placement.csv',tags)
 budgets=[]
 for T in [80,300]:
  for eta in [.001,.003,.01,.03,.1,.25,.5,.75,.9,.99]:
   r=sm.solve(1.5,T,eta=eta);budgets.append(r);allrows.append(r)
 writecsv('fixed_budget_allocation.csv',budgets)
 decay=[]
 for T in [80,300]:
  for xi in [.2,.5,1.,2.]:
   for v in ['direct','explicit']:
    r=sm.solve(1.5,T,v,decay_nm=xi);decay.append(r);allrows.append(r)
 writecsv('distance_attenuation_controls.csv',decay)
 noC=[]
 for T in [80,300]:
  for F in [1.,1.5,2.,3.]:
   for v in ['colocated','direct','explicit']:
    r=sm.solve(F,T,v,epsr=None);noC.append(r);allrows.append(r)
 writecsv('zero_coulomb_controls.csv',noC)
 drives=[]
 for T in [80,300]:
  for F in [0.,1.5]:
   for dm in [-1.5,0.,1.5,1.8]:
    for v in ['direct','explicit']:
     r=sm.solve(F,T,v,delta_mu=dm,mean_mu=.023);drives.append(r);allrows.append(r)
 writecsv('independent_drive_equilibrium.csv',drives)
 # Eight selected states get independent gauge/precision/quadrature checks; no broad unverifiable claim.
 check=[]
 for T,F,v,s in [(80,.5,'explicit',6.),(80,1.5,'explicit',6.),(80,3.,'explicit',6.),(300,1.5,'explicit',6.),(80,1.5,'direct',6.),(300,1.5,'direct',6.),(80,1.5,'explicit',5.5),(300,1.5,'explicit',5.5),(80,1.,'explicit',5.25),(80,1.,'explicit',5.5)]:
  a=sm.solve(F,T,v,s=s,detail=True);b=sm.solve(F,T,v,s=s,gauge=.417,detail=True);c=sm.solve(F,T,v,s=s,order=48,detail=True);d=sm.solve(F,T,v,s=s,dps=350,detail=True)
  check.append(dict(T_K=T,F_MVcm=F,variant=v,s_nm=s,gauge_R_relative=compare(a[0]['R_s'],b[0]['R_s']),gauge_probability_relative=float(max(abs(a[1][i]/b[1][i]-1) for i in range(len(a[1])))),gauge_heat_relative=compare(a[0]['total_bath_heat_eV_s'],b[0]['total_bath_heat_eV_s']),quadrature_R_relative=compare(a[0]['R_s'],c[0]['R_s']),quadrature_probability_relative=float(max(abs(a[1][i]/c[1][i]-1) for i in range(len(a[1])))),precision_R_relative=compare(a[0]['R_s'],d[0]['R_s']),precision_probability_relative=float(max(abs(a[1][i]/d[1][i]-1) for i in range(len(a[1]))))))
 writecsv('gauge_quadrature_precision.csv',check)
 sh=[schur(F,T,s) for F,T,s in [(1.5,80,6.),(.5,80,6.),(1.5,300,6.),(1.5,80,5.5)]]
 (OUT/'exact_steady_elimination.json').write_text(json.dumps(sh,indent=2))
 edgeRows=[];stateRows=[]
 for T in [80,300]:
  for v,s in [('colocated',6.),('direct',6.),('explicit',6.),('explicit',5.5),('left_only',6.),('right_only',6.)]:
   r,p,Q,rt,de,E,mu,g,states,edges=sm.solve(1.5,T,v,s=s,detail=True)
   for z in de:edgeRows.append(dict(T_K=T,F_MVcm=1.5,variant=v,s_nm=s,**z))
   for i,st in enumerate(states):stateRows.append(dict(T_K=T,F_MVcm=1.5,variant=v,s_nm=s,state=''.join(map(str,st)),probability=float(p[i]),energy_eV=float(E[i])))
 writecsv('representative_edges.csv',edgeRows);writecsv('representative_states.csv',stateRows)
 validation=dict(number_main_and_control_points=len(allrows),max_node_continuity_relative=max(r['max_node_continuity_relative'] for r in allrows),max_heat_chemical_relative=max(r['heat_chemical_relative'] for r in allrows),max_heat_partition_relative=max(r['heat_partition_relative'] for r in allrows),max_entropy_heat_relative=max(r['entropy_heat_relative'] for r in allrows),max_displacement_relative=max(r['displacement_relative'] for r in allrows),min_probability=min(r['min_probability'] for r in allrows),max_direct_DB_log_error=max(r['direct_DB_max_log_error'] for r in allrows),max_direct_mean_energy_error_eV=max(r['direct_mean_energy_error_eV'] for r in allrows),max_common_mu_probability_relative=max(r['common_mu_probability_relative'] or 0 for r in allrows),max_common_mu_edge_flux_over_gross=max(r['common_mu_edge_flux_over_gross'] or 0 for r in allrows),max_gauge_probability_relative=max(r['gauge_probability_relative'] for r in check),max_quadrature_R_relative=max(r['quadrature_R_relative'] for r in check),max_precision_probability_relative=max(r['precision_probability_relative'] for r in check),max_elimination_probability_relative=max(r['probability_max_relative'] for r in sh),runtime_s=time.time()-start)
 for key in ['max_node_continuity_relative','max_heat_chemical_relative','max_heat_partition_relative','max_displacement_relative','max_common_mu_probability_relative']:
  assert validation[key]<1e-65,(key,validation[key])
 assert validation['max_direct_DB_log_error']<1e-10
 assert validation['max_quadrature_R_relative']<1e-9
 assert validation['min_probability']>0
 (OUT/'validation_summary.json').write_text(json.dumps(validation,indent=2));print(json.dumps(validation,indent=2),flush=True)
 prov=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__,platform=platform.platform(),input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'inputs').rglob('*') if p.is_file()},code_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'spatial_model.py',ROOT/'run_audit.py']})
 (OUT/'provenance.json').write_text(json.dumps(prov,indent=2))
if __name__=='__main__':main()
