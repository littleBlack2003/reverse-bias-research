"""Numerical confidence probe only. Classical rates are not validated at low T."""
from pathlib import Path
import sys,json,time,traceback,numpy as np
ROOT=Path(__file__).resolve().parents[1]
REFERENCE=ROOT/'reference_solver'
sys.path.insert(0,str(REFERENCE/'code'))
from qf_solver import QFDevice,PFParams,BASE,qsolve,audit,ControlledPF
ROOT=Path(__file__).resolve().parents[1]
def save(path,obj):path.write_text(json.dumps(obj,indent=2,default=float,allow_nan=False))
def detailed(d,z,V):
 o=d.evaluate(z,V);y=o['y'];ch,r=__import__('rate_audit').channels(d,y);En=np.array([sum(x[2] for x in c) for c in ch['n']]);Ep=np.array([sum(x[2] for x in c) for c in ch['p']]);An=np.array([sum(x[1] for x in c) for c in ch['n']]);Ap=np.array([sum(x[1] for x in c) for c in ch['p']]);S=An+Ap+En+Ep
 pair=(En*Ep-An*Ap)/S;pairA=d.q*d.p.Nt*d.h*d.p.d*sum(pair[1:-1]);bimA=d.q*d.h*d.p.d*sum(o['bim'][1:-1]);phys=d.physical(z,V);nsrc=phys['electron_source_Acm2'];psrc=phys['hole_source_Acm2'];expected=pairA-bimA
 full_error=max(abs(phys['electron_out_Acm2']-expected),abs(phys['hole_out_Acm2']-expected));den=max(abs(expected),1e-300)
 J=d.jacobian(z,V,0);sc=np.asarray(abs(J).max(axis=1).toarray()).ravel();rates=np.concatenate([En[1:-1],Ep[1:-1],An[1:-1],Ap[1:-1]])
 affinity=[];balanced=[]
 for carrier,t,end,inc,em,energy,barrier in r[5]:
  conc=np.exp(y[end,1 if carrier=='n' else 2])*d.n0;valid=(inc>0)&(em>0)
  if np.any(valid):affinity.extend(abs(np.log(em[valid]/inc[valid])+np.log(conc[valid]/d.p.Nc)+energy[valid]/d.vt).tolist())
  f,fb=r[2],r[-1];u,v=(inc*fb[t],em*f[t]) if carrier=='n' else(inc*f[t],em*fb[t]);balanced.extend((abs(u-v)/np.maximum(u+v,1e-300)).tolist())
 return dict(V=V,T=d.p.temperature,depth=d.trap_e_depth,nodes=d.nodes,actual_distance_nm=d.actual_escape_nm,**phys,pair_cycle_net_Acm2=float(pairA),pair_minus_bim_Acm2=float(expected),independent_cycle_error_Acm2=float(full_error),relative_cycle_error=float(full_error/den),source_n_error_vs_cycles_Acm2=float(nsrc-expected),source_p_error_vs_cycles_Acm2=float(psrc-expected),min_positive_rate_s=float(min(rates[rates>0])),zero_aggregate_rate_count=int(sum(rates==0)),min_row_magnitude=float(min(sc)),rows_below_1e20=int(sum(sc<1e-20)),endpoint_log_ratio_error=float(max(affinity,default=0)),equilibrium_max_edge_relative_imbalance=float(max(balanced,default=0)) if V==0 else None,finite=bool(np.all(np.isfinite(z)) and np.all(np.isfinite(rates))),physical_validity='Classical fixed-parameter rate model; no low-temperature experimental validity asserted')
def run(T,depth,nodes=81):
 name=f'T{T:g}_D{depth:g}_N{nodes}';kw=dict(BASE);kw['temperature']=T;d=QFDevice(PFParams(**kw),nodes,trap_e_depth=depth)
 z=np.c_[d.x*d.vbi/d.vt,np.zeros(nodes),np.zeros(nodes)];trace=[];rows=[];failure=None;start=time.monotonic()
 try:
  z,info=qsolve(d,0,z,gate=True);rows.append(dict(**detailed(d,z,0),solver_status=info.get('status','converged'),scaled_residual=info['scaled_residual']));np.savez_compressed(ROOT/'data'/f'{name}_0.npz',z=z,y=d.y_from_z(z,0),V=0)
  if info.get('status'):raise RuntimeError('Initial equilibrium '+info['status'])
  prev=0.
  for target in [-5.,-20.]:
   lastV=prev
   for V in np.linspace(prev,target,int(np.ceil(abs(target-prev)/.25))+1)[1:]:
    shift=-(float(V)-lastV)/d.vt;z[:,0]+=shift*d.x;z[:,1]+=shift*(1-d.x);z[:,2]+=shift*d.x
    z,info=qsolve(d,float(V),z,gate=True);trace.append(dict(V=float(V),iterations=info['iterations'],scaled_residual=info['scaled_residual'],gate_passed=info['gate_passed'],status=info.get('status','converged')));lastV=float(V)
    if info.get('status'):
     np.savez_compressed(ROOT/'data'/f'{name}_failed.npz',z=z,y=d.y_from_z(z,V),V=V);save(ROOT/'data'/f'{name}_failed.json',dict(details=detailed(d,z,V),iteration_history=info['history']));raise RuntimeError(f'Continuation unresolved at {V}: '+info['status'])
   prev=target;row=dict(**detailed(d,z,target),solver_status='converged',scaled_residual=info['scaled_residual']);rows.append(row);np.savez_compressed(ROOT/'data'/f'{name}_{int(target)}.npz',z=z,y=d.y_from_z(z,target),V=target);print(name,row,flush=True)
 except Exception as exc:failure=dict(error=str(exc),traceback=traceback.format_exc());print(name,failure['error'],flush=True)
 save(ROOT/'data'/f'{name}.json',dict(name=name,T=T,depth=depth,nodes=nodes,rows=rows,trace=trace,failure=failure,seconds=time.monotonic()-start))
if __name__=='__main__':run(float(sys.argv[1]),float(sys.argv[2]),int(sys.argv[3]) if len(sys.argv)>3 else 81)
