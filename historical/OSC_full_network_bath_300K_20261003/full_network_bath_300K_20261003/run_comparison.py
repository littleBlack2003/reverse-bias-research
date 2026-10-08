#!/usr/bin/env python3
"""One fixed physical point; numerical refinement and common-mu validation only.
Copied source kernels are read unchanged. No fit, clipping or rate renormalization.
"""
from pathlib import Path
import sys, json, time, hashlib, csv, platform
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import logsumexp
from scipy.integrate import quad
import mpmath as mp
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'inputs'))
import spatial_model as sm
import kernel_benchmark as kb
OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
original=sm.core.logkernel

def savecsv(name,rows):
 with (OUT/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

class FullKernel:
 def __init__(self,step,L,tol):
  self.L=L;self.tol=tol
  x=np.linspace(-9,9,round(18/step)+1)
  y=np.array([kb.ohmic_logp(float(e),300.,.02,.05,tol) for e in x])
  self.x=x;self.y=y;self.spline=CubicSpline(x,y,extrapolate=False)
  ls,self.lp=kb.thermal_weights(300.,.15,1.,L);self.energy=ls*.15
  # Off-grid checks include tails. Direct Fourier integration, both signs independently.
  checks=np.linspace(-8.85,8.85,121)+step*.317
  self.interp_error=max(abs(self.spline(e)-kb.ohmic_logp(float(e),300.,.02,.05,tol)) for e in checks)
 def __call__(self,dg,T,kind='quantum',H=1e-5,L=60):
  assert T==300
  dg=np.asarray(dg);e=-dg[...,None]-self.energy
  assert np.max(e)<=9 and np.min(e)>=-9
  ans=np.log(2*np.pi/kb.HBAR)+2*np.log(H)+logsumexp(self.lp+self.spline(e),axis=-1)
  assert np.all(np.isfinite(ans))
  return ans

def evaluate(label,order,dps,eq=False):
 sm.old.exchange.cache_clear()
 z=sm.solve(1.5,300,order=order,dps=dps,detail=True,delta_mu=0 if eq else None,mean_mu=.023 if eq else 0)
 row,p,Q,rates,details,E,mu,g,states,edges=z
 row['case']=label;row['quadrature_order']=order
 row['maximum_exit_rate_s']=float(max(-Q[i,i] for i in range(Q.rows)))
 # Independently assembled rates: do not enforce detailed balance for this extra check.
 D=mp.zeros(len(states));raw=[]
 for (i,j,ty,k),r in zip(edges,rates):
  if ty=='hop':
   h=np.sqrt(g['bondH2'][k]);lf=float(sm.core.logkernel(float(r['a']),300,H=h));lr=float(sm.core.logkernel(-float(r['a']),300,H=h))
  else:
   fac=np.log(g['bathH2'][k]/sm.B);lf=r['diag']['log_add']+fac;lr=r['diag']['log_out']+fac
  f=mp.exp(mp.mpf(str(lf)));b=mp.exp(mp.mpf(str(lr)));raw.append((f,b))
  D[j,i]+=f;D[i,j]+=b;D[i,i]-=f;D[j,j]-=b
 directp=sm.core.stationary(D)
 directR=sum(directp[i]*f-directp[j]*b for (i,j,ty,k),(f,b),r in zip(edges,raw,rates) if r['name']=='V_D')
 row['direct_assembly_R_s']=float(directR)
 row['direct_assembly_probability_max_abs_delta']=float(max(abs(a-b) for a,b in zip(p,directp)))
 row['direct_assembly_R_difference_s']=float(directR-mp.mpf(str(row['R_s'])))
 row['direct_assembly_eq_max_edge_flux_over_gross']=float(max(abs(directp[i]*f-directp[j]*b)/(directp[i]*f+directp[j]*b) for (i,j,ty,k),(f,b) in zip(edges,raw))) if eq else None
 st=[]
 for i,s in enumerate(states):st.append(dict(state=''.join(map(str,s)),probability=float(p[i]),energy_eV=float(E[i]),exit_rate_s=float(-Q[i,i]),mean_wait_s=float(-1/Q[i,i])))
 for d,(f,b),r in zip(details,raw,rates):
  d['raw_forward_rate_s']=float(f);d['raw_reverse_rate_s']=float(b)
  d['raw_DB_log_error']=float(mp.log(f/b)+r['tdg']/(mp.mpf(str(sm.KB))*300))
  d['mean_exchange_energy_eV']=float(r['mean']) if r['mean'] is not None else ''
 savecsv(label+'_states.csv',st);savecsv(label+'_edges.csv',details)
 print(label,row['R_s'],row['direct_DB_max_log_error'],flush=True)
 return row

def main():
 rows=[];validation={}
 for label,step,L,tol,order,dps in [('main',.002,30,2e-10,24,250),('refined',.001,40,1e-12,48,300)]:
  sm.core.logkernel=original
  rows.append(evaluate('baseline_'+label,order,dps))
  if label=='main':rows.append(evaluate('baseline_equilibrium',order,dps,True))
  t=time.time();kernel=FullKernel(step,L,tol);print('tabulation seconds',time.time()-t,flush=True)
  sm.core.logkernel=kernel
  rows.append(evaluate('full_quantum_'+label,order,dps))
  rows.append(evaluate('full_quantum_equilibrium_'+label,order,dps,True))
  v=dict(grid_step_eV=step,sideband_L=L,Fourier_tolerance=tol,off_grid_max_log_error=float(kernel.interp_error),sideband_mass=float(np.sum(np.exp(kernel.lp))),slow_density_mass=float(quad(lambda e:np.exp(kernel.spline(e)),-2,4,epsabs=1e-11,epsrel=1e-11)[0]))
  moments=[quad(lambda e:e**n*np.exp(kernel.spline(e)),-2,4,epsabs=1e-11,epsrel=1e-11)[0] for n in [1,2]]
  v.update(slow_density_mean_eV=moments[0],slow_density_variance_eV2=moments[1]-moments[0]**2,exact_variance_eV2=float(kb.exact_ohmic_var(300,.02)))
  validation[label]=v
 savecsv('summary.csv',rows)
 a={r['case']:r for r in rows}
 validation['rate_ratio']=a['full_quantum_refined']['R_s']/a['baseline_refined']['R_s']
 validation['full_quantum_R_refinement_relative']=abs(a['full_quantum_main']['R_s']/a['full_quantum_refined']['R_s']-1)
 validation['baseline_R_refinement_relative']=abs(a['baseline_main']['R_s']/a['baseline_refined']['R_s']-1)
 validation['budget_equal']=len(set(r['total_coupling_H2_eV2'] for r in rows))==1
 validation['all_positive']=all(r['min_probability']>0 for r in rows)
 validation['max_raw_DB_log_error']=max(r['direct_DB_max_log_error'] for r in rows)
 assert validation['max_raw_DB_log_error']<1e-6
 assert validation['full_quantum_R_refinement_relative']<1e-6
 assert validation['all_positive'] and validation['budget_equal']
 assert all(abs(validation[k]['slow_density_mass']-1)<1e-7 for k in ['main','refined'])
 assert all(r['max_node_continuity_relative']<1e-15 and r['heat_chemical_relative']<1e-15 for r in rows)
 validation['assertions_passed']=True
 (OUT/'validation.json').write_text(json.dumps(validation,indent=2))
 prov={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((ROOT/'inputs').rglob('*')) if p.is_file()}
 (ROOT/'source_sha256.json').write_text(json.dumps(prov,indent=2))
 print(json.dumps(validation,indent=2),flush=True)
if __name__=='__main__':main()
