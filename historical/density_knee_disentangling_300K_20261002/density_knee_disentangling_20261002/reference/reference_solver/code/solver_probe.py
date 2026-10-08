"""Isolated diagnostic variants. All stored original files are unchanged."""
from pathlib import Path
import sys,json,time,traceback,hashlib,warnings
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
ROOT=Path(__file__).resolve().parents[1]
ORIGINAL=ROOT/'reference_audit'
sys.path.insert(0,str(ORIGINAL/'code'))
from rate_audit import ControlledPF,device,audit,PFParams,BASE,channels
ROOT=Path(__file__).resolve().parents[1]
def save(p,o):p.write_text(json.dumps(o,indent=2,allow_nan=False,default=float),encoding='utf-8')

class ProbeDevice(ControlledPF):
 def __init__(self,*a,arithmetic='original',**kw):
  self.arithmetic=arithmetic;super().__init__(*a,**kw)
 def transport(self,y):
  if self.arithmetic=='original':return super().transport(y)
  yy=y.astype(np.clongdouble if np.iscomplexobj(y) else np.longdouble) if self.arithmetic=='extended' else y
  n=np.exp(yy[:,1]);p=np.exp(yy[:,2]);dp=np.diff(yy[:,0]);bp=self.bern(dp);bm=self.bern(-dp)
  # Exact SG identity B(-z)=exp(z)B(z), with differences in log densities.
  jn=n[:-1]*bm*np.expm1(np.diff(yy[:,1])-dp)/self.h
  jp=-p[:-1]*bp*np.expm1(np.diff(yy[:,2])+dp)/self.h*(self.p.mu_p/self.p.mu_n)
  bim=self.p.beta*self.n0**2*(n*p-self.ni_s**2)
  return n,p,jn,jp,bim
 def trap_rates(self,y):
  yy=y.astype(np.clongdouble if np.iscomplexobj(y) else np.longdouble) if self.arithmetic=='extended' else y
  return super().trap_rates(yy)

def diagnostics(d,y,V):
 o,a=audit(d,y,V,0)
 return {k:o[k] for k in ['J_mAcm2','pair_net_mAcm2','bim_net_mAcm2','electron_continuity_error_Acm2','hole_continuity_error_Acm2','current_spread_Acm2']}

def solve_probe(d,V,y,case,linear_scaled=False,tolerance=2e-12,maxiter=100):
 y=y.copy();y[0]=[0,*d.bc[0]];y[-1]=[(d.vbi-V)/d.vt,*d.bc[1]];history=[]
 for it in range(maxiter):
  r=np.asarray(d.residual(y,V,0),float);J=d.jacobian(y,V,0)
  scale=np.maximum(np.asarray(abs(J).max(axis=1).toarray()).ravel(),1e-20);norm=float(max(abs(r)/scale))
  if norm<tolerance:return y,dict(iterations=it,residual_scaled=norm,history=history)
  Js=diags(1/scale)@J if linear_scaled else J;rs=r/scale if linear_scaled else r
  step=spsolve(Js,-rs).reshape(y.shape)
  linear_relative=float(np.linalg.norm((J@step.ravel()+r)/scale)/max(np.linalg.norm(r/scale),1e-300))
  alpha=min(1.,3./max(np.max(abs(step[:,1:])),1e-30),20./max(np.max(abs(step[:,0])),1e-30))
  trial_norms=[]
  for ls in range(25):
   trial=y+alpha*step;rn=d.residual(trial,V,0);nn=float(max(abs(rn)/scale));trial_norms.append([float(alpha),nn])
   if np.all(np.isfinite(rn)) and nn<norm:y=trial;break
   alpha*=.5
  else:
   if norm<2e-8 and tolerance==2e-12:return y,dict(iterations=it,residual_scaled=norm,history=history,loose_line_search_fallback=True)
   np.savez_compressed(ROOT/'data'/f'{case}_failed_newton.npz',y=y,V=V,residual=r,jacobian=J.toarray(),scale=scale,step=step)
   history.append(dict(iteration=it,norm=norm,step_max=float(max(abs(step.ravel()))),linear_relative_residual=linear_relative,trials=trial_norms))
   save(ROOT/'data'/f'{case}_failed_newton.json',dict(V=V,history=history,diagnostics=diagnostics(d,y,V)))
   raise RuntimeError(f'line search failed V={V},norm={norm},it={it}')
  history.append(dict(iteration=it,norm=norm,step_max=float(max(abs(step.ravel()))),linear_relative_residual=linear_relative,trials=trial_norms))
 raise RuntimeError(f'maxiter V={V} norm={norm}')

def run(name,arith,scaled,tol,predictor=False):
 d=ProbeDevice(device().p,81,trap_e_depth=.164,arithmetic=arith)
 phi=d.x*d.vbi/d.vt;y=np.c_[phi,d.bc[0,0]+phi,d.bc[0,1]-phi];trace=[];rows=[];start=time.monotonic()
 try:
  y,info=solve_probe(d,0,y,name,scaled,tol);rows.append(dict(V=0,**diagnostics(d,y,0),**{k:v for k,v in info.items() if k!='history'}));prev=0
  for target in [-5,-10,-15,-20]:
   lastV=prev
   for V in np.linspace(prev,target,int(np.ceil(abs(target-prev)/.1))+1)[1:]:
    if predictor:y[:,0]-=(float(V)-lastV)/d.vt*d.x
    y,info=solve_probe(d,float(V),y,name,scaled,tol);trace.append(dict(V=float(V),**{k:v for k,v in info.items() if k!='history'}));lastV=float(V)
   prev=target;out=dict(V=target,**diagnostics(d,y,target),**{k:v for k,v in info.items() if k!='history'});rows.append(out);np.savez_compressed(ROOT/'data'/f'{name}_{target}.npz',y=y,V=target);print(name,out,flush=True)
  status='completed';error=None
 except Exception as exc:status='failed';error=dict(error=str(exc),traceback=traceback.format_exc());print(name,error['error'],flush=True)
 save(ROOT/'data'/f'{name}.json',dict(name=name,arithmetic=arith,linear_row_scaling=scaled,tolerance=tol,linear_potential_predictor=predictor,rows=rows,trace=trace,status=status,error=error,seconds=time.monotonic()-start))
if __name__=='__main__':
 name=sys.argv[1];cases={'original':('original',False,2e-12),'scaled':('original',True,2e-12),'expm1_scaled':('expm1',True,2e-12),'extended_scaled':('extended',True,2e-12),'tight_scaled':('original',True,2e-15),'tight_expm1':('expm1',True,2e-15),'tight_extended':('extended',True,2e-15),'pred_original':('original',False,2e-12,True),'pred_tight':('expm1',True,2e-15,True),'pred_extended':('extended',True,2e-15,True)}
 run(name,*cases[name])
