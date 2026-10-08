"""Equivalent-contact-referenced quasi-Fermi coordinates; no rate/physics changes.
The physical budget supplements, rather than replaces, the old scaled test.
"""
from solver_probe import *
from dataclasses import asdict
class QFDevice(ControlledPF):
 def y_from_z(self,z,V):
  psi=z[:,0];psiR=(self.vbi-V)/self.vt
  return np.c_[psi,self.bc[1,0]+psi-psiR+z[:,1],self.bc[0,1]-psi+z[:,2]]
 def z_from_y(self,y,V):
  psiR=(self.vbi-V)/self.vt
  return np.c_[y[:,0],y[:,1]-self.bc[1,0]-y[:,0]+psiR,y[:,2]-self.bc[0,1]+y[:,0]]
 def evaluate(self,z,V):
  y=self.y_from_z(z,V);ns=np.exp(y[:,1]);ps=np.exp(y[:,2]);dp=np.diff(z[:,0])
  jn=ns[:-1]*self.bern(-dp)*np.expm1(np.diff(z[:,1]))/self.h*self.j0
  jp=-ps[:-1]*self.bern(dp)*np.expm1(np.diff(z[:,2]))/self.h*self.j0*self.p.mu_p/self.p.mu_n
  rn,rp,f,E,delta,links,fb=self.trap_rates(y)
  bim=self.p.beta*self.n0**2*(ns*ps-self.ni_s**2)
  return dict(y=y,ns=ns,ps=ps,Jn=jn,Jp=jp,Rn=rn,Rp=rp,f=f,fb=fb,bim=bim)
 def residual(self,z,V,light=0):
  o=self.evaluate(z,V);r=np.zeros_like(z);rho=o['ps']-o['ns']+self.p.Nt/self.n0*(.5-o['f']);G=light*self.p.Jgen/(self.q*self.p.d)
  r[1:-1,0]=np.diff(z[:,0],n=2)/self.h**2+rho[1:-1]
  r[1:-1,1]=np.diff(o['Jn']/self.j0)/self.h-(o['Rn']+o['bim']-G)[1:-1]/self.r0
  r[1:-1,2]=np.diff(o['Jp']/self.j0)/self.h+(o['Rp']+o['bim']-G)[1:-1]/self.r0
  r[0]=z[0]-np.array([0.,-V/self.vt,0.]);r[-1]=z[-1]-np.array([(self.vbi-V)/self.vt,0.,-V/self.vt]);return r.ravel()
 def physical(self,z,V):
  o=self.evaluate(z,V);dx=self.h*self.p.d
  nsrc=-self.q*dx*sum((o['Rn']+o['bim'])[1:-1]);psrc=-self.q*dx*sum((o['Rp']+o['bim'])[1:-1]);Jn=o['Jn'];Jp=o['Jp']
  nout=Jn[0]-Jn[-1];pout=Jp[-1]-Jp[0]
  jnl=np.r_[0.,-self.q*dx*np.cumsum((o['Rn']-o['Rp'])[1:-1])];jt=Jn+Jp+jnl
  gross_scale=max(abs(nsrc),abs(psrc),1e-300);budget=1e-17+1e-5*gross_scale
  error=max(abs(nout-nsrc),abs(pout-psrc),np.ptp(jt))
  return dict(J_Acm2=float((Jn[0]+Jp[0]+Jn[-1]+Jp[-1])/2),electron_source_Acm2=float(nsrc),hole_source_Acm2=float(psrc),electron_out_Acm2=float(nout),hole_out_Acm2=float(pout),electron_error_Acm2=float(nout-nsrc),hole_error_Acm2=float(pout-psrc),total_current_spread_Acm2=float(np.ptp(jt)),budget_Acm2=float(budget),error_Acm2=float(error),gate_passed=bool(error<=budget))

def qsolve(d,V,z,gate=True):
 z=z.copy();z[0]=[0.,-V/d.vt,0.];z[-1]=[(d.vbi-V)/d.vt,0.,-V/d.vt];hist=[]
 for it in range(100):
  r=d.residual(z,V);J=d.jacobian(z,V,0);scale=np.maximum(np.asarray(abs(J).max(axis=1).toarray()).ravel(),1e-20)
  norm=float(max(abs(r)/scale));phys=d.physical(z,V)
  hist.append(dict(it=it,scaled_residual=norm,**phys))
  if norm<2e-12 and (phys['gate_passed'] or not gate):return z,dict(iterations=it,scaled_residual=norm,history=hist,**phys)
  step=spsolve(diags(1/scale)@J,-r/scale).reshape(z.shape);alpha=min(1.,3/max(np.max(abs(step[:,1:])),1e-30),20/max(np.max(abs(step[:,0])),1e-30))
  for ls in range(30):
   zz=z+alpha*step;rr=d.residual(zz,V);nn=float(max(abs(rr)/scale))
   # Near floating-point Poisson floor, accept an improving physical budget
   # only while remaining below the same strict scaled-residual threshold.
   pp=d.physical(zz,V)
   if np.all(np.isfinite(rr)) and (nn<norm or (norm<2e-12 and nn<2e-12 and pp['error_Acm2']<phys['error_Acm2'])):z=zz;break
   alpha*=.5
  else:return z,dict(iterations=it,scaled_residual=norm,history=hist,status='line_search_unresolved',**phys)
 return z,dict(iterations=it,scaled_residual=norm,history=hist,status='maxiter',**phys)

def run_qf(name,depth=.164,nodes=81,gate=True,dv=.1):
 d=QFDevice(device(nodes).p,nodes,trap_e_depth=depth);z=np.c_[d.x*d.vbi/d.vt,np.zeros(nodes),np.zeros(nodes)];trace=[];rows=[];prev=0.
 z,info=qsolve(d,0,z,gate);rows.append(dict(V=0,**{k:v for k,v in info.items() if k!='history'}));saved=[]
 for target in [-5,-10,-15,-20]:
  lastV=prev
  for V in np.linspace(prev,target,int(np.ceil(abs(target-prev)/dv))+1)[1:]:
   shift=-(float(V)-lastV)/d.vt
   z[:,0]+=shift*d.x;z[:,1]+=shift*(1-d.x);z[:,2]+=shift*d.x
   z,info=qsolve(d,float(V),z,gate);trace.append(dict(V=float(V),**{k:v for k,v in info.items() if k!='history'}));lastV=float(V)
   if info.get('status'):
    np.savez_compressed(ROOT/'data'/f'{name}_unresolved_{V:.4f}.npz',z=z,y=d.y_from_z(z,V),V=V);save(ROOT/'data'/f'{name}_unresolved_{V:.4f}.json',info)
    print(name,'UNRESOLVED',V,info['scaled_residual'],info['error_Acm2'],info['budget_Acm2'],flush=True)
    # Continue permitted mathematical guesses, but all unresolved gates retained.
   if info.get('status')=='maxiter':raise RuntimeError(f'Maxiter at actual voltage {V}; no later target may be labeled solved')
  prev=target;out=dict(V=target,**{k:v for k,v in info.items() if k!='history'});rows.append(out)
  np.savez_compressed(ROOT/'data'/f'{name}_{target}.npz',z=z,y=d.y_from_z(z,target),V=target);print(name,out,flush=True)
  save(ROOT/'data'/f'{name}.json',dict(name=name,depth=depth,nodes=nodes,gate_enabled=gate,voltage_step_V=dv,parameters=asdict(d.p),rows=rows,trace=trace))
if __name__=='__main__':
 from dataclasses import asdict
 name=sys.argv[1]
 if name=='qf_asym':run_qf(name)
 elif name=='qf_no_gate':run_qf(name,gate=False)
 elif name=='qf_symmetric':run_qf(name,depth=.65)
 elif name=='qf_asym_161':run_qf(name,nodes=161)
