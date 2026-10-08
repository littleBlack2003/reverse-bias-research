"""Algebraically equivalent cycle-source evaluation in QF coordinates.
Uses Kn*Kp=exp(V/Vt), so equilibrium is not a subtraction of gross rates.
No changed barriers, attempt rates, band energies or physical pathways.
"""
from temperature_probe import *
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
class StableCycleQF(QFDevice):
 def channel_arrays(self,z,V):
  p=self.p;N=self.nodes;psi=z[:,0];phi=psi*self.vt;y=self.y_from_z(z,V);dtype=z.dtype
  ends=np.c_[np.arange(N),np.clip(np.arange(N)-self.jump,0,N-1),np.clip(np.arange(N)+self.jump,0,N-1)]
  En=np.zeros((N,3),dtype=dtype);Ep=np.zeros_like(En);dn=self.trap_e_depth;dp=p.Eg-dn;c0=p.capture*(1-p.nonlocal_fraction)
  En[:,0]=c0*p.Nc*np.exp(-dn/self.vt);Ep[:,0]=c0*p.Nc*np.exp(-dp/self.vt)
  field=-np.gradient(phi,self.h)/p.d;delta=p.pf_strength*np.sqrt(self.q*np.sqrt(field*field+1.)/(np.pi*self.eps));pref=p.capture*p.Nc*p.nonlocal_fraction/2;w=p.barrier_smoothing_eV
  for ci,(t,end) in enumerate(self.pairs,1):
   du=phi[end]-phi[t]
   for E,depth,energy,mult in [(En,dn,dn-du,self.branch_n),(Ep,dp,dp+du,self.branch_p)]:
    base=self.smoothmax(0.,depth-mult*delta,w);barrier=self.smoothmax(base[t],energy,w);E[t,ci]=pref*np.exp(-barrier/self.vt)
  L=V/self.vt;Kn=-p.contact_barrier/self.vt-(self.vbi-V)/self.vt+psi+dn/self.vt;Kp=L-Kn
  un=z[ends,1];vp=z[ends,2]
  An=En*np.exp(Kn[:,None]+un);Ap=Ep*np.exp(Kp[:,None]+vp)
  total=An.sum(1)+Ap.sum(1)+En.sum(1)+Ep.sum(1);f=(An.sum(1)+Ep.sum(1))/total;fb=(Ap.sum(1)+En.sum(1))/total
  return ends,En,Ep,An,Ap,total,f,fb,Kn,Kp,un,vp
 def evaluate(self,z,V):
  y=self.y_from_z(z,V);ns=np.exp(y[:,1]);ps=np.exp(y[:,2]);dpsi=np.diff(z[:,0])
  jn=ns[:-1]*self.bern(-dpsi)*np.expm1(np.diff(z[:,1]))/self.h*self.j0
  jp=-ps[:-1]*self.bern(dpsi)*np.expm1(np.diff(z[:,2]))/self.h*self.j0*self.p.mu_p/self.p.mu_n
  ends,En,Ep,An,Ap,S,f,fb,Kn,Kp,un,vp=self.channel_arrays(z,V)
  rn=np.zeros(self.nodes,dtype=z.dtype);rp=np.zeros_like(rn);pair=np.zeros_like(rn);mask=np.ones(self.nodes);mask[[0,-1]]=0
  # Every pair event contributes equally to total n and p source.
  for j in range(3):
   for k in range(3):
    G=En[:,j]*Ep[:,k]/S*(-np.expm1(V/self.vt+un[:,j]+vp[:,k]))*mask
    pair+=G
    np.add.at(rn,ends[:,j],-self.p.Nt*G);np.add.at(rp,ends[:,k],-self.p.Nt*G)
  # Unordered same-carrier cycles: net transported particles, no pair source.
  for i in range(3):
   for j in range(i+1,3):
    tn=En[:,j]*An[:,i]/S*(-np.expm1(un[:,j]-un[:,i]))*mask
    tp=Ep[:,j]*Ap[:,i]/S*(-np.expm1(vp[:,j]-vp[:,i]))*mask
    np.add.at(rn,ends[:,i],self.p.Nt*tn);np.add.at(rn,ends[:,j],-self.p.Nt*tn)
    np.add.at(rp,ends[:,i],self.p.Nt*tp);np.add.at(rp,ends[:,j],-self.p.Nt*tp)
  bim=self.p.beta*self.ni**2*np.expm1(V/self.vt+z[:,1]+z[:,2])
  return dict(y=y,ns=ns,ps=ps,Jn=jn,Jp=jp,Rn=rn,Rp=rp,f=f,fb=fb,bim=bim,pair=pair,ends=ends,En_channels=En,Ep_channels=Ep)
 def physical(self,z,V):
  o=self.evaluate(z,V);dx=self.h*self.p.d;Jn=o['Jn'];Jp=o['Jp']
  pairA=self.q*self.p.Nt*dx*np.sum(o['pair'][1:-1]);bimA=self.q*dx*np.sum(o['bim'][1:-1]);expected=pairA-bimA
  nout=Jn[0]-Jn[-1];pout=Jp[-1]-Jp[0];jnl=np.r_[0.,-self.q*dx*np.cumsum((o['Rn']-o['Rp'])[1:-1])];jt=Jn+Jp+jnl
  error=max(abs(nout-expected),abs(pout-expected),np.ptp(jt),abs(jnl[-1]))
  scale=max(abs(pairA),abs(bimA),abs(Jn[0]),abs(Jn[-1]),abs(Jp[0]),abs(Jp[-1]),1e-300)
  budget=(1e-17 if V==0 else 1e-280)+1e-8*scale
  return dict(J_Acm2=float((Jn[0]+Jp[0]+Jn[-1]+Jp[-1])/2),pair_source_Acm2=float(pairA),expected_particle_output_Acm2=float(expected),electron_out_Acm2=float(nout),hole_out_Acm2=float(pout),electron_cycle_error_Acm2=float(nout-expected),hole_cycle_error_Acm2=float(pout-expected),total_current_spread_Acm2=float(np.ptp(jt)),nonlocal_boundary_error_Acm2=float(abs(jnl[-1])),budget_Acm2=float(budget),error_Acm2=float(error),gate_passed=bool(error<=budget),relative_error=float(error/scale))

def solve_stable(d,V,z):
 z=z.copy();z[0]=[0,-V/d.vt,0];z[-1]=[(d.vbi-V)/d.vt,0,-V/d.vt];hist=[]
 for it in range(100):
  r=d.residual(z,V);J=d.jacobian(z,V,0);scale=np.maximum(np.asarray(abs(J).max(axis=1).toarray()).ravel(),np.finfo(float).tiny)
  norm=float(max(abs(r)/scale));phys=d.physical(z,V);hist.append(dict(it=it,scaled_residual=norm,**phys))
  if norm<2e-12 and phys['gate_passed']:return z,dict(iterations=it,scaled_residual=norm,history=hist,**phys)
  step=spsolve(diags(1/scale)@J,-r/scale).reshape(z.shape)
  if not np.all(np.isfinite(step)):return z,dict(status='nonfinite_step',iterations=it,scaled_residual=norm,history=hist,**phys)
  alpha=min(1.,3/max(np.max(abs(step[:,1:])),1e-30),20/max(np.max(abs(step[:,0])),1e-30))
  for ls in range(30):
   zz=z+alpha*step;rr=d.residual(zz,V);nn=float(max(abs(rr)/scale));pp=d.physical(zz,V)
   if np.all(np.isfinite(rr)) and (nn<norm or (norm<2e-12 and nn<2e-12 and pp['relative_error']<phys['relative_error'])):z=zz;break
   alpha*=.5
  else:return z,dict(status='line_search_unresolved',iterations=it,scaled_residual=norm,history=hist,**phys)
 return z,dict(status='maxiter',iterations=it,scaled_residual=norm,history=hist,**phys)

def refine_existing(T,depth,nodes=81):
 name=f'T{T:g}_D{depth:g}_N{nodes}';p=dict(BASE);p['temperature']=T;d=StableCycleQF(PFParams(**p),nodes,trap_e_depth=depth);rows=[]
 for V in [0,-5,-20]:
  z=np.load(ROOT/'data'/f'{name}_{V}.npz')['z'];before=d.physical(z,V);z,info=solve_stable(d,V,z);row=dict(T=T,depth=depth,nodes=nodes,V=V,before=before,**{k:v for k,v in info.items() if k!='history'});rows.append(row)
  np.savez_compressed(ROOT/'data'/f'stable_{name}_{V}.npz',z=z,y=d.y_from_z(z,V),V=V);save(ROOT/'data'/f'stable_{name}_{V}.json',info);print(name,V,row,flush=True)
 save(ROOT/'data'/f'stable_{name}.json',dict(rows=rows))
if __name__=='__main__':refine_existing(float(sys.argv[1]),float(sys.argv[2]),int(sys.argv[3]) if len(sys.argv)>3 else 81)
