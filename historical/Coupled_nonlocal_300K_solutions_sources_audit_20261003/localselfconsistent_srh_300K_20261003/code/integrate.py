"""Local four-rate adapter in recovered self-consistent SG/QF device solver.
No material fit. All units inherited: cm, s, V, A, eV. Original files untouched.
"""
from run_density import *
class LocalDevice(DensityDevice):
 def __init__(self,*a,single_pf=False,**kw):
  self.single_pf=single_pf
  super().__init__(*a,**kw)
  assert self.p.nonlocal_fraction==0 and self.p.enhancement==0
 def local(self,z,V):
  y=self.y_from_z(z,V);n=np.exp(y[:,1])*self.n0;p=np.exp(y[:,2])*self.n0
  field=-np.gradient(z[:,0]*self.vt,self.h*self.p.d)
  # Analytic smoothing inherited from baseline; permits complex-step derivatives.
  delta=np.sqrt(self.q*np.sqrt(field*field+1.)/(np.pi*self.eps))
  gn=np.exp(delta/self.vt) if self.single_pf else np.ones_like(field)
  cn=self.p.capture*gn;cp=self.p.capture
  n1=self.ni;p1=self.ni;en=cn*n1;ep=cp*p1
  A=cn*n+ep;B=en+cp*p;D=A+B;f=A/D;fb=B/D
  # Affinity is determined by local carrier quasi-Fermi variables, no new mu ports.
  affinity=V/self.vt+z[:,1]+z[:,2]
  U=self.p.Nt*en*ep*np.expm1(affinity)/D
  return dict(n=n,p=p,f=f,fb=fb,U=U,cn=cn,cp=cp,en=en,ep=ep,D=D,gn=gn,field=field)
 def evaluate(self,z,V):
  a=self.local(z,V);y=self.y_from_z(z,V);ns=a['n']/self.n0;ps=a['p']/self.n0;dp=np.diff(z[:,0])
  jn=ns[:-1]*self.bern(-dp)*np.expm1(np.diff(z[:,1]))/self.h*self.j0
  jp=-ps[:-1]*self.bern(dp)*np.expm1(np.diff(z[:,2]))/self.h*self.j0*self.p.mu_p/self.p.mu_n
  U=a['U'].copy();U[[0,-1]]=0
  bim=self.p.beta*self.ni**2*np.expm1(V/self.vt+z[:,1]+z[:,2])
  return dict(y=y,ns=ns,ps=ps,Jn=jn,Jp=jp,Rn=U,Rp=U.copy(),f=a['f'],fb=a['fb'],bim=bim,pair=-U/self.p.Nt)

def make(N,pf=False):return LocalDevice(PFParams(**dict(BASE,nonlocal_fraction=0.)),N,trap_e_depth=.65,single_pf=pf)
def solve_case(N,pf=False):
 d=make(N,pf);z=np.c_[d.x*d.vbi/d.vt,np.zeros(N),np.zeros(N)];z,info=solve_stable(d,0,z)
 assert not info.get('status'),info
 rows=[];trace=[];old=0.
 for target in [0.,-1.,-5.,-15.]:
  for v in np.linspace(old,target,max(1,int(np.ceil(abs(target-old)/.25)))+1)[1:]:z,info=advance(d,z,old,float(v),trace);old=float(v)
  row,a=details(d,z,target,info);rows.append(row)
  np.savez_compressed(ROOT/'data'/f'{"singlePF" if pf else "baseline"}_N{N}_V{target:g}.npz',z=z,**a)
  print('POINT',N,pf,target,row['J_Acm2'],row['relative_error'],flush=True)
 save(ROOT/'data'/f'{"singlePF" if pf else "baseline"}_N{N}.json',dict(parameters=asdict(d.p),single_pf=pf,rows=rows,trace=trace))
 return rows
if __name__=='__main__':
 mode=sys.argv[1] if len(sys.argv)>1 else 'baseline'
 if mode not in ['baseline','singlePF']:raise ValueError('unknown mode')
 for N in [81,161,321]:solve_case(N,mode=='singlePF')
