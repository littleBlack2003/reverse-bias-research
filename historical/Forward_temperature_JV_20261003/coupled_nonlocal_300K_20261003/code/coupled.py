"""Conditional one-level FD DOS + lattice-gas SG + sequential Marcus device.
All quantities cm,eV,V,s,A. No calibration; fixed-temperature diagnostic only.
"""
import sys,json,numpy as np
from pathlib import Path
from scipy.sparse import csr_matrix,coo_matrix
ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT.parent/'localselfconsistent_srh_300K_20261003'
sys.path.insert(0,str(OLD/'code'))
from integrate import LocalDevice,PFParams,BASE,solve_stable,advance
HBAR=6.582119569e-16
class Coupled(LocalDevice):
 def __init__(self,N,distance_nm=1.,t0=.001,lam=.3,xi_nm=.5):
  super().__init__(PFParams(**dict(BASE,nonlocal_fraction=0.)),N,trap_e_depth=.65)
  if not np.isfinite(distance_nm) or distance_nm<0 or not np.all(np.isfinite([t0,lam,xi_nm])) or min(t0,lam,xi_nm)<=0: raise ValueError('positive finite kernel parameters and nonnegative distance required')
  self.distance=distance_nm*1e-7;self.t0=t0;self.lam=lam;self.xi=xi_nm*1e-7
  dx=self.h*self.p.d;x=self.x*self.p.d;lo=2e-7;hi=self.p.d-2e-7
  vols=np.maximum(0,np.minimum(x+dx/2,hi)-np.maximum(x-dx/2,lo));ids=np.flatnonzero(vols)
  self.xt=np.repeat(x[ids],2);self.xH=self.xt+np.tile([-1,1],len(ids))*self.distance;self.xL=2*self.xt-self.xH
  self.a=np.repeat(self.p.Nt*vols[ids]/2,2)
  def W(xx):
   if np.any(xx<x[1]) or np.any(xx>x[-2]): raise ValueError('reaction endpoint must resolve on interior nodes; refine grid or reduce distance')
   rr=[];cc=[];vv=[]
   for i,v in enumerate(xx):
    k=min(int(np.searchsorted(x,v)),N-1)
    if x[k]==v: rr+=[i];cc+=[k];vv+=[1.]
    else:
     f=(v-x[k-1])/dx;rr += [i,i];cc +=[k-1,k];vv +=[1-f,f]
   return csr_matrix((vv,(rr,cc)),shape=(len(xx),N))
  self.H=W(self.xH);self.T=W(self.xt);self.L=W(self.xL)
  radius=int(np.ceil(2*self.distance/dx))+2;rr=[];cc=[]
  for i in range(N):
   for j in range(max(0,i-radius),min(N,i+radius+1)):
    for a in range(3):
     for b in range(3):rr.append(3*i+a);cc.append(3*j+b)
  self.rr=np.array(rr);self.cc=np.array(cc);self.ncolors=3*(2*radius+1);self.colors=self.cc%self.ncolors
 def kernel(self,de):
  return 2*np.pi/HBAR*self.t0**2*np.exp(-2*self.distance/self.xi-(de+self.lam)**2/(4*self.lam*self.vt))/np.sqrt(4*np.pi*self.lam*self.vt)
 def evaluate(self,z,V):
  y=self.y_from_z(z,V);an=np.exp(y[:,1])*self.n0/self.p.Nc;ap=np.exp(y[:,2])*self.n0/self.p.Nc
  ns=self.p.Nc/self.n0*an/(1+an);ps=self.p.Nc/self.n0*ap/(1+ap);dp=np.diff(z[:,0])
  jn=np.exp(y[:-1,1])*self.bern(-dp)*np.expm1(np.diff(z[:,1]))/self.h*self.j0/((1+an[:-1])*(1+an[1:]))
  jp=-np.exp(y[:-1,2])*self.bern(dp)*np.expm1(np.diff(z[:,2]))/self.h*self.j0*self.p.mu_p/self.p.mu_n/((1+ap[:-1])*(1+ap[1:]))
  phi=z[:,0]*self.vt;muH=self.p.contact_barrier-self.vt*z[:,2];muL=self.p.contact_barrier+V+self.vt*z[:,1]
  ph=self.H@phi;pt=self.T@phi;pl=self.L@phi;mh=self.H@muH;ml=self.L@muL
  eh=-ph;et=self.p.Eg/2-pt;el=self.p.Eg-pl
  ah=np.exp((eh-mh)/self.vt);al=np.exp((ml-el)/self.vt)
  FH=1/(1+ah);VH=ah/(1+ah);FL=al/(1+al);VL=1/(1+al)
  aH=self.kernel(et-eh)*FH;bH=self.kernel(eh-et)*VH;aL=self.kernel(et-el)*FL;bL=self.kernel(el-et)*VL
  D=aH+bH+aL+bL;f=(aH+aL)/D
  # Exact LDB cycle ratio, stable at common chemical potential.
  affinity=-V/self.vt-self.H@z[:,2]-self.L@z[:,1];j=-aH*bL/D*np.expm1(-affinity)
  dx=self.h*self.p.d;Sn=self.L.T@(self.a*j)/dx;Sp=self.H.T@(self.a*j)/dx
  rhot=self.T.T@(self.a*(.5-f))/dx;feff=.5-rhot/self.p.Nt
  pair=(self.T.T@(self.a*j))/dx/self.p.Nt
  bim=self.p.beta*self.ni**2*np.expm1(V/self.vt+z[:,1]+z[:,2])/((1+an)*(1+ap))
  return dict(y=y,ns=ns,ps=ps,Jn=jn,Jp=jp,Rn=-Sn,Rp=-Sp,f=feff,fb=1-feff,bim=bim,pair=pair,actual_f=f,actual_fb=(bH+bL)/D,j=j,muH=muH,muL=muL,affinity=affinity,eh=eh,et=et,el=el,mh=mh,ml=ml,rhot=rhot,aH=aH,bH=bH,aL=aL,bL=bL)
 def ledger(self,z,V):
  o=self.evaluate(z,V);p=self.physical(z,V);dx=self.h*self.p.d;phi=z[:,0]*self.vt
  E=-np.diff(phi)/dx;charge=self.q*dx*np.sum((o['ps']*self.n0-o['ns']*self.n0+o['rhot'])[1:-1]);fieldcharge=self.eps*(E[-1]-E[0]);gross=self.q*dx*np.sum((o['ps']*self.n0+o['ns']*self.n0+abs(o['rhot']))[1:-1])
  transport=np.sum(o['Jn']*np.diff(o['muL'])+o['Jp']*np.diff(o['muH']))
  reaction=self.q*np.sum(self.a*o['j']*self.vt*o['affinity'])
  bim=self.q*dx*np.sum((o['bim']*(o['muL']-o['muH']))[1:-1])
  terminal=V*p['J_Acm2'];total=transport+reaction+bim
  bath=self.q*np.sum(self.a*o['j']*(o['eh']-o['el']))
  reservoir=self.q*np.sum(self.a*o['j']*((o['mh']-o['eh'])+(o['el']-o['ml'])))
  p.update(V=V,nodes=self.nodes,distance_nm=self.distance/1e-7,charge_relative_error=float(abs(charge-fieldcharge)/max(gross,1e-300)),terminal_Wcm2=float(terminal),transport_dissipation_Wcm2=float(transport),reaction_dissipation_Wcm2=float(reaction),bim_dissipation_Wcm2=float(bim),total_dissipation_Wcm2=float(total),energy_relative_error=float(abs(total-terminal)/max(abs(total),abs(terminal),1e-300)),reaction_bath_heat_Wcm2=float(bath),reaction_reservoir_heat_Wcm2=float(reservoir),reaction_heat_identity_error_Wcm2=float(bath+reservoir-reaction),occupancy_min=float(min(o['actual_f'])),occupancy_max=float(max(o['actual_f'])),max_LDB_cycle_error=float(max(abs(np.log(o['aH'])+np.log(o['bL'])-np.log(o['aL'])-np.log(o['bH'])-o['affinity']))),trap_areal_cm2=float(sum(self.a)),max_density_fraction=float(max(max(o['ns']),max(o['ps']))*self.n0/self.p.Nc))
  return p,o

def run(N,distance):
 d=Coupled(N,distance);z=np.load(OLD/'data'/f'baseline_N{N}_V0.npz')['z'];z,info=solve_stable(d,0,z);assert not info.get('status'),info
 rows=[];trace=[];old=0.
 for target in [0.,-1.,-5.,-15.]:
  for V in np.linspace(old,target,max(1,int(np.ceil(abs(target-old)/.25)))+1)[1:]:z,info=advance(d,z,old,float(V),trace);old=float(V)
  r,o=d.ledger(z,target);r.update(iterations=info['iterations'],scaled_residual=info['scaled_residual']);rows.append(r)
  np.savez_compressed(ROOT/'data'/f'd{distance:g}_N{N}_V{target:g}.npz',z=z,**o)
  print(json.dumps(r),flush=True)
 (ROOT/'data'/f'd{distance:g}_N{N}.json').write_text(json.dumps(dict(rows=rows,trace=trace),indent=2))
 return rows
if __name__=='__main__':run(int(sys.argv[1]),float(sys.argv[2]))
