"""Unipolar FD-Gaussian Poisson/DD SCLC forward solver.
Symmetric ideal DOS-center contacts c=.5, no injection barrier, no series/shunt.
Current uses the same thermodynamic SG discretization as bipolar device.
Not an inferred experimental contact model; report this boundary assumption.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
import numpy as np
from scipy.sparse import coo_matrix,diags
from scipy.sparse.linalg import spsolve
from device import KB,Q,EPS0,sg_flux,log_bernoulli
from fd_eos import GaussianDOS
from transport import mobility_factor,dilute_shape
class SCLC:
 def __init__(self,T=300,sigma=.074,d_cm=170e-7,mu300=1.3e-4,gamma=.001145,mode='density_only',N=161,c_cut=.1,field_cap=235000):
  self.T=T;self.sigma=sigma;self.vt=KB*T;self.d=d_cm;self.N0=2.4e20;self.eps=3.5*EPS0;self.a=self.N0**(-1/3);self.mu=mu300*dilute_shape(sigma,T,mode);self.gamma=gamma;self.mode=mode;self.N=N;self.c_cut=c_cut;self.field_cap=field_cap
  t=np.linspace(-1,1,N);self.x=(1+np.tanh(5*t)/np.tanh(5))/2;self.dx=np.diff(self.x)*self.d;self.vol=(self.dx[:-1]+self.dx[1:])/2;self.eos=GaussianDOS(sigma,T)
  rr=[];cc=[]
  for i in range(N):
   for a in range(2):
    for j in range(max(0,i-1),min(N,i+2)):
     for b in range(2):rr.append(2*i+a);cc.append(2*j+b)
  self.rr=np.array(rr);self.cc=np.array(cc);self.colors=self.cc%6
 def evaluate(self,z,V):
  psi,v=z.T;eta=v-psi;ln=self.eos.logc(eta);c=np.exp(ln)
  # Piecewise contact-referenced unknowns remove large-offset cancellation
  # at BOTH highly conducting majority contacts. The physical fields are
  # psi+offset and v+offset; eta is exactly v-psi.
  offset=-V/self.vt*(self.x>.5);dp=np.diff(psi)+np.diff(offset);dq=np.diff(v)+np.diff(offset)
  F=abs(dp)*self.vt/self.dx;de=np.diff(eta);dl=np.diff(ln)
  g=np.divide(de,dl,out=1/self.eos.dlogc((eta[:-1]+eta[1:])/2),where=abs(de)>1e-7)
  affinity=dq/g;az=abs(affinity);le=np.full_like(az,-np.inf);nz=az>0
  le[nz]=np.maximum(affinity[nz],0)+np.log(-np.expm1(-az[nz]))
  lj=np.log(Q*self.mu*self.N0*self.vt)+np.log(g)-np.log(self.dx)+ln[:-1]+log_bernoulli(dp/g)+le
  J=-np.sign(affinity)*np.exp(lj)
  ce=np.exp(self.eos.logc((eta[:-1]+eta[1:])/2))
  fac=mobility_factor(ce,F,self.sigma,self.T,self.a,self.mode,self.gamma,self.c_cut,self.field_cap)
  return dict(J=J*fac,c=c,c_edge=ce,F=F,mobility=self.mu*fac,g=g,dp=dp)
 def residual(self,z,V):
  o=self.evaluate(z,V);r=np.zeros_like(z);r[1:-1,0]=np.diff(o['dp']/self.dx)+Q*self.N0/(self.eps*self.vt)*self.vol*o['c'][1:-1];r[1:-1,1]=np.diff(o['J']);r[0,:]=z[0,:];r[-1,:]=z[-1,:];return r.ravel()
 def jacobian(self,z,V):
  vals=np.zeros(len(self.rr));h=1e-5
  for col in range(6):
   zp=z.copy().ravel();zm=z.copy().ravel();zp[col::6]+=h;zm[col::6]-=h;der=(self.residual(zp.reshape(z.shape),V)-self.residual(zm.reshape(z.shape),V))/(2*h);m=self.colors==col;vals[m]=der[self.rr[m]]
  return coo_matrix((vals,(self.rr,self.cc)),shape=(z.size,z.size)).tocsc()
 def solve(self,V,z=None):
  if z is None:z=np.zeros((self.N,2))
  else:z=z.copy()
  for it in range(200):
   r=self.residual(z,V);A=self.jacobian(z,V);sc=np.maximum(abs(A).max(axis=1).toarray().ravel(),1e-250);norm=max(abs(r)/sc);o=self.evaluate(z,V);spread=np.ptp(o['J']);j=np.mean(o['J']);gate=spread<max(abs(j)*1e-7,1e-15)
   if norm<1e-9 and gate:return z,dict(success=True,iterations=it,norm=norm,J=float(j),spread=float(spread),density_cut_fraction=float(np.sum(self.dx*(o['c_edge']>self.c_cut))/self.d),field_cut_fraction=float(np.sum(self.dx*(o['F']>self.field_cap))/self.d),mean_c=float(np.dot(self.vol,o['c'][1:-1])/sum(self.vol)),max_F=float(max(o['F'])))
   step=spsolve(diags(1/sc)@A,-r/sc).reshape(z.shape);alpha=min(1.,4/max(np.max(abs(step)),1e-30))
   for ls in range(35):
    try:norm2=max(abs(self.residual(z+alpha*step,V))/sc)
    except (ValueError,FloatingPointError):norm2=np.inf
    if norm2<norm or (norm<1e-9 and norm2<1e-9):z+=alpha*step;break
    alpha*=.5
   else:return z,dict(success=False,norm=norm,reason='line search')
  return z,dict(success=False,norm=norm,reason='iterations')
 def advance(self,z,V0,V1,depth=0):
  guess=z.copy();guess-=((V1-V0)/self.vt)*(self.x-(self.x>.5))[:,None];zz,r=self.solve(V1,guess)
  if r['success']:return zz,r
  if depth>=10:raise RuntimeError(str(r))
  Vm=(V0+V1)/2;zm,_=self.advance(z,V0,Vm,depth+1);return self.advance(zm,Vm,V1,depth+1)
 def sweep(self,Vs):
  z,r=self.solve(0)
  if not r['success']:raise RuntimeError(str(r))
  old=0.;rows=[]
  for V in sorted(Vs):
   z,r=self.advance(z,old,V);old=V;rows.append(dict(V=float(V),**r))
  return rows
if __name__=='__main__':
 import json
 for mode in ['baseline','density_only','egdm']:
  d=SCLC(mode=mode);r=d.sweep([.5,1,2,4]);print(json.dumps(dict(mode=mode,rows=r)),flush=True)
