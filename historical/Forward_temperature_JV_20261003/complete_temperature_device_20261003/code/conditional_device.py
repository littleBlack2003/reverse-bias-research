"""Explicitly conditional Gaussian-DOS extension, original code preserved.
All state variables remain contact-referenced dimensionless quasi-Fermi potentials.
Gaussian local-equilibrium baths replace single-energy endpoint occupations.
Transport mobility and reaction Marcus parameters are distinct physical objects.
"""
from pathlib import Path
import sys,numpy as np
from scipy.special import expit
from scipy.sparse import coo_matrix
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent/'coupled_nonlocal_300K_20261003/code'))
from coupled import Coupled,LocalDevice,PFParams,BASE,solve_stable,advance
from gaussian_transport import EGDM,gaussian_eos,inverse_activity_flux,normal_rule,KB

class GaussianDevice(Coupled):
 def __init__(self,N,T=300.,mode='gaussian',sigma_n=.1,sigma_p=.1,a_n=1e-7,a_p=1e-7,order=128):
  if a_n!=a_p:raise ValueError('this first integration assumes equal Ns; closure library supports independent a')
  super().__init__(N)
  rr,cc,ncolors,colors=self.rr,self.cc,self.ncolors,self.colors
  # In the strict lattice model Ns=a^-3; density changes must be explicit.
  Nc=a_n**-3 if mode.startswith('gaussian') else BASE['Nc']
  LocalDevice.__init__(self,PFParams(**dict(BASE,temperature=T,Nc=Nc,nonlocal_fraction=0.)),N,trap_e_depth=.65)
  self.rr,self.cc,self.ncolors,self.colors=rr,cc,ncolors,colors
  self.mode=mode;self.order=order;self.sigma_n=sigma_n;self.sigma_p=sigma_p
  # Match only the ZERO-DENSITY, ZERO-FIELD mobility at 300 K, not all states.
  self.mn=EGDM(sigma_n,a_n,BASE['mu_n']*np.exp(.42*(sigma_n/(KB*300))**2))
  self.mp=EGDM(sigma_p,a_p,BASE['mu_p']*np.exp(.42*(sigma_p/(KB*300))**2))
 def jacobian(self,z,V,light):
  data=np.zeros(len(self.rr));h=2e-5
  for color in range(self.ncolors):
   zp=z.copy().ravel();zm=zp.copy();zp[color::self.ncolors]+=h;zm[color::self.ncolors]-=h
   deriv=(self.residual(zp.reshape(z.shape),V)-self.residual(zm.reshape(z.shape),V))/(2*h)
   mask=self.colors==color;data[mask]=deriv[self.rr[mask]]
  return coo_matrix((data,(self.rr,self.cc)),shape=(z.size,z.size)).tocsc()
 def evaluate(self,z,V):
  if not self.mode.startswith('gaussian'):
   o=super().evaluate(z,V)
   if self.mode=='mobility_only':
    # Diagnostic only: old one-level EOS, rates, contacts and finite-occupancy flux.
    # T-only factor is exactly 1 at 300 K; no claim this is complete EGDM.
    rn=np.exp(-.42*(self.sigma_n/KB)**2*(1/self.p.temperature**2-1/300**2))
    rp=np.exp(-.42*(self.sigma_p/KB)**2*(1/self.p.temperature**2-1/300**2))
    o['Jn']*=rn;o['Jp']*=rp
   elif self.mode!='legacy':raise ValueError('unknown mode')
   return o
  T=self.p.temperature;vt=self.vt;Nc=self.p.Nc;dx=self.h*self.p.d
  y=self.y_from_z(z,V);etan=y[:,1]+np.log(self.n0/Nc);etap=y[:,2]+np.log(self.n0/Nc)
  cn,dcn,g3n=gaussian_eos(etan,self.sigma_n/vt,self.order);cp,dcp,g3p=gaussian_eos(etap,self.sigma_p/vt,self.order)
  nm=Nc*gaussian_eos((etan[:-1]+etan[1:])/2,self.sigma_n/vt,self.order)[0]
  pm=Nc*gaussian_eos((etap[:-1]+etap[1:])/2,self.sigma_p/vt,self.order)[0]
  F=abs(np.diff(z[:,0]))*vt/dx
  Tmu=300. if self.mode=='gaussian_mu300' else T
  mn=self.mn.evaluate(Tmu,nm,F);mp=self.mp.evaluate(Tmu,pm,F)
  if self.mode=='gaussian_constmu':
   mn['mu']=np.full_like(F,BASE['mu_n']);mp['mu']=np.full_like(F,BASE['mu_p'])
  jn=inverse_activity_flux(etan,z[:,0],T,Nc,self.sigma_n,dx,mn['mu'],'n',np.diff(z[:,1]),order=self.order)
  jp=inverse_activity_flux(etap,z[:,0],T,Nc,self.sigma_p,dx,mp['mu'],'p',np.diff(z[:,2]),order=self.order)
  phi=z[:,0]*vt;muH=self.p.contact_barrier-vt*z[:,2];muL=self.p.contact_barrier+V+vt*z[:,1]
  ph=self.H@phi;pt=self.T@phi;pl=self.L@phi;mh=self.H@muH;ml=self.L@muL
  eh=-ph;et=self.p.Eg/2-pt;el=self.p.Eg-pl
  nodes,w=normal_rule(self.order)
  eH=eh[:,None]+self.sigma_p*nodes;eL=el[:,None]+self.sigma_n*nodes
  fH=expit((mh[:,None]-eH)/vt);bHocc=expit((eH-mh[:,None])/vt)
  fL=expit((ml[:,None]-eL)/vt);bLocc=expit((eL-ml[:,None])/vt)
  aH=(self.kernel(et[:,None]-eH)*fH)@w;bH=(self.kernel(eH-et[:,None])*bHocc)@w
  aL=(self.kernel(et[:,None]-eL)*fL)@w;bL=(self.kernel(eL-et[:,None])*bLocc)@w
  D=aH+bH+aL+bL;f=(aH+aL)/D
  affinity=-V/vt-self.H@z[:,2]-self.L@z[:,1]
  j=-aH*bL/D*np.expm1(-affinity)
  Sn=self.L.T@(self.a*j)/dx;Sp=self.H.T@(self.a*j)/dx
  rhot=self.T.T@(self.a*(.5-f))/dx;feff=.5-rhot/self.p.Nt
  pair=(self.T.T@(self.a*j))/dx/self.p.Nt
  A=V/vt+z[:,1]+z[:,2]
  # Reversible phenomenological bimolecular law in chemical-potential form.
  # beta remains constant to isolate closure changes; not claimed Langevin calibration.
  bim=self.p.beta*Nc*Nc*np.exp(np.log(cn)+np.log(cp)-A)*np.expm1(A)
  return dict(y=y,ns=Nc*cn/self.n0,ps=Nc*cp/self.n0,Jn=jn,Jp=jp,Rn=-Sn,Rp=-Sp,f=feff,fb=1-feff,bim=bim,pair=pair,actual_f=f,actual_fb=(bH+bL)/D,j=j,muH=muH,muL=muL,affinity=affinity,eh=eh,et=et,el=el,mh=mh,ml=ml,rhot=rhot,aH=aH,bH=bH,aL=aL,bL=bL,g3n=g3n,g3p=g3p,mu_n_face=mn['mu'],mu_p_face=mp['mu'],reduced_field_n=mn['reduced_field'],reduced_field_p=mp['reduced_field'],transport_window_n=mn['comparison_window'],transport_window_p=mp['comparison_window'],density_n=mn['occupation'],density_p=mp['occupation'])
 def ledger(self,z,V):
  p,o=super().ledger(z,V)
  if self.mode.startswith('gaussian'):
   # Center-energy heat decomposition is not valid for DOS-integrated transitions.
   for k in ['reaction_bath_heat_Wcm2','reaction_reservoir_heat_Wcm2','reaction_heat_identity_error_Wcm2']:p.pop(k,None)
   p.update(mu_n_min=float(min(o['mu_n_face'])),mu_n_max=float(max(o['mu_n_face'])),mu_p_min=float(min(o['mu_p_face'])),mu_p_max=float(max(o['mu_p_face'])),max_reduced_field=float(max(o['reduced_field_n'].max(),o['reduced_field_p'].max())),fraction_faces_high_density=float(np.mean((o['density_n']>1e-2)|(o['density_p']>1e-2))),fraction_faces_low_density=float(np.mean((o['density_n']<1e-6)|(o['density_p']<1e-6))),fraction_faces_outside_comparison_window=float(np.mean(~(o['transport_window_n']&o['transport_window_p']))),max_density_fraction_n=float(max(o['density_n'])),max_density_fraction_p=float(max(o['density_p'])))
  p.update(T=self.p.temperature,mode=self.mode)
  return p,o
