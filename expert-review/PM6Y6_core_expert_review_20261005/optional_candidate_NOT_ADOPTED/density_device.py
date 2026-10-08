"""Frozen-parameter bipolar DD with an audited warm finite-density candidate.
No beta/Eg/G/contact/Rsh tuning. Low-T extrapolation opt-in is deliberately absent.
The baseline FD generalized diffusion is retained exactly once.
"""
import numpy as np
from device import Device
from extend_device import parameters,ANCHOR
from field_device import FieldDevice
from transport import mobility_factor,dilute_shape
class DensityDevice(FieldDevice):
 def __init__(self,p,N=321,mode='density_only',amplitude_n=1.,amplitude_p=1.,c_cut=.1,field_cap=None):
  _,gn,gp=parameters(p.T)
  self.mode=mode;self.amplitude_n=amplitude_n;self.amplitude_p=amplitude_p;self.c_cut=c_cut
  super().__init__(p,N,gamma_n=gn,gamma_p=gp,field_cap=ANCHOR['field_cap'] if field_cap is None else field_cap)
  self.mn=p.mu_n_300*dilute_shape(p.sigma_n_eV,p.T,mode)*amplitude_n
  self.mp=p.mu_p_300*dilute_shape(p.sigma_p_eV,p.T,mode)*amplitude_p
 def evaluate(self,z,V):
  o=Device.evaluate(self,z,V);F=abs(np.diff(z[:,0]))*self.vt/self.dx
  cn=np.exp(self.en.logc((o['eta_n'][:-1]+o['eta_n'][1:])/2));cp=np.exp(self.ep.logc((o['eta_p'][:-1]+o['eta_p'][1:])/2))
  a=self.p.N0_cm3**(-1/3)
  fn=mobility_factor(cn,F,self.p.sigma_n_eV,self.p.T,a,self.mode,self.gamma_n,self.c_cut,self.field_cap)
  fp=mobility_factor(cp,F,self.p.sigma_p_eV,self.p.T,a,self.mode,self.gamma_p,self.c_cut,self.field_cap)
  o['Jn']*=fn;o['Jp']*=fp
  o.update(field_abs_Vcm=F,mobility_factor_n=fn,mobility_factor_p=fp,c_edge_n=cn,c_edge_p=cp)
  return o
 def ledger(self,z,V,L):
  r=super().ledger(z,V,L);o=self.evaluate(z,V)
  r.update(transport_mode=self.mode,density_cutoff=self.c_cut,density_capped_length_fraction_n=float(np.sum(self.dx*(o['c_edge_n']>self.c_cut))/self.p.d_cm),density_capped_length_fraction_p=float(np.sum(self.dx*(o['c_edge_p']>self.c_cut))/self.p.d_cm),field_law_status='warm-window model candidate; SCLC fit quality and microscopic model limitations separately reported')
  return r
