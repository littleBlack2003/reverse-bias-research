"""Independent closure and reaction audit, no source files modified.
Run: python independent_review/review.py. Results include nonzero discretization errors.
"""
from pathlib import Path
import sys, json, hashlib
import numpy as np
from scipy.integrate import quad
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from gaussian_transport import *
from gaussian_device import GaussianDevice
out={'tests':[], 'finite_gradient_errors':[], 'saved_state_reactions':[]}
def record(name,value,limit):
 value=float(value);assert np.isfinite(value) and value<limit,(name,value,limit)
 out['tests'].append(dict(name=name,value=value,limit=limit))
for car in ['n','p']:
 psi=np.linspace(-20,20,161);eta=-60+(psi if car=='n' else -psi)
 j=inverse_activity_flux(eta,psi,300,1e21,.1,1e-7,1e-4,car)
 record('default_affinity_equilibrium_'+car,np.max(abs(j)),1e-25)
 eta=np.array([-100.,-99.]);psi=np.array([0.,.6]);s=.1/(KB*300);c=np.exp(eta+s*s/2);pref=Q*1e-4*KB*300*1e21/1e-7
 exact=pref*(c[1]*bernoulli(.6)-c[0]*bernoulli(-.6)) if car=='n' else pref*(c[0]*bernoulli(.6)-c[1]*bernoulli(-.6))
 val=inverse_activity_flux(eta,psi,300,1e21,.1,1e-7,1e-4,car)[0]
 record('dilute_SG_'+car,abs(val/exact-1),1e-11)
for h in [2.,1.,.5,.25]:
 eta=np.array([-9-h/2,-9+h/2]);approx=gaussian_eos(-9,4)[0]*2*np.sinh(h/2)
 exact=quad(lambda x:float(gaussian_eos(x,4,1024)[0]),*eta,epsabs=1e-13)[0]
 out['finite_gradient_errors'].append(dict(delta_eta=h,relative_error=float(approx/exact-1)))
errors=np.array([r['relative_error'] for r in out['finite_gradient_errors']]);orders=np.log(errors[:-1]/errors[1:])/np.log(2)
record('finite_gradient_second_order',max(abs(orders-2)),.04)
for s in [2.,4.,6.]:
 eta=np.linspace(-60,30,601);a,b,_=gaussian_eos(eta,s,256);aa,bb,_=gaussian_eos(eta,s,2048)
 record('EOS_256_2048_density_s'+str(s),max(abs(a-aa)/np.maximum(aa,1e-290)),5e-6)
 record('EOS_256_2048_derivative_s'+str(s),max(abs(b-bb)/np.maximum(bb,1e-290)),5e-6)
# Audit saved solutions, independently checking each reservoir detailed-balance ratio,
# cycle identity, direct master-equation current, positive entropy, and quadrature.
for T in [300,340,380,420]:
 for V in [0,-5,-15]:
  f=ROOT/'data'/f'gaussian_T{T}_N81_V{V}.npz'
  if not f.exists():continue
  z=np.load(f)['z'];d=GaussianDevice(81,T=T);o=d.evaluate(z,V)
  lh=np.log(o['aH']/o['bH']);ll=np.log(o['aL']/o['bL'])
  errH=max(abs(lh-(o['mh']-o['et'])/d.vt));errL=max(abs(ll-(o['ml']-o['et'])/d.vt))
  record(f'port_LDB_H_{T}_{V}',errH,1e-10);record(f'port_LDB_L_{T}_{V}',errL,1e-10)
  raw=(o['aH']*o['bL']-o['aL']*o['bH'])/(o['aH']+o['bH']+o['aL']+o['bL'])
  scale=np.maximum((o['aH']*o['bL']+o['aL']*o['bH'])/(o['aH']+o['bH']+o['aL']+o['bL']),1e-290)
  record(f'current_direct_master_{T}_{V}',max(abs(raw-o['j'])/scale),1e-10)
  A=V/d.vt+z[:,1]+z[:,2]
  record(f'bim_entropy_sign_{T}_{V}',max(0.,float(-np.min(o['bim']*A))),1e-20)
  ref=GaussianDevice(81,T=T,order=512).evaluate(z,V)
  qerr=max(float(np.max(abs(o[k]-ref[k])/np.maximum(abs(ref[k]),1e-290))) for k in ['aH','bH','aL','bL'])
  record(f'rate_quad128_512_{T}_{V}',qerr,1e-4)
  # A=0 requires common local chemical potential, regardless of bent bands.
  zeq=z.copy();zeq[:,1:]=0;eq=d.evaluate(zeq,0)
  record(f'equilibrium_bim_{T}_{V}',max(abs(eq['bim'])),1e-25)
  record(f'equilibrium_generation_{T}_{V}',max(abs(eq['j'])),1e-25)
  out['saved_state_reactions'].append(dict(T=T,V=V,port_H_error=float(errH),port_L_error=float(errL),rate_quadrature_error=qerr))
out['source_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'code'/'gaussian_transport.py',ROOT/'code'/'gaussian_device.py']}
out['scope']='Formula, equilibrium, quadrature and saved-state reaction audit; not material calibration or full temperature curve validation.'
(ROOT/'independent_review'/'results.json').write_text(json.dumps(out,indent=2))
print(len(out['tests']),'checks passed; finite-gradient error retained as an explicit discretization limitation')
