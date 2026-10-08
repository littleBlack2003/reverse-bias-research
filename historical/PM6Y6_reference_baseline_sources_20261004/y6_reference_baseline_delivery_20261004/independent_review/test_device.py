"""Adversarial finite-volume implementation checks (independent reviewer)."""
from pathlib import Path
import sys,json,hashlib,time
from dataclasses import replace
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'code'))
from device import *
sys.path.insert(0,str(ROOT/'independent_review'))
# Do not import compare_fluxes: importing it executes its test program.
out={'tests':[],'solved_equilibria':[],'scope':'Implementation checks, not material validation'}
def check(name,val,tol):
 val=float(val);out['tests'].append({'name':name,'error':val,'tolerance':tol});assert np.isfinite(val) and val<tol,(name,val,tol)
for T in [100,200,300]:
 for contact in ['blocking','equilibrium_reservoir']:
  p=Parameters(T=T,G_cm3s=1e22,minority_contact=contact,allow_mobility_extrapolation=True)
  d=Device(p,nodes=41)
  z=d.initial();o=d.evaluate(z,0)
  check(f'exact_equilibrium_J_{T}_{contact}',max(abs(o['Jn']).max(),abs(o['Jp']).max()),1e-300)
  check(f'exact_equilibrium_R_{T}_{contact}',max(abs(o['R'])),1e-300)
  # Arbitrary smooth, finite nonequilibrium state: verify telescoping ledger
  # identity WITH residual correction; it must hold before Newton converges.
  x=d.x;V=.37;z=d.initial(V);z[:,0]+=.1*np.sin(np.pi*x)
  z[:,1]+=1.3*np.sin(np.pi*x);z[:,2]+=.8*np.sin(2*np.pi*x)
  oo=d.evaluate(z,V);r=d.residual(z,V,.8).reshape(z.shape)
  td=np.dot(oo['Jn'],np.diff(z[:,1])*d.vt)-np.dot(oo['Jp'],np.diff(z[:,2])*d.vt)
  rd=Q*np.dot(d.vol,(oo['R']*oo['A']*d.vt)[1:-1]);ld=Q*.8*p.G_cm3s*np.dot(d.vol,oo['A'][1:-1]*d.vt)
  un=z[:,1]+(p.barrier_p_eV+V)/d.vt; vp=z[:,2]-p.barrier_p_eV/d.vt
  boundary=d.vt*(oo['Jn'][-1]*un[-1]-oo['Jn'][0]*un[0]-oo['Jp'][-1]*vp[-1]+oo['Jp'][0]*vp[0])
  residual_work=d.vt*(-np.dot(un[1:-1],r[1:-1,1])+np.dot(vp[1:-1],r[1:-1,2]))
  check(f'energy_telescope_{T}_{contact}',abs(td+rd-ld-boundary-residual_work)/max(abs(td),abs(rd),abs(ld),abs(boundary),abs(residual_work),1e-100),3e-14)
  check(f'transport_entropy_{T}_{contact}',max(0.,-td),1e-300)
  check(f'reaction_entropy_{T}_{contact}',max(0.,-rd),1e-300)
  # Manufactured smooth continuum flux, both species and high-disorder low T.
  for car in ['n','p']:
   sgn=1 if car=='n' else -1;dos=d.en if car=='n' else d.ep;mu=d.mn if car=='n' else d.mp
   errors=[]
   for n in [81,161,321,641]:
    xx=np.linspace(0,1,n);xm=(xx[:-1]+xx[1:])/2;eta=-30+8*np.sin(2*np.pi*xx);psi=3*np.cos(2*np.pi*xx)
    qf=eta-sgn*psi
    cur,_=sg_flux(eta,psi,qf,dos.logc,dos.dlogc,mu,p.N0_cm3,d.vt,1/(n-1),car)
    expected=sgn*Q*mu*p.N0_cm3*d.vt*dos.evaluate(-30+8*np.sin(2*np.pi*xm))[0]*(16*np.pi*np.cos(2*np.pi*xm)+sgn*6*np.pi*np.sin(2*np.pi*xm))
    errors.append(np.linalg.norm(cur-expected)/np.linalg.norm(expected))
   check(f'sg_continuum_order_{T}_{contact}_{car}',abs(np.log2(errors[-2]/errors[-1])-2),.02)
  # Derivative coloring includes precisely the nearest-neighbor dependence.
  rng=np.random.default_rng(712);v=rng.normal(size=z.shape);v/=np.linalg.norm(v)
  eps=2e-5;jac=d.jacobian(z,V,.8)
  direct=(d.residual(z+eps*v,V,.8)-d.residual(z-eps*v,V,.8))/(2*eps)
  approx=jac@v.ravel();check(f'colored_jacobian_{T}_{contact}',np.linalg.norm(direct-approx)/np.linalg.norm(direct),2e-7)
  zeq,info=solve(d,0,d.initial(),0,maxiter=120)
  assert info['success'],('equilibrium solve failed',T,contact,info)
  assert info['gate_passed']
  out['solved_equilibria'].append({k:v for k,v in info.items() if k!='history'})
  check(f'equilibrium_charge_{T}_{contact}',info['charge_relative'],1e-9)
  np.savez_compressed(ROOT/'independent_review'/f'equilibrium_T{T}_{contact}.npz',z=zeq)
# Unmeasured mobility range must require explicit opt-in.
try: Device(Parameters(T=200))
except ValueError: pass
else: raise AssertionError('unapproved mobility extrapolation accepted')
# Uniform mesh API and invalid physical scales.
assert np.allclose(np.diff(Device(nodes=9,mesh_strength=0).x),1/8)
for kwargs in [{'T':0},{'d_cm':0},{'N0_cm3':0},{'mu_n_300':-1},{'epsilon_r':0},{'beta_cm3s':0},{'G_cm3s':-1}]:
 try:Device(replace(Parameters(),**kwargs))
 except ValueError:pass
 else:raise AssertionError(('invalid physical scale accepted',kwargs))
out['sha256']={n:hashlib.sha256((ROOT/'code'/n).read_bytes()).hexdigest() for n in ['device.py','fd_eos.py']}
(ROOT/'independent_review'/'device_results.json').write_text(json.dumps(out,indent=2))
print('PASS',len(out['tests']),'device checks and',len(out['solved_equilibria']),'dark equilibria')
