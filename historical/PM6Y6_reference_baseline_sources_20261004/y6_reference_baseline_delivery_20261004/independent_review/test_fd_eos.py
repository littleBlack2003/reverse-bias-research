"""Independent Gauss-Fermi EOS tests, including integration and inversion.

Independent reference: adaptive QUADPACK integration over [-20,20], split at
Fermi edge and Gaussian moment saddle; no interpolant used by the reference.
"""
from pathlib import Path
import sys,json,hashlib,time
import numpy as np
from scipy.special import expit, roots_hermitenorm
from scipy.integrate import quad
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from fd_eos import GaussianDOS,KB,direct_eos
out={'tests':[],'domain':'100–300 K, sigma_e=60 meV, sigma_h=74 meV; auxiliary sigma=0 and numerical s=8.7 edge'}
def check(name,val,tol):
 val=float(val);out['tests'].append({'name':name,'error':val,'tolerance':tol});assert np.isfinite(val) and val<tol,(name,val,tol)
def independent(eta,s):
 def f(z,derivative=False):
  x=eta-s*z
  return np.exp(-z*z/2)/np.sqrt(2*np.pi)*expit(x)*(expit(-x) if derivative else 1.)
 pts=sorted(set([-20.,20.,0.]+([max(-19.99,min(19.99,eta/s)),-s] if s else [])))
 return np.array([sum(quad(lambda z:f(z,d),lo,hi,epsabs=1e-280,epsrel=3e-12,limit=300)[0] for lo,hi in zip(pts[:-1],pts[1:])) for d in [False,True]])
start=time.time()
for T in [100.,120.,150.,200.,223.,250.,300.]:
 for sigma in [.060,.074]:
  dos=GaussianDOS(sigma,T);s=dos.s
  etas=np.array([-300.,-239.993,-160.013,-100.007,-s*s-0.003,-60.013,-40.017,-20.001,-7.009,-2.003,-0.013,0.,.013,7.009,20.001,60.013])
  refs=np.array([independent(eta,s) for eta in etas]);c,dc,g=dos.evaluate(etas)
  check(f'c_adaptive_{T}_{sigma}',np.max(abs(c-refs[:,0])/np.maximum(refs[:,0],1e-290)),2e-10)
  check(f'dc_adaptive_{T}_{sigma}',np.max(abs(dc-refs[:,1])/np.maximum(refs[:,1],1e-290)),1e-7)
  check(f'particle_hole_{T}_{sigma}',np.max(abs(dos.evaluate(-etas)[0]+c-1)),3e-16)
  # Test inverse by direct independent quadrature, rather than same-table roundtrip.
  targets=np.r_[np.logspace(-16,-2,15),.1,.5,.9,1-1e-10]
  roots=dos.eta_from_c(targets)
  reference_c=np.array([independent(eta,s)[0] for eta in roots])
  check(f'inverse_direct_{T}_{sigma}',np.max(abs(reference_c-targets)/targets),2e-10)
  check(f'log_roundtrip_{T}_{sigma}',np.max(abs(dos.logc(dos.eta_from_logc(np.array([-1000.,-100.,-20.,-2.,-np.log(2),-.01,-1e-16])))-np.array([-1000.,-100.,-20.,-2.,-np.log(2),-.01,-1e-16]))),1e-12)
  grid=np.linspace(-240,120,48013);lc=dos.logc(grid);dlc=dos.dlogc(grid)
  assert np.all(np.diff(lc)>0),('strict monotonicity',T,sigma)
  assert np.all(dlc>0) and np.all(dlc<=1+1e-12),('positive susceptibility',T,sigma)
  # Check derivative consistency; points avoid numerical differencing in saturated c.
  et=np.array([-80.003,-40.007,-20.017,-5.003,.007,20.013]);h=1e-4
  dfd=(dos.logc(et+h)-dos.logc(et-h))/(2*h)
  check(f'derivative_consistency_{T}_{sigma}',np.max(abs(dfd/dos.dlogc(et)-1)),1e-7)
  same=dos.secant_g(et,et);_,_,g=dos.evaluate(et)
  check(f'secant_equal_state_{T}_{sigma}',np.max(abs(same/g-1)),1e-15)
  check(f'secant_symmetry_{T}_{sigma}',np.max(abs(dos.secant_g(et-.03,et+.03)/dos.secant_g(et+.03,et-.03)-1)),1e-15)
for s in [0.,8.7]:
 dos=GaussianDOS(s*KB,1.);etas=np.array([-220.013,-70.017,-35.013,-2.003,0.,2.003,35.013])
 c,dc,g=dos.evaluate(etas);refs=np.array([independent(eta,s) for eta in etas])
 check(f'edge_c_{s}',np.max(abs(c-refs[:,0])/np.maximum(refs[:,0],1e-290)),2e-10)
 check(f'edge_dc_{s}',np.max(abs(dc-refs[:,1])/np.maximum(refs[:,1],1e-290)),1e-7)
 if s==0:check('single_level_fermi',np.max(abs(c-expit(etas))),2e-16)
for bad in [0.,1.,-1.,np.nan,np.inf]:
 try:GaussianDOS(.060,100.).eta_from_c(bad)
 except ValueError:pass
 else:raise AssertionError('invalid occupation accepted')
try:GaussianDOS(.080,100.)
except ValueError:pass
else:raise AssertionError('unvalidated EOS range accepted')
out['wall_seconds']=time.time()-start
out['source_sha256']=hashlib.sha256((ROOT/'code'/'fd_eos.py').read_bytes()).hexdigest()
(ROOT/'independent_review'/'fd_eos_results.json').write_text(json.dumps(out,indent=2))
print('PASS',len(out['tests']),'tests;',out['wall_seconds'],'seconds')
for prefix in ['c_adaptive','dc_adaptive','inverse_direct','derivative_consistency']:
 print(prefix,max(r['error'] for r in out['tests'] if r['name'].startswith(prefix)))
