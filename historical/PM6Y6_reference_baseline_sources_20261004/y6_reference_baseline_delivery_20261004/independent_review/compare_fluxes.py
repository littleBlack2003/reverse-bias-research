"""Independent normalized flux tests: secant-Einstein SG vs inverse activity.
The two schemes need not match at finite mesh. Both must recover the SAME
Gaussian-FD local-equilibrium continuum law and zero equilibrium current.
Currents here are divided by q*mu*kT*N; x is dimensionless.
"""
from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from fd_eos import GaussianDOS,KB

def logB(x):
 x=np.asarray(x);ans=np.empty_like(x);pos=x>50;neg=x<-50;small=abs(x)<1e-5;mid=~(pos|neg|small)
 ans[pos]=np.log(x[pos])-x[pos]-np.log1p(-np.exp(-x[pos]))
 ans[neg]=np.log(-x[neg])-np.log1p(-np.exp(x[neg]))
 ans[small]=np.log1p(-x[small]/2+x[small]**2/12-x[small]**4/720)
 ans[mid]=np.log(x[mid]/np.expm1(x[mid]));return ans

def flux(eta,psi,dx,dos,car='n',scheme='secant',affinity=None):
 de=np.diff(eta);dp=np.diff(psi);sgn=1 if car=='n' else -1
 A=de-sgn*dp if affinity is None else np.asarray(affinity)
 if scheme=='secant':
  g=dos.secant_g(eta[:-1],eta[1:]);a=A/g
  lpre=np.log(g)+dos.logc(eta[:-1])+logB(-sgn*dp/g)
 else:
  a=A;lpre=dos.logc((eta[:-1]+eta[1:])/2)-de/2+logB(-sgn*dp)
 aa=abs(a);lexp=np.full_like(a,-np.inf);nonzero=aa>0
 lexp[nonzero]=np.maximum(a[nonzero],0)+np.log(-np.expm1(-aa[nonzero]))
 return sgn*np.sign(a)*np.exp(lpre+lexp-np.log(dx))

out={'finite_mesh':[],'tests':[]}
def check(name,val,tol):
 val=float(val);out['tests'].append({'name':name,'error':val,'tolerance':tol});assert np.isfinite(val) and val<tol,(name,val,tol)
for T in [100,150,200,300]:
 for sig in [.060,.074]:
  d=GaussianDOS(sig,T)
  for car in ['n','p']:
   sgn=1 if car=='n' else -1
   for scheme in ['secant','inverse_activity']:
    psi=np.linspace(-80,80,1001);eta=-100+sgn*psi
    j=flux(eta,psi,.001,d,car,scheme,affinity=np.zeros(1000))
    check(f'exact_equilibrium_{T}_{sig}_{car}_{scheme}',np.max(abs(j)),1e-300)
    # Arbitrary far-from-equilibrium faces must obey flux-affinity sign.
    rng=np.random.default_rng(528);eta=rng.uniform(-160,30,1001);psi=rng.uniform(-100,100,1001)
    A=np.diff(eta)-sgn*np.diff(psi)
    j=flux(eta,psi,.001,d,car,scheme)
    check(f'flux_entropy_{T}_{sig}_{car}_{scheme}',max(0.,float(-np.min(sgn*j*A))),1e-300)
  for scheme in ['secant','inverse_activity']:
   rows=[]
   for n in [41,81,161,321,641]:
    x=np.linspace(0,1,n);h=1/(n-1);xm=(x[:-1]+x[1:])/2
    eta=-30+8*np.sin(2*np.pi*x);psi=3*np.cos(2*np.pi*x)
    c=d.evaluate(-30+8*np.sin(2*np.pi*xm))[0]
    exact=c*(16*np.pi*np.cos(2*np.pi*xm)+6*np.pi*np.sin(2*np.pi*xm))
    j=flux(eta,psi,h,d,scheme=scheme)
    err=np.linalg.norm(j-exact)/np.linalg.norm(exact)
    rows.append({'N':n,'relative_L2_error':float(err)})
   orders=np.log2(np.array([z['relative_L2_error'] for z in rows[:-1]])/np.array([z['relative_L2_error'] for z in rows[1:]]))
   check(f'second_order_{T}_{sig}_{scheme}',abs(orders[-1]-2),.02)
   out['finite_mesh'].append({'T':T,'sigma_eV':sig,'scheme':scheme,'rows':rows,'orders':orders.tolist()})
(ROOT/'independent_review'/'flux_comparison_results.json').write_text(json.dumps(out,indent=2))
print('PASS',len(out['tests']),'flux checks; both schemes converge at order2')
for r in out['finite_mesh']:
 if r['T']==100 and r['sigma_eV']==.074: print(r)
