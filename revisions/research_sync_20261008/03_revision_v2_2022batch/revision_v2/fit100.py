"""fit100.py: at 100 K only, free (beta scale, mobility multiplier) with G, Eg fixed at the 300 K refit; targets = 2022 batch 100 K
(Jsc 17.30, Voc 0.959, FF 40.2 %). Three targets, two parameters -> a consistency test of the closure."""
import sys,json,time,dataclasses;sys.path.insert(0,'.')
from lib_v2 import *
from scipy.optimize import least_squares
tag=sys.argv[1];TG=np.array([17.30,.959,.402]);SC=np.array([.2,.005,.01])
zo=orig_states([275,250,225,200,175,150,125,100])
d0,z0=to_final(tag,100,zo[100])
def mk(x):
    p=dataclasses.replace(d0.p,beta_cm3s=d0.p.beta_cm3s*np.exp(x[0]))
    return ClosureDevice(p,321,amp_n=d0.amp_n*np.exp(x[1]),amp_p=d0.amp_p*np.exp(x[1]),lam_d=d0.lam_d,lam_b=d0.lam_b,zeta=(None if d0.zeta is None else d0.zeta*np.exp(x[0])))
st=dict(x=np.zeros(2),z=z0,n=0)
def at(x):
    x=np.asarray(x,float);xc=st['x'];z=st['z'];t=0.;dt=1.
    while t<1-1e-12:
        s=min(dt,1-t);d=mk(xc+(t+s)*(x-xc));zz,info=solve(d,0.,z,1.,maxiter=250)
        if info['success']:z=zz;t+=s;dt=min(dt*2,1.)
        else:
            dt=s/2
            if dt<1e-5:raise RuntimeError('stalled')
    st['x']=x.copy();st['z']=z;return d,z
def resid(x):
    d,z=at(x);vc,FF,Jsc=jv(d,z,dv=.04);m=np.array([Jsc,vc,FF]);st['n']+=1
    r=np.r_[(m-TG)/SC,x[1]/1.5]
    print(f'  eval {st["n"]:2d} lnS={x[0]:+.3f} lnM={x[1]:+.3f} Jsc={Jsc:.3f} Voc={vc:.4f} FF={100*FF:.1f}% |r|={np.linalg.norm(r[:3]):.2f}',flush=True);return r
t0=time.time();sol=least_squares(resid,np.array([5.0,0.5]),bounds=([-8,-8],[9,8]),x_scale=[1.,1.],diff_step=5e-2,max_nfev=14,xtol=1e-3,ftol=1e-4)
d,z=at(sol.x);vc,FF,Jsc=jv(d,z,dv=.02);o=d.evaluate(z,0.);bulk=(d.x>.2)&(d.x<.8)
mun=float(np.median((d.mn*o['mobility_factor_n'])[bulk[:-1]]));mup=float(np.median((d.mp*o['mobility_factor_p'])[bulk[:-1]]))
out=dict(tag=tag,x=list(sol.x),beta_scale=float(np.exp(sol.x[0])),beta100=float(d.p.beta_cm3s),M100=float(np.exp(sol.x[1])),Jsc=Jsc,Voc=vc,FF=FF,mu_n_bulk=mun,mu_p_bulk=mup,
 zeta_eff=float(d.p.beta_cm3s*d.eps/(Q*(mun+mup))),target=list(TG),seconds=time.time()-t0)
json.dump(out,open(f'results/fit100_{tag}.json','w'),indent=1);print('DONE',json.dumps(out),flush=True)
