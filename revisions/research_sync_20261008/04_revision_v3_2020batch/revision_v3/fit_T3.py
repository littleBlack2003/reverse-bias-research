"""fit_T3.py (3 params incl. generation scale) -- derived from fit_T.py: TAG T1,T2,...  : per-temperature inversion on the 2020 batch. G, Eg fixed at the 300 K refit (tag);
free: ln(beta scale vs. Arrhenius extrapolation), ln(mobility multiplier). 3 targets (Jsc,Voc,FF), 2 params + weak ridge."""
import sys,json,time,dataclasses,csv;sys.path.insert(0,'.')
from lib_v2 import *
from scipy.optimize import least_squares
tag=sys.argv[1];Tl=[float(x) for x in sys.argv[2].split(',')]
D={float(r['T']):r for r in csv.DictReader(open('data2020.csv'))}
SC=np.array([.2,.005,.01])
ALL=[280,250,230,200,175,150,125,100]
zo=orig_states(ALL)
def run(T):
    TG=np.array([float(D[T]['Jsc']),float(D[T]['Voc']),float(D[T]['FF'])/100])
    d0,z0=to_final(tag,T,zo[int(T)])
    def mk(x):
        p=dataclasses.replace(d0.p,beta_cm3s=d0.p.beta_cm3s*np.exp(x[0]),G_cm3s=d0.p.G_cm3s*np.exp(x[2]))
        return ClosureDevice(p,321,amp_n=d0.amp_n*np.exp(x[1]),amp_p=d0.amp_p*np.exp(x[1]),lam_d=d0.lam_d,lam_b=d0.lam_b,zeta=None)
    st=dict(x=np.zeros(3),z=z0,n=0)
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
        r=np.r_[(m-TG)/SC,x[1]/3.]
        print(f'  T={T:g} eval {st["n"]:2d} lnS={x[0]:+.3f} lnM={x[1]:+.3f} lnG={x[2]:+.3f} Jsc={Jsc:.3f} Voc={vc:.4f} FF={100*FF:.1f}% |r|={np.linalg.norm(r[:3]):.2f}',flush=True);return r
    t0=time.time()
    def pj(Tq):
        o=json.load(open(f'results/fitT_{tag}_{int(Tq)}.json'));return np.array(o['x']),o['Jsc']
    if T==200:
        xa,_=pj(175);xb,_=pj(230);xp=(xa+xb)/2;Jp=22.6
    else:xp,Jp=pj(T)
    x0=np.r_[xp,np.log(TG[0]/Jp)]
    sol=least_squares(resid,x0,bounds=([-8,-8,-1.5],[10,9,.3]),x_scale=[1.,1.,.3],diff_step=5e-2,max_nfev=14,xtol=1e-3,ftol=1e-4)
    d,z=at(sol.x);vc,FF,Jsc=jv(d,z,dv=.02);o=d.evaluate(z,0.);bulk=(d.x>.2)&(d.x<.8)
    mun=float(np.median((d.mn*o['mobility_factor_n'])[bulk[:-1]]));mup=float(np.median((d.mp*o['mobility_factor_p'])[bulk[:-1]]))
    out=dict(tag=tag,T=T,x=list(sol.x),beta_scale=float(np.exp(sol.x[0])),beta=float(d.p.beta_cm3s),G=float(d.p.G_cm3s),G_ratio=float(np.exp(sol.x[2])),M=float(np.exp(sol.x[1])),Jsc=Jsc,Voc=vc,FF=FF,
        mu_n_bulk=mun,mu_p_bulk=mup,zeta_eff=float(d.p.beta_cm3s*d.eps/(Q*(mun+mup))),target=list(TG),seconds=time.time()-t0)
    json.dump(out,open(f'results/fitT3_{tag}_{int(T)}.json','w'),indent=1);print('DONE',json.dumps(out),flush=True)
for T in Tl:
    try:run(T)
    except Exception as e:print('FAIL',T,repr(e),flush=True)
