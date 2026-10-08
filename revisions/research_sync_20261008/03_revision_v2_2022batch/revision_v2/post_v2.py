"""post_v2.py: after the 300 K refit, evaluate Jsc/Voc/FF versus T for one variant (refit parameters applied at ALL T)."""
import sys,json,dataclasses,time;sys.path.insert(0,'.')
from run_v2 import *
from extend_device import temperature_continue
from scipy.optimize import brentq,minimize_scalar
tag=sys.argv[1];R=json.load(open(f'results/recal_{tag}.json'));lb=1. if tag.startswith('V2') else 0.;amp=tag.split('_')[1]
lnG,dEg,lnS,lnM=R['x'];zeta0=R['zeta0']
def dev(T,t):
    p0=parameters(float(T))[0]
    p=dataclasses.replace(p0,G_cm3s=p0.G_cm3s*np.exp(t*lnG),Eg_eV=p0.Eg_eV+t*dEg,beta_cm3s=p0.beta_cm3s*(1. if lb else np.exp(t*lnS)))
    a=AMP[amp]
    return ClosureDevice(p,321,amp_n=a[0]**t*np.exp(t*lnM),amp_p=a[1]**t*np.exp(t*lnM),lam_d=t,lam_b=lb*t,zeta=(zeta0*np.exp(t*lnS) if lb else None))
def to_final(T,z):
    t=0.;dt=.25
    while t<1-1e-12:
        st=min(dt,1-t);d=dev(T,t+st);zz,info=solve(d,0.,z,1.,maxiter=300)
        if info['success']:z=zz;t+=st;dt=min(dt*1.5,.5)
        else:
            dt=st/2
            if dt<1e-4:raise RuntimeError('stalled')
    return d,z
def jv(d,z0):
    o=d.evaluate(z0,0.);Jsc=-float(np.mean(o['Jn']+o['Jp']));states={0.:z0};z=z0;V=0.;br=None
    while V<1.45:
        V1=round(V+.05,6);z,_=advance(d,z,V,V1,1.,1.);states[V1]=z.copy();o=d.evaluate(z,V1);J=float(np.mean(o['Jn']+o['Jp']))
        if J>0:br=(V,V1);break
        V=V1
    cache={}
    def Jt(v):
        k=round(v,10)
        if k in cache:return cache[k]
        lo=max(q for q in states if q<=v+1e-12);zz,_=advance(d,states[lo],lo,v,1.,1.);o=d.evaluate(zz,v);cache[k]=float(np.mean(o['Jn']+o['Jp']));return cache[k]
    vc=brentq(Jt,*br,xtol=2e-8);r=minimize_scalar(lambda v:v*Jt(float(v)),bounds=(.3*vc,vc),method='bounded',options={'xatol':1e-6})
    return vc,-r.fun/(Jsc*vc)
Ts=[300,275,250,225,200,175,150,125,100]
z=saved(300)['sc'][2];trace=[];zo={300:z};prev=300
for T in Ts[1:]:
    z=temperature_continue(z,prev,T,0.,1.,321,trace);zo[T]=z;prev=T
out={}
for T in Ts:
    row={}
    try:
        d,z=to_final(T,zo[T]);m=metrics_at(d,z);row.update(Jsc=m['Jsc_mAcm2'],R_over_G=m['R_over_G'],n_bulk=m['n_bulk_median'],mu_n=m['mu_n_bulk_med'],mu_p=m['mu_p_bulk_med'],gate=m['gate'])
        try:
            vc,FF=jv(d,z);row.update(Voc=vc,FF=FF)
        except Exception as e:row.update(Voc=None,FF=None,jv_error=repr(e)[:80])
    except Exception as e:row.update(error=repr(e)[:100])
    out[T]=row;print(tag,T,{k:(round(v,4) if isinstance(v,float) else v) for k,v in row.items()},flush=True)
json.dump(out,open(f'results/post_{tag}.json','w'),indent=1)
