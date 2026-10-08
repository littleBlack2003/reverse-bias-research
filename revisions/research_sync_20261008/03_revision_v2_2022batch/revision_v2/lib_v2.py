"""lib_v2.py: refit-parameter devices and a robust J-V (fine voltage grid, fallback stepping)."""
import sys,json,dataclasses;sys.path.insert(0,'.')
from run_v2 import *
from extend_device import temperature_continue
from scipy.optimize import brentq,minimize_scalar
def load(tag):
    R=json.load(open(f'results/recal_{tag}.json'));return R,tag.startswith('V2'),tag.split('_')[1]
def dev(tag,T,t=1.):
    R,lb,amp=load(tag);lnG,dEg,lnS,lnM=R['x'];zeta0=R['zeta0'];p0=parameters(float(T))[0]
    p=dataclasses.replace(p0,G_cm3s=p0.G_cm3s*np.exp(t*lnG),Eg_eV=p0.Eg_eV+t*dEg,beta_cm3s=p0.beta_cm3s*(1. if lb else np.exp(t*lnS)))
    a=AMP[amp];return ClosureDevice(p,321,amp_n=a[0]**t*np.exp(t*lnM),amp_p=a[1]**t*np.exp(t*lnM),lam_d=t,lam_b=float(lb)*t,zeta=(zeta0*np.exp(t*lnS) if lb else None))
def orig_states(Ts):
    z=saved(300)['sc'][2];tr=[];prev=300;out={300:z}
    for T in Ts:
        z=temperature_continue(z,prev,T,0.,1.,321,tr);out[T]=z;prev=T
    return out
def to_final(tag,T,z):
    t=0.;dt=.25
    while t<1-1e-12:
        st=min(dt,1-t);d=dev(tag,T,t+st);zz,info=solve(d,0.,z,1.,maxiter=300)
        if info['success']:z=zz;t+=st;dt=min(dt*1.5,.5)
        else:
            dt=st/2
            if dt<1e-4:raise RuntimeError('stalled')
    return d,z
import device as _dev
_strict_solve=_dev.solve
def loose_solve(device,V,z,light=0.,maxiter=160):
    """same Newton iteration as device.solve, but accepts on the scaled-residual norm alone (<3e-10),
    i.e. WITHOUT the relative current-accuracy gate, which cannot be met when J -> 0 at open circuit."""
    z=z.copy()
    for it in range(maxiter):
        r=device.residual(z,V,light);jac=device.jacobian(z,V,light)
        scale=np.maximum(np.asarray(abs(jac).max(axis=1).toarray()).ravel(),1e-250);norm=float(max(abs(r)/scale))
        if norm<3e-10:
            audit=device.ledger(z,V,light);return z,dict(success=True,iterations=it,norm=norm,**audit)
        from scipy.sparse import diags;from scipy.sparse.linalg import spsolve
        step=spsolve(diags(1/scale)@jac,-r/scale).reshape(z.shape)
        if not np.all(np.isfinite(step)):return z,dict(success=False,reason='nonfinite',norm=norm)
        alpha=min(1.,4/max(np.max(abs(step[:,1:])),1e-30),8/max(np.max(abs(step[:,0])),1e-30))
        for ls in range(32):
            trial=z+alpha*step
            try:nn=float(max(abs(device.residual(trial,V,light))/scale))
            except (ValueError,FloatingPointError):nn=np.inf
            if nn<norm:z=trial;break
            alpha*=.5
        else:return z,dict(success=False,reason='line search stalled',norm=norm)
    return z,dict(success=False,reason='iteration limit',norm=norm)
def Jmean(d,z,V):
    o=d.evaluate(z,V);return float(np.mean(o['Jn']+o['Jp']))
def adv(d,z,V0,V1,sub=4):
    try:
        _dev.solve=_strict_solve;return advance(d,z,V0,V1,1.,1.)[0]
    except RuntimeError:
        _dev.solve=loose_solve
        try:return advance(d,z,V0,V1,1.,1.)[0]
        finally:_dev.solve=_strict_solve
def jv(d,z0,dv=.02,vmax=1.5):
    Jsc=-Jmean(d,z0,0.);states={0.:z0};z=z0;V=0.;br=None
    while V<vmax:
        V1=round(V+dv,6);z=adv(d,z,V,V1);states[V1]=z.copy();J=Jmean(d,z,V1)
        if J>0:br=(V,V1);break
        V=V1
    cache={}
    def Jt(v):
        k=round(v,10)
        if k in cache:return cache[k]
        lo=max(q for q in states if q<=v+1e-12);zz=adv(d,states[lo],lo,v);cache[k]=Jmean(d,zz,v);return cache[k]
    vc=brentq(Jt,*br,xtol=2e-8)
    r=minimize_scalar(lambda v:v*Jt(float(v)),bounds=(.3*vc,vc),method='bounded',options={'xatol':1e-6})
    return vc,-r.fun/(Jsc*vc),Jsc*1e3
