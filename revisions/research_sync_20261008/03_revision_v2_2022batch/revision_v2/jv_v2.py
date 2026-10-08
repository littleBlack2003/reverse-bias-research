import sys,time,json;sys.path.insert(0,'.')
from run_v2 import *
def jv(T,amp,ld,lb,zeta,N=321,vmax=1.45,dv=.05):
    d,z0=solve_state(T,amp,ld,lb,zeta)
    states={0.:z0};rows=[];z=z0;vprev=0.
    V=0.;Jprev=None;bracket=None
    while V<vmax:
        V1=round(V+dv,6)
        z,info=advance(d,z,V,V1,1.,1.);states[V1]=z.copy()
        o=d.evaluate(z,V1);J=float(np.mean(o['Jn']+o['Jp']));rows.append((V1,J))
        if J>0 and bracket is None:bracket=(V,V1);break
        V=V1
    o=d.evaluate(z0,0.);rows.insert(0,(0.,float(np.mean(o['Jn']+o['Jp']))))
    cache={}
    def Jt(v):
        k=round(v,9)
        if k in cache:return cache[k][0]
        lo=max(s for s in states if s<=v+1e-12);zz,_=advance(d,states[lo],lo,v,1.,1.);o=d.evaluate(zz,v);cache[k]=(float(np.mean(o['Jn']+o['Jp'])),zz);return cache[k][0]
    vc=brentq(Jt,*bracket,xtol=1e-7)
    vs=np.linspace(0,vc,41);P=[-v*Jt(float(v)) if v>0 else 0. for v in vs]
    # refine Pmax
    i=int(np.argmax(P));lo,hi=vs[max(i-1,0)],vs[min(i+1,40)]
    opt=minimize_scalar(lambda v:v*Jt(float(v)),bounds=(lo,hi),method='bounded',options={'xatol':1e-6})
    Jsc=-rows[0][1];Pm=-opt.fun
    return dict(T=T,Jsc_mAcm2=Jsc*1e3,Voc_V=float(vc),FF=float(Pm/(Jsc*vc)),Vmp=float(opt.x)),rows
out={}
for T in (300,100):
  for amp in ('cal','one'):
    zeta=calibrate_zeta(amp)
    for name,(ld,lb) in {'V1':(1.,0.),'V2':(1.,1.)}.items():
        t=time.time()
        try:
            m,rows=jv(T,amp,ld,lb,zeta);out[f'{name}_{amp}_{T}']=m;print(name,amp,T,{k:round(v,4) for k,v in m.items()},f'{time.time()-t:.0f}s',flush=True)
        except Exception as e:print(name,amp,T,'FAIL',repr(e)[:200],flush=True)
json.dump(out,open('results/jv_metrics.json','w'),indent=1)
