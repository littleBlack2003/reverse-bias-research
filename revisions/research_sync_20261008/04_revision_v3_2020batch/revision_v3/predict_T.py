"""predict_T.py TAG : smooth (quadratic in w=(300-T)/200) functions for ln beta, ln M, ln(G/G300fit) trained on the 3-parameter
per-T inversions at TRAIN temperatures only; forward J-V at ALL 8 temperatures; held-out T = 175 (interpolation), 100 (extrapolation)."""
import sys,json,time,dataclasses,csv;sys.path.insert(0,'.')
from lib_v2 import *
tag=sys.argv[1];HOLD=[175,100];ALL=[280,250,230,200,175,150,125,100]
D={float(r['T']):r for r in csv.DictReader(open('data2020.csv'))}
F={T:json.load(open(f'results/fitT3_{tag}_{T}.json')) for T in ALL}
TR=[T for T in ALL if T not in HOLD]
w=lambda T:(300-np.asarray(T,float))/200.
def poly(key,deg=2):
    y=np.array([key(F[T]) for T in TR]);c=np.polyfit(w(TR),y,deg);return c
cb=poly(lambda o:np.log(o['beta']));cM=poly(lambda o:o['x'][1]);cG=poly(lambda o:o['x'][2])
zo=orig_states(ALL);res={}
for T in ALL:
    t0=time.time();d0,z0=to_final(tag,T,zo[T]);ww=w(T)
    xb=np.polyval(cb,ww)-np.log(d0.p.beta_cm3s);xM=np.polyval(cM,ww);xG=np.polyval(cG,ww);x=np.array([xb,xM,xG])
    # move from the base state to x (adaptive homotopy)
    def mk(x):
        p=dataclasses.replace(d0.p,beta_cm3s=d0.p.beta_cm3s*np.exp(x[0]),G_cm3s=d0.p.G_cm3s*np.exp(x[2]))
        return ClosureDevice(p,321,amp_n=d0.amp_n*np.exp(x[1]),amp_p=d0.amp_p*np.exp(x[1]),lam_d=d0.lam_d,lam_b=d0.lam_b,zeta=None)
    z=z0;xc=np.zeros(3);tt=0.;dt=.5
    while tt<1-1e-12:
        s=min(dt,1-tt);d=mk(xc+(tt+s)*(x-xc));zz,info=solve(d,0.,z,1.,maxiter=250)
        if info['success']:z=zz;tt+=s;dt=min(dt*2,1.)
        else:
            dt=s/2
            if dt<1e-5:raise RuntimeError('stalled')
    try:
        vc,FF,Jsc=jv(d,z,dv=.02);res[T]=dict(Jsc=Jsc,Voc=vc,FF=FF,heldout=T in HOLD,beta=float(d.p.beta_cm3s),Gratio=float(np.exp(xG)),M=float(np.exp(xM)))
    except Exception as e:res[T]=dict(error=repr(e),heldout=T in HOLD)
    r=res[T];tg=D[float(T)];print(T,'HELD' if T in HOLD else 'train',r if 'error' in r else 'Jsc %.2f/%s Voc %.4f/%s FF %.1f/%s'%(r['Jsc'],tg['Jsc'],r['Voc'],tg['Voc'],100*r['FF'],tg['FF']),'%.0fs'%(time.time()-t0),flush=True)
json.dump(dict(coef_lnbeta=list(cb),coef_lnM=list(cM),coef_lnG=list(cG),res=res),open(f'results/predict_{tag}.json','w'),indent=1)
