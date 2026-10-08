"""recal_v2.py: refit 300 K (G, Eg, beta-scale | zeta-scale) so a closure variant reproduces the frozen 300 K
Jsc / Voc / FF of the ORIGINAL model (targets = original model's own 300 K outputs, results/expected_metrics.json).
Original solver only (device.solve / advance)."""
import sys,json,time,dataclasses;sys.path.insert(0,'.')
from run_v2 import *
from scipy.optimize import least_squares,brentq,minimize_scalar
t=json.load(open('reference/expected_metrics.json'))['300']
TGT=np.array([t['Jsc_mAcm2'],t['Voc_V'],t['FF']])
SCALE=np.array([.001*TGT[0],.0005,.001])       # 0.1 % Jsc, 0.5 mV, 0.1 FF-point

def build_dev(amp,ld,lb,zeta0,x):
    p0=parameters(300.)[0];lnG,dEg,lnS,lnM=x
    p=dataclasses.replace(p0,G_cm3s=p0.G_cm3s*np.exp(lnG),Eg_eV=p0.Eg_eV+dEg,beta_cm3s=p0.beta_cm3s*(1. if lb else np.exp(lnS)))
    a=AMP[amp];return ClosureDevice(p,321,amp_n=a[0]*np.exp(lnM),amp_p=a[1]*np.exp(lnM),lam_d=ld,lam_b=lb,zeta=(zeta0*np.exp(lnS) if lb else None))

class Fit:
    def __init__(s,amp,ld,lb):
        s.amp,s.ld,s.lb=amp,ld,lb;s.zeta0=calibrate_zeta(amp) if lb else None
        s.x=np.zeros(4);d,z=solve_state(300,amp,ld,lb,s.zeta0);s.z=z;s.n=0
    def at(s,x):
        """SC state at parameters x, by adaptive interpolation from the current state"""
        x=np.asarray(x,float);xc=s.x.copy();z=s.z;tt=0.;dt=1.
        while tt<1-1e-12:
            st=min(dt,1-tt);xi=xc+(tt+st)*(x-xc);d=build_dev(s.amp,s.ld,s.lb,s.zeta0,xi)
            zz,info=solve(d,0.,z,1.,maxiter=200)
            if info['success']:z=zz;tt+=st;dt=min(dt*2,1.);xc2=xi
            else:
                dt=st/2
                if dt<1e-5:raise RuntimeError('param move stalled')
        s.x=x.copy();s.z=z;return build_dev(s.amp,s.ld,s.lb,s.zeta0,x),z
    def metrics(s,x):
        d,z0=s.at(x);o=d.evaluate(z0,0.);Jsc=-float(np.mean(o['Jn']+o['Jp']))
        states={0.:z0};z=z0;V=0.;br=None
        while V<1.3:
            V1=round(V+.05,6);z,_=advance(d,z,V,V1,1.,1.);states[V1]=z.copy();o=d.evaluate(z,V1);J=float(np.mean(o['Jn']+o['Jp']))
            if J>0:br=(V,V1);break
            V=V1
        cache={}
        def Jt(v):
            k=round(v,10)
            if k in cache:return cache[k]
            lo=max(q for q in states if q<=v+1e-12);zz,_=advance(d,states[lo],lo,v,1.,1.);o=d.evaluate(zz,v);cache[k]=float(np.mean(o['Jn']+o['Jp']));return cache[k]
        vc=brentq(Jt,*br,xtol=2e-8)
        r=minimize_scalar(lambda v:v*Jt(float(v)),bounds=(.55*vc,vc),method='bounded',options={'xatol':1e-6})
        FF=-r.fun/(Jsc*vc);s.n+=1
        return np.array([Jsc*1e3,vc,FF])
    def resid(s,x):
        m=s.metrics(x);r=(m-TGT)/SCALE
        pr=np.array([x[3]/.35,x[1]/.02,x[0]/.10])      # weak priors: mobility within SCLC reading tolerance, Eg, G
        print(f'  eval {s.n:3d} x={np.round(x,5)} Jsc={m[0]:.4f} Voc={m[1]:.5f} FF={m[2]:.4f} |r_data|={np.linalg.norm(r):.2f}',flush=True);return np.r_[r,pr]

if __name__=='__main__':
    amp,ld,lb=sys.argv[1],float(sys.argv[2]),float(sys.argv[3]);tag=f'{"V2" if lb else "V1"}_{amp}'
    f=Fit(amp,ld,lb);t0=time.time()
    m0=f.metrics(np.zeros(4));print(tag,'before refit: Jsc %.4f Voc %.5f FF %.4f'%tuple(m0),flush=True)
    sol=least_squares(f.resid,np.array([0.,0.,float(sys.argv[4]) if len(sys.argv)>4 else 0.,0.]),method='trf',bounds=([-.3,-.06,-6.,-.7],[.3,.06,3.,.7]),x_scale=[.05,.005,.5,.3],diff_step=1e-3,xtol=1e-7,ftol=1e-8,max_nfev=25)
    m=f.metrics(sol.x);p0=parameters(300.)[0]
    out=dict(tag=tag,before=list(m0),after=list(m),target=list(TGT),x=list(sol.x),G_scale=float(np.exp(sol.x[0])),dEg_V=float(sol.x[1]),
             S=float(np.exp(sol.x[2])),M=float(np.exp(sol.x[3])),beta300_new=(None if lb else float(p0.beta_cm3s*np.exp(sol.x[2]))),zeta0=f.zeta0,
             zeta_new=(None if not lb else float(f.zeta0*np.exp(sol.x[2]))),nfev=f.n,seconds=time.time()-t0)
    json.dump(out,open(f'results/recal_{tag}.json','w'),indent=1);print('DONE',json.dumps(out),flush=True)
