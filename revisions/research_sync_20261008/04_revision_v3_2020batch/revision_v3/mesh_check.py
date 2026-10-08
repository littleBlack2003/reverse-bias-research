"""mesh_check.py TAG T... : same fitted parameters (per-T 3-param inversion) on N=161/321/641 nodes."""
import sys,json,time,dataclasses;sys.path.insert(0,'.')
from lib_v2 import *
tag=sys.argv[1];Ts=[int(t) for t in sys.argv[2:]];zo=orig_states([280,250,230,200,175,150,125,100]);out={}
for T in Ts:
    o=json.load(open(f'results/fitT3_{tag}_{T}.json'));x=np.array(o['x']);d0,z0=to_final(tag,T,zo[T])
    def mk(x,N=321):
        p=dataclasses.replace(d0.p,beta_cm3s=d0.p.beta_cm3s*np.exp(x[0]),G_cm3s=d0.p.G_cm3s*np.exp(x[2]))
        return ClosureDevice(p,N,amp_n=d0.amp_n*np.exp(x[1]),amp_p=d0.amp_p*np.exp(x[1]),lam_d=d0.lam_d,lam_b=d0.lam_b,zeta=None)
    z=z0;xc=np.zeros(3);tt=0.;dt=.5
    while tt<1-1e-12:
        s=min(dt,1-tt);d=mk(xc+(tt+s)*(x-xc));zz,info=solve(d,0.,z,1.,maxiter=250)
        if info['success']:z=zz;tt+=s;dt=min(dt*2,1.)
        else:dt=s/2
    out[T]={}
    for N in [161,321,641]:
        dN=mk(x,N)
        if N==321:zN=z
        else:
            zN=np.column_stack([np.interp(dN.x,d.x,z[:,k]) for k in range(z.shape[1])])
            zz,info=solve(dN,0.,zN,1.,maxiter=300)
            if not info['success']:print(T,N,'direct solve failed',info.get('reason'),flush=True);continue
            zN=zz
        vc,FF,Jsc=jv(dN,zN,dv=.02);out[T][N]=(Jsc,vc,FF);print(T,N,'Jsc %.4f Voc %.5f FF %.3f%%'%(Jsc,vc,100*FF),flush=True)
json.dump(out,open(f'results/mesh_{tag}.json','w'),indent=1)
