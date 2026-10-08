import json,sys,time
from gaussian_device import *
N=int(sys.argv[1]);T=float(sys.argv[2]);mode=sys.argv[3];targets=[float(s) for s in sys.argv[4:]] or [0,-5,-10,-15]
d=GaussianDevice(N,T,mode);z=np.c_[d.x*d.vbi/d.vt,np.zeros(N),np.zeros(N)]
z,info=solve_stable(d,0,z);assert not info.get('status'),info
rows=[];trace=[];old=0.;start=time.monotonic()
for target in targets:
 for V in np.linspace(old,target,max(1,int(np.ceil(abs(target-old)/.5)))+1)[1:]:
  z,info=advance(d,z,old,float(V),trace);old=float(V)
 p,o=d.ledger(z,target);p['iterations']=info['iterations'];p['scaled_residual']=info['scaled_residual'];rows.append(p)
 tag=f'{mode}_T{T:g}_N{N}_V{target:g}'
 np.savez_compressed(ROOT/'data'/f'{tag}.npz',z=z,**o)
 print(json.dumps(p),flush=True)
 (ROOT/'data'/f'{mode}_T{T:g}_N{N}.json').write_text(json.dumps(dict(rows=rows,trace=trace,seconds=time.monotonic()-start),indent=2))
