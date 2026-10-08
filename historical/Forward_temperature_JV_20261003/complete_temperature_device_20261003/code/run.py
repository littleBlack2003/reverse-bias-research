from pathlib import Path
import sys,json,time
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent/'gaussian_transport_20261003/code'))
import conditional_device as gd
import numpy as np
N=int(sys.argv[1]);T=float(sys.argv[2]);barrier=float(sys.argv[3]);mode=sys.argv[4] if len(sys.argv)>4 else 'gaussian'
gd.BASE=dict(gd.BASE,contact_barrier=barrier)
d=gd.GaussianDevice(N,T,mode)
z=np.c_[d.x*d.vbi/d.vt,np.zeros(N),np.zeros(N)]
z,info=gd.solve_stable(d,0,z);assert not info.get('status'),info
rows=[];trace=[];old=0.;start=time.monotonic()
for target in [0.,-5.,-10.,-15.]:
 for V in np.linspace(old,target,max(1,int(np.ceil(abs(target-old)/.5)))+1)[1:]:
  z,info=gd.advance(d,z,old,float(V),trace);old=float(V)
 p,o=d.ledger(z,target);p.update(iterations=info['iterations'],scaled_residual=info['scaled_residual'],barrier_eV=barrier,vbi_V=d.vbi)
 rows.append(p);tag=f'{mode}_b{barrier:g}_T{T:g}_N{N}'
 np.savez_compressed(ROOT/'data'/f'{tag}_V{target:g}.npz',z=z,**o)
 print(json.dumps(p),flush=True)
 (ROOT/'data'/f'{tag}.json').write_text(json.dumps(dict(rows=rows,trace=trace,seconds=time.monotonic()-start),indent=2))
