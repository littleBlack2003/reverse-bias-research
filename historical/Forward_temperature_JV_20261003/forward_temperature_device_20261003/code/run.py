from pathlib import Path
import sys,json,time,os
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent/'gaussian_transport_20261003/code'))
sys.path.insert(0,str(ROOT.parent/'complete_temperature_device_20261003/code'))
import conditional_device as gd
import numpy as np
N=int(sys.argv[1]);T=float(sys.argv[2]);barrier=float(sys.argv[3]);gd.BASE=dict(gd.BASE,contact_barrier=barrier)
d=gd.GaussianDevice(N,T,'gaussian');z=np.c_[d.x*d.vbi/d.vt,np.zeros(N),np.zeros(N)]
z,info=gd.solve_stable(d,0,z);assert not info.get('status'),info
rows=[];trace=[];old=0.;start=time.monotonic();tag=f'b{barrier:g}_T{T:g}_N{N}'
for target in np.round(np.arange(0,1.20001,.05),8):
 try:
  if target: z,info=gd.advance(d,z,old,float(target),trace)
  old=float(target);p,o=d.ledger(z,float(target));p.update(iterations=info['iterations'],scaled_residual=info['scaled_residual'],barrier_eV=barrier,vbi_V=d.vbi)
  rows.append(p);np.savez_compressed(ROOT/'data'/f'{tag}_V{target:g}.npz',z=z,**o)
  print(target,p['J_Acm2'],p['max_density_fraction'],p['max_reduced_field'],flush=True)
  (ROOT/'data'/f'{tag}.json').write_text(json.dumps(dict(rows=rows,trace=trace,seconds=time.monotonic()-start),indent=2))
  if p['max_density_fraction']>.01 or p['max_reduced_field']>3:break
 except Exception as e:
  (ROOT/'data'/f'{tag}_failure.txt').write_text(repr(e));raise
