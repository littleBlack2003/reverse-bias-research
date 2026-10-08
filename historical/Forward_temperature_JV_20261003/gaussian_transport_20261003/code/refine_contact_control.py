import json,numpy as np
from gaussian_device import *
T=420.;V=-5.;mode='gaussian_constmu'
saved=np.load(ROOT/'data/gaussian_constmu_T420_N81_V-5.npz')['z'];sx=np.linspace(0,1,len(saved));rows=[]
for N in [161,321,641]:
 d=GaussianDevice(N,T,mode);z=np.column_stack([np.interp(d.x,sx,saved[:,i]) for i in range(3)])
 z,info=solve_stable(d,V,z)
 if info.get('status'):raise RuntimeError((N,info.get('status'),info['scaled_residual']))
 p,o=d.ledger(z,V);rows.append(p);saved=z;sx=d.x
 np.savez_compressed(ROOT/'data'/f'contact_control_T420_N{N}.npz',z=z,**o)
 (ROOT/'data/contact_control.json').write_text(json.dumps(rows,indent=2))
 print(N,p['J_Acm2'],p['relative_error'],flush=True)
