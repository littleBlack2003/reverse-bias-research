from pathlib import Path
import sys,json
import numpy as np
from scipy.interpolate import PchipInterpolator
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'code'))
from device import *
f=ROOT/'independent_review/mesh_refine_T100';j=json.loads(f.with_suffix('.json').read_text());a=np.load(f.with_suffix('.npz'));p=Parameters(**j['parameters']);d=Device(p,1281)
i=int(np.argmin(abs(a['V'])));z=PchipInterpolator(a['x_cm']/p.d_cm,a['z'][i],axis=0)(d.x)
z,info=solve(d,0.,z,1.,maxiter=200)
if not info['success']:raise RuntimeError(info['reason'])
row=d.ledger(z,0.,1.);prev=next(r['J_Acm2'] for r in j['rows'] if r['V']==0)
row.update(Jsc641_Acm2=prev,relative_641_to_1281_change=abs(row['J_Acm2']/prev-1),scope='One selected100K short-circuit spatial refinement; conditional extrapolative material model')
(ROOT/'independent_review/jsc_100K_N1281.json').write_text(json.dumps(row,indent=2));np.savez_compressed(ROOT/'independent_review/jsc_100K_N1281.npz',z=z,V=0.,x_cm=d.x*p.d_cm)
print(json.dumps(row,indent=2))
