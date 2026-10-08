from pathlib import Path
import sys,json,time,numpy as np
R=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(R.parent/'gaussian_transport_20261003/code'))
import conditional_device as gd
gd.BASE=dict(gd.BASE,contact_barrier=.4)
rows=[]
for mode,V in [('gaussian',-5.),('gaussian',-10.),('gaussian',-15.),('gaussian_mu300',-15.)]:
 d=gd.GaussianDevice(641,420,mode);a=np.load(R/'data'/f'{mode}_b0.4_T420_N321_V{V:g}.npz')['z']
 z=np.column_stack([np.interp(d.x,np.linspace(0,1,len(a)),a[:,i]) for i in range(3)])
 t=time.monotonic();z,info=gd.solve_stable(d,V,z);assert not info.get('status'),info
 p,o=d.ledger(z,V);p.update(iterations=info['iterations'],scaled_residual=info['scaled_residual'],barrier_eV=.4,vbi_V=d.vbi,seconds=time.monotonic()-t)
 np.savez_compressed(R/'data'/f'{mode}_b0.4_T420_N641_V{V:g}.npz',z=z,**o)
 rows.append(p);print(json.dumps(p),flush=True)
 (R/'data'/f'refine641.json').write_text(json.dumps({'rows':rows},indent=2))
