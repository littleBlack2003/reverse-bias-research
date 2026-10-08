from pathlib import Path
import sys,json,numpy as np
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R.parent/'gaussian_transport_20261003/code'));sys.path.insert(0,str(R.parent/'complete_temperature_device_20261003/code'));import conditional_device as gd
out=[]
for b,T,V in [(.4,300,.3),(.4,420,1.2)]:
 gd.BASE=dict(gd.BASE,contact_barrier=b);d=gd.GaussianDevice(161,T,'gaussian',order=256)
 z=np.load(R/'data'/f'b{b:g}_T{T:g}_N161_V{V:g}.npz')['z'];z,info=gd.solve_stable(d,V,z);assert not info.get('status'),info
 p,o=d.ledger(z,V);p.update(barrier_eV=b,order=256);out.append(p)
(R/'data'/'quadrature256.json').write_text(json.dumps(out,indent=2));print(out)
