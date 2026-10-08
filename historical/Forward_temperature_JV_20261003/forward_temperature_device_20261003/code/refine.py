from pathlib import Path
import sys,json,numpy as np
R=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(R.parent/'gaussian_transport_20261003/code'));sys.path.insert(0,str(R.parent/'complete_temperature_device_20261003/code'))
import conditional_device as gd
N=int(sys.argv[1]);T=float(sys.argv[2]);b=float(sys.argv[3]);V=float(sys.argv[4]);gd.BASE=dict(gd.BASE,contact_barrier=b)
d=gd.GaussianDevice(N,T,'gaussian')
seed=np.load(R/'data'/f'b{b:g}_T{T:g}_N161_V{V:g}.npz')['z'];z=np.array([np.interp(d.x,np.linspace(0,1,len(seed)),seed[:,i]) for i in range(3)]).T
z,info=gd.solve_stable(d,V,z);assert not info.get('status'),info
p,o=d.ledger(z,V);p.update(barrier_eV=b,scaled_residual=info['scaled_residual'],iterations=info['iterations'])
np.savez_compressed(R/'data'/f'b{b:g}_T{T:g}_N{N}_V{V:g}.npz',z=z,**o)
(R/'data'/f'refine_b{b:g}_T{T:g}_N{N}_V{V:g}.json').write_text(json.dumps(p,indent=2));print(p,flush=True)
