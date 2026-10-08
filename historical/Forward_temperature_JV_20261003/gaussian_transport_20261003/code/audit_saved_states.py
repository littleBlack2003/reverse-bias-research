import json,re,numpy as np
from gaussian_device import *
rows=[]
for path in sorted((ROOT/'data').glob('*_T*_N*_V*.npz')):
 m=re.fullmatch(r'(legacy|mobility_only|gaussian|gaussian_constmu)_T([0-9.]+)_N([0-9]+)_V(-?[0-9.]+).npz',path.name)
 if not m:continue
 mode,t,n,v=m.groups();d=GaussianDevice(int(n),float(t),mode);z=np.load(path)['z'];p,o=d.ledger(z,float(v));p['file']=path.name
 if not p['gate_passed']:
  z,info=solve_stable(d,float(v),z)
  assert not info.get('status'),info
  p,o=d.ledger(z,float(v));p['file']=path.name;p['refreshed_after_final_numerical_guards']=True
  np.savez_compressed(path,z=z,**o)
 assert p['gate_passed'],p
 assert p['max_LDB_cycle_error']<1e-10,p
 assert p['charge_relative_error']<1e-8,p
 assert p['energy_relative_error']<1e-8,p
 assert min(p['transport_dissipation_Wcm2'],p['reaction_dissipation_Wcm2'],p['bim_dissipation_Wcm2'])>=-1e-25,p
 rows.append(p)
(ROOT/'data/final_source_state_audit.json').write_text(json.dumps(rows,indent=2))
print(len(rows),'saved-state ledger audits pass; this does NOT establish spatial convergence')
