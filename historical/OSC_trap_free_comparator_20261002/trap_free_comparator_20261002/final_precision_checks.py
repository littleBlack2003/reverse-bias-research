from pathlib import Path
import json
import comparator_model as m
r=[]
for F in [.5,3.]:
 for topology in ['direct','trap']:
  a=m.solve(F,80,topology)
  for order,dps in [(48,220),(24,320)]:
   b=m.solve(F,80,topology,order=order,dps=dps)
   r.append(dict(F_MVcm=F,T_K=80,topology=topology,order=order,dps=dps,relative_rate_error=abs(b['R_pair_s']/a['R_pair_s']-1)))
assert max(x['relative_rate_error'] for x in r)<1e-8
p=Path(__file__).parent/'review/final_precision_checks.json';p.write_text(json.dumps(r,indent=2));print('maximum relative rate error',max(x['relative_rate_error'] for x in r))
