"""Re-evaluate saved relative-quasi-Fermi states under final closure gates."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'code'))
from device import Device,Parameters
rows=[]
for f in sorted((ROOT/'data').glob('T*_N*_I*_*.json')):
 try:
  run=json.loads(f.read_text());params=dict(run['parameters']);params['allow_mobility_extrapolation']=True;p=Parameters(**params);a=np.load(f.with_suffix('.npz'))
 except (OSError,ValueError,KeyError): continue
 d=Device(p,nodes=a['z'].shape[1])
 for z,(L,V) in zip(a['z'],a['labels']):
  info=d.ledger(z,float(V),float(L));r=d.residual(z,float(V),float(L)).reshape(z.shape)
  o=d.evaluate(z,float(V))
  info.update(tag=run['tag'],source_file=f.name,local_continuity_L1_Acm2=float(np.sum(abs(r[1:-1,1:]))),minimum_transport_sign=float(min(np.min(o['Jn']*np.diff(z[:,1])),-np.max(o['Jp']*np.diff(z[:,2])))),minimum_reaction_sign=float(np.min(o['R']*o['A'])))
  rows.append(info)
out={'state_coordinates':'contact-relative quasi Fermi; absolute-QF archived attempts excluded','sha256':{s:hashlib.sha256((ROOT/'code'/s).read_bytes()).hexdigest() for s in ['device.py','fd_eos.py']},'rows':rows,'count':len(rows),'gate_failed':sum(not x['gate_passed'] for x in rows),'unresolved_sign_count':sum(not x['numerical_sign_resolved'] for x in rows if x['V'] or x['light'])}
(ROOT/'independent_review'/'saved_state_audit.json').write_text(json.dumps(out,indent=2))
print('States',len(rows),'failed current/charge gate',out['gate_failed'],'unresolved sign',out['unresolved_sign_count'])
for tag in sorted({x['tag'] for x in rows}):
 rr=[x for x in rows if x['tag']==tag];lit=[x for x in rr if x['light']];dark=[x for x in rr if not x['light'] and x['V']]
 print(tag,'N',len(rr),'light N',len(lit),'max current-relative light spread',max([x['current_spread_relative_to_terminal'] for x in lit],default=None),'min dark J',min([x['J_Acm2'] for x in dark],default=None),'failed',sum(not x['gate_passed'] for x in rr))
