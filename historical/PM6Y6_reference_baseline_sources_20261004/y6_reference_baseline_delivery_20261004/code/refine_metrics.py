from pathlib import Path
import json,glob,sys
import numpy as np
from scipy.optimize import brentq,minimize_scalar
from device import *
ROOT=Path(__file__).resolve().parents[1]
def main(name):
 path=Path(name);data=json.loads(path.read_text());rows=sorted(data['curves'].get('light',[]),key=lambda x:x['V'])
 if not rows:return
 store=np.load(path.with_suffix('.npz'));mask=store['labels'][:,0]>0;vol=store['labels'][mask,1];states=store['z'][mask];idx=np.argsort(vol);vol=vol[idx];states=states[idx]
 d=Device(Parameters(**dict(data['parameters'],allow_mobility_extrapolation=data['parameters']['T']<223)),int(rows[0]['nodes']));I=data['parameters'].get('light',rows[0]['light']);evaluated=[]
 def evaluate(V):
  k=int(np.argmin(abs(vol-V)));z,info=advance(d,states[k],float(vol[k]),float(V),I,I)
  r=d.ledger(z,V,I);evaluated.append(r);return r['J_Acm2']
 current=np.array([r['J_Acm2'] for r in rows]);V=np.array([r['V'] for r in rows]);cross=np.flatnonzero(current[:-1]*current[1:]<0)
 out=dict(tag=data['tag'],T=data['parameters']['T'],I=I,contact=data['parameters']['minority_contact'],Jsc_Acm2=None,Voc_V=None,FF=None,Vmp_V=None,Jmp_Acm2=None,voltage_resolution_V=1e-5,source='independently converged N321 PDE states',status='conditional prediction, not an experimental fit')
 if V[0]==0:out['Jsc_Acm2']=float(current[0])
 if len(cross):
  k=int(cross[0]);vc=brentq(evaluate,V[k],V[k+1],xtol=1e-6)
  m=minimize_scalar(lambda x:x*evaluate(x),bounds=(0,vc),method='bounded',options={'xatol':1e-5})
  out.update(Voc_V=vc,Vmp_V=float(m.x),Jmp_Acm2=float(m.fun/m.x),FF=float(-m.fun/(vc*abs(out['Jsc_Acm2']))) if out['Jsc_Acm2'] else None)
 else:out['unresolved']='No signed zero crossing within saved voltage window'
 out['refined_rows']=evaluated
 (ROOT/'data'/f'{data["tag"]}_metrics.json').write_text(json.dumps(out,indent=2))
 print({k:v for k,v in out.items() if k!='refined_rows'},flush=True)
if __name__=='__main__':main(sys.argv[1])
