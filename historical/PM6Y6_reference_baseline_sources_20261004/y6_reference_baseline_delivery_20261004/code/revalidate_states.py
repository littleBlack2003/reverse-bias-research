"""Recompute every stored ledger under final source; repair warm states if needed."""
from pathlib import Path
import json,sys,hashlib
import numpy as np
from device import *
ROOT=Path(__file__).resolve().parents[1]
def main(path):
 f=Path(path);a=json.loads(f.read_text());p=dict(a['parameters']);p['allow_mobility_extrapolation']=not 223<=p['T']<=328;a['parameters']=p
 store=np.load(f.with_suffix('.npz'));states=store['z'].copy();labels=store['labels'].copy();d=Device(Parameters(**p),len(states[0]));repaired=[];failed=[]
 for k,(L,V) in enumerate(labels):
  out=d.ledger(states[k],V,L)
  if not out['gate_passed']:
   states[k],info=solve(d,V,states[k],L,maxiter=240);out=d.ledger(states[k],V,L)
   (repaired if out['gate_passed'] else failed).append(dict(index=k,V=float(V),light=float(L),info={x:y for x,y in info.items() if x!='history'}))
  key='light' if L else'dark'
  for i,r in enumerate(a['curves'][key]):
   if r['V']==V:a['curves'][key][i]=out;break
 a['state_coordinates']='contact_relative_qf_v2';a['final_source_hashes']={x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in [ROOT/'code/device.py',ROOT/'code/fd_eos.py']};a['final_revalidation']=dict(repaired=repaired,failed=failed,total_states=len(states),all_stored_states_pass=not failed)
 f.write_text(json.dumps(a,indent=2));np.savez_compressed(f.with_suffix('.npz'),z=states,labels=labels,x_cm=store['x_cm'])
 print(f.name,'states',len(states),'repaired',len(repaired),'failed',len(failed),flush=True)
if __name__=='__main__':main(sys.argv[1])
