from pathlib import Path
import json,sys
import numpy as np
from device import *
f=Path(sys.argv[1]);a=json.loads(f.read_text());s=np.load(f.with_suffix('.npz'));lab=s['labels'].copy();states=s['z'].copy();p=Parameters(**dict(a['parameters'],allow_mobility_extrapolation=a['parameters']['T']<223));d=Device(p,len(states[0]));I=1.
ids=np.flatnonzero(lab[:,0]==I);idx=ids[np.argmax(lab[ids,1])];old=float(lab[idx,1]);z=states[idx]
for V in [1.125,1.15]:
 if V<=old:continue
 z,info=advance(d,z,old,V,I,I);r=d.ledger(z,V,I);r['range_note']='Extended solely to bracket conditional model Voc above reference window';a['curves']['light'].append(r);states=np.concatenate([states,z[None,...]]);lab=np.concatenate([lab,[[I,V]]]);old=V
 print(V,r['J_Acm2'],flush=True)
a['curves']['light'].sort(key=lambda x:x['V']);f.write_text(json.dumps(a,indent=2));np.savez_compressed(f.with_suffix('.npz'),z=states,labels=lab,x_cm=s['x_cm'])
