"""Sign-resolved independent N641 Voc brackets; interpolation is not certification."""
from pathlib import Path
import sys,json,time,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'code'))
from device import *
T=float(sys.argv[1]);base=ROOT/'independent_review'/f'mesh_refine_T{T:g}'
info=json.loads(base.with_suffix('.json').read_text());a=np.load(base.with_suffix('.npz'));info['parameters']['allow_mobility_extrapolation']=True;p=Parameters(**info['parameters']);d=Device(p,641)
vs=list(a['V']);zs=list(a['z']);rows=[]
def at(V):
 i=int(np.argmin(abs(np.array(vs)-V)));z,_=advance(d,zs[i],vs[i],V,1.,1.)
 row=d.ledger(z,V,1.);rows.append(row);vs.append(V);zs.append(z)
 return row,z
if T==100:lo,hi=1.1,1.125
else:lo,hi=info['rows'][-2]['V'],info['rows'][-1]['V']
left,zl=at(lo);right,zr=at(hi)
if not left['J_Acm2']+left['current_closure_uncertainty_Acm2']<0 or not right['J_Acm2']-right['current_closure_uncertainty_Acm2']>0:
 raise RuntimeError(('no resolved Voc bracket',T,lo,left['J_Acm2'],hi,right['J_Acm2']))
while hi-lo>1e-5:
 mid=(lo+hi)/2;row,zm=at(mid)
 if row['J_Acm2']+row['current_closure_uncertainty_Acm2']<0:lo=mid;left=row;zl=zm
 elif row['J_Acm2']-row['current_closure_uncertainty_Acm2']>0:hi=mid;right=row;zr=zm
 else:
  # A numerically unresolved sign is an interval, never silently a zero.
  if hi-lo<2e-5:break
  for candidate in [mid-(hi-lo)/8,mid+(hi-lo)/8]:
   cr,cz=at(candidate)
   if cr['J_Acm2']+cr['current_closure_uncertainty_Acm2']<0 and candidate>lo:lo=candidate;left=cr;zl=cz
   elif cr['J_Acm2']-cr['current_closure_uncertainty_Acm2']>0 and candidate<hi:hi=candidate;right=cr;zr=cz
slope=(right['J_Acm2']-left['J_Acm2'])/(hi-lo)
root=lo-left['J_Acm2']/slope
res={'T':T,'nodes':641,'Voc_bracket_V':[lo,hi],'Voc_interpolated_inside_bracket_V':root,'bracket_width_V':hi-lo,'slope_Acm2_per_V':slope,'max_closure_uncertainty_as_V':max(left['current_closure_uncertainty_Acm2'],right['current_closure_uncertainty_Acm2'])/slope,'endpoint_signs_resolved':True,'rows':rows,'parameters':info['parameters'],'source_sha256':hashlib.sha256((ROOT/'code'/'device.py').read_bytes()).hexdigest(),'scope':'Numerical bracket for conditional ideal-contact extrapolative device; excludes model/calibration and mesh uncertainty'}
(ROOT/'independent_review'/f'voc_refine_T{T:g}.json').write_text(json.dumps(res,indent=2));np.savez_compressed(ROOT/'independent_review'/f'voc_refine_T{T:g}.npz',z=np.array([zl,zr]),V=np.array([lo,hi]),x_cm=d.x*p.d_cm)
print('VOC',T,root,[lo,hi],'closureV',res['max_closure_uncertainty_as_V'],flush=True)
