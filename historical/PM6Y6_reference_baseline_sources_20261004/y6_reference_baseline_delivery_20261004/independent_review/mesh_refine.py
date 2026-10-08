"""Selected refined-mesh solves, owned by independent numerical review."""
from pathlib import Path
from dataclasses import asdict
import sys,json,time,traceback
import numpy as np
from scipy.interpolate import PchipInterpolator
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'code'))
from device import *
from calibration import beta,G0,bulk_voc
T=float(sys.argv[1]);N=641
p=Parameters(T=T,G_cm3s=G0,beta_cm3s=float(beta(T)),allow_mobility_extrapolation=True)
d=Device(p,N);rows=[];states=[];traces=[];failed=[]
source=ROOT/'data'/f'T{T:g}_N321_I1_blocking_bn0_bp0.npz'
if source.exists():
 a=np.load(source);light_indices=np.flatnonzero(a['labels'][:,0]==1)
else:light_indices=[]
# Pick a sign-bracketing neighborhood from coarse J-V, if available.
vbulk=float(bulk_voc(T));vnear=vbulk
jf=source.with_suffix('.json')
if jf.exists():
 run=json.loads(jf.read_text());light=sorted(run.get('curves',{}).get('light',[]),key=lambda r:r['V'])
 for l,r in zip(light[:-1],light[1:]):
  if l['J_Acm2']*r['J_Acm2']<=0:
   vnear=l['V']-l['J_Acm2']*(r['V']-l['V'])/(r['J_Acm2']-l['J_Acm2']);break
# Exact same voltages allow direct comparison; near-root pair additionally checks sign.
targets=[0.,.8,round(vnear-.0125,8),round(vnear+.0125,8)]
for target in targets:
 start=time.time()
 try:
  if len(light_indices):
   idx=min(light_indices,key=lambda i:abs(a['labels'][i,1]-target));v0=float(a['labels'][idx,1])
   seed=PchipInterpolator(a['x_cm']/p.d_cm,a['z'][idx],axis=0)(d.x)
   seed,_info=solve(d,v0,seed,1.,maxiter=200)
   if not _info['success']:raise RuntimeError('refinement seed '+_info['reason'])
  else:
   v0=vbulk;seed=poisson_seed(d,v0);seed,_info=solve(d,v0,seed,1.,maxiter=240)
   if not _info['success']:raise RuntimeError('Poisson seed '+_info['reason'])
  z,info=advance(d,seed,v0,target,1.,1.,traces)
  info=d.ledger(z,target,1.);info['wall_seconds']=time.time()-start
  rows.append(info);states.append(z)
  print('PASS',T,target,info['J_Acm2'],'relative spread',info['current_spread_relative_to_terminal'],'closure bound',info['current_closure_uncertainty_Acm2'],flush=True)
 except Exception as e:
  failed.append({'V':target,'error':str(e),'traceback':traceback.format_exc()});print('FAIL',T,target,str(e),flush=True)
 result={'T':T,'nodes':N,'parameters':asdict(p),'rows':rows,'failure':failed,'traces':traces,'state_coordinates':'contact-relative quasi-Fermi','purpose':'selected 321-to-641 spatial refinement; no new parameter calibration'}
 (ROOT/'independent_review'/f'mesh_refine_T{T:g}.json').write_text(json.dumps(result,indent=2))
 if states:np.savez_compressed(ROOT/'independent_review'/f'mesh_refine_T{T:g}.npz',z=np.array(states),V=np.array([x['V'] for x in rows]),x_cm=d.x*p.d_cm)
