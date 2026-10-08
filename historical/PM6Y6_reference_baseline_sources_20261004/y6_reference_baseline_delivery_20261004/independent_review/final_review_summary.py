"""Aggregate final reference-kernel checks and re-evaluate selected refined states."""
from pathlib import Path
import json,sys,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[1];R=ROOT/'independent_review';sys.path.insert(0,str(ROOT/'code'))
from device import *
summary={'scope':'Independent thermodynamic, quadrature and numerical review; not a material-calibration certificate','source_sha256':{name:hashlib.sha256((ROOT/'code'/name).read_bytes()).hexdigest() for name in ['device.py','fd_eos.py']},'unit_checks':{},'selected_refinement':[],'voc_refinement':[],'final_saved_state_checks':[]}
for f in ['fd_eos_results.json','flux_comparison_results.json','device_results.json']:
 d=json.loads((R/f).read_text());summary['unit_checks'][f]={'checks':len(d['tests']),'all_passed':all(x['error']<x['tolerance'] for x in d['tests'])}
for T in [100,200,300]:
 f=R/f'mesh_refine_T{T}';info=json.loads(f.with_suffix('.json').read_text());a=np.load(f.with_suffix('.npz'));p=Parameters(**info['parameters']);dev=Device(p,641)
 rows=[]
 for z,V in zip(a['z'],a['V']):
  row=dev.ledger(z,float(V),1.);assert row['gate_passed'];assert row['transport_dissipation_Wcm2']>=0 and row['recombination_dissipation_Wcm2']>=0
  assert row['energy_relative']<1e-7;rows.append(row)
 root=R/f'voc_refine_T{T}';vr=json.loads(root.with_suffix('.json').read_text());va=np.load(root.with_suffix('.npz'))
 left,right=[dev.ledger(z,float(V),1.) for z,V in zip(va['z'],va['V'])]
 assert left['J_Acm2']+left['current_closure_uncertainty_Acm2']<0
 assert right['J_Acm2']-right['current_closure_uncertainty_Acm2']>0
 assert left['gate_passed'] and right['gate_passed']
 summary['final_saved_state_checks'].append({'T':T,'selected_points_passed':len(rows),'root_endpoints_passed':2,'max_energy_relative':max(x['energy_relative'] for x in rows),'max_charge_relative':max(x['charge_relative'] for x in rows)})
 tag=f'T{T}_N321_I1_blocking_bn0_bp0'+('_lightonly' if T==100 else '')
 met=json.loads((ROOT/'data'/f'{tag}_metrics.json').read_text());coarse=json.loads((ROOT/'data'/f'{tag}.json').read_text())
 for V in [0.,.8]:
  hi=next(r for r in rows if r['V']==V);lo=next(r for r in coarse['curves']['light'] if r['V']==V)
  summary['selected_refinement'].append({'T':T,'V':V,'J321_Acm2':lo['J_Acm2'],'J641_Acm2':hi['J_Acm2'],'relative_mesh_change':abs(hi['J_Acm2']/lo['J_Acm2']-1),'closure_relative_bound641':hi['current_closure_uncertainty_Acm2']/abs(hi['J_Acm2'])})
 summary['voc_refinement'].append({'T':T,'Voc321_V':met['Voc_V'],'Voc641_V':vr['Voc_interpolated_inside_bracket_V'],'central_mesh_shift_V':vr['Voc_interpolated_inside_bracket_V']-met['Voc_V'],'certified_641_bracket_V':vr['Voc_bracket_V'],'closure_uncertainty_as_V':vr['max_closure_uncertainty_as_V'],'resolution_note':'Mesh shift not resolved at approximately10 microvolt root precision; central digits are not experimental accuracy'})
q=json.loads((R/'jsc_100K_N1281.json').read_text());summary['Jsc100K_1281']=q
p=Parameters(T=100,beta_cm3s=json.loads((R/'mesh_refine_T100.json').read_text())['parameters']['beta_cm3s'],G_cm3s=json.loads((R/'mesh_refine_T100.json').read_text())['parameters']['G_cm3s'],allow_mobility_extrapolation=True)
d=Device(p,1281);z=np.load(R/'jsc_100K_N1281.npz')['z'];check=d.ledger(z,0.,1.);assert check['gate_passed'] and check['energy_relative']<1e-7
summary['final_saved_state_checks'].append({'T':100,'nodes':1281,'selected_points_passed':1,'final_source_gate_passed':True})
summary['remaining_limits']=['Below223K mobility requires explicit unvalidated extrapolation opt-in','Contacts, epsilon and absolute G are scenario-dependent inputs','PIA n(T) calibrates beta/G and cannot also validate it','Published fitted Eg cannot turn the sameVoc data into a holdout','Some preliminary main-directory dark states have stale or failed final gates; use saved_state_audit and producer reconciliation before release','Electrical/chemical work balance is not full optical detailed balance','Full-device results are not experimental fits; 100K current collapse contradicts interpreting lowT GDM as a calibrated transport law']
(R/'FINAL_VALIDATION.json').write_text(json.dumps(summary,indent=2));print(json.dumps({k:v for k,v in summary.items() if k not in ['Jsc100K_1281','remaining_limits']},indent=2))
