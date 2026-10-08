from pathlib import Path
import json,csv,hashlib,sys
import numpy as np
import spatial_model as sm
from run_audit import writecsv
ROOT=Path(__file__).parent
rows=[]
for T in [80,300]:
 for F in np.round(np.arange(1.4,2.2001,.05),8):
  for cap in [False,True]:rows.append(sm.solve(float(F),T,single_pair=cap))
writecsv('multiple_occupancy_controls.csv',rows)
kr=[]
for T in [80,300]:
 for F in [1.5,2.,3.]:
  for kind in ['classical','quantum']:kr.append(sm.solve(F,T,kind=kind))
writecsv('classical_quantum_controls.csv',kr)
eq=[]
for T in [80,300]:
 for F in [0.,1.5,3.]:eq.append(sm.solve(F,T,delta_mu=0.,mean_mu=.023,single_pair=True))
writecsv('restricted_graph_equilibrium.csv',eq)
# Preserve representative high-field full ledgers, where the simple pair bottleneck fails.
er=[];sr=[]
for F in [2.,3.]:
 for cap in [False,True]:
  r,p,Q,rt,de,E,mu,g,states,edges=sm.solve(F,80,single_pair=cap,detail=True)
  for z in de:er.append(dict(T_K=80,F_MVcm=F,single_pair_constraint=cap,**z))
  for st,pp,ee in zip(states,p,E):sr.append(dict(T_K=80,F_MVcm=F,single_pair_constraint=cap,state=''.join(map(str,st)),probability=float(pp),energy_eV=float(ee)))
writecsv('high_field_edges.csv',er);writecsv('high_field_states.csv',sr)
# Bare instantaneous formation at the fully neutral reference; not a steady source.
bare=[]
for F in [.5,1.,1.5,2.,3.]:
 for T in [80,300]:
  r,p,Q,rt,de,E,mu,g,states,edges=sm.solve(F,T,detail=True)
  neutral=''.join(map(str,g['b']))
  for z in de:
   if z['source']==neutral and z['process']=='V_D':bare.append(dict(T_K=T,F_MVcm=F,neutral_state=neutral,bare_first_transfer_s=z['forward_rate_s'],stationary_cycle_s=r['R_s'],ratio_bare_over_stationary=z['forward_rate_s']/r['R_s']))
writecsv('instantaneous_vs_steady.csv',bare)
validation=dict(points=len(rows)+len(kr)+len(eq),max_node_continuity_relative=max(r['max_node_continuity_relative'] for r in rows+kr+eq),max_heat_chemical_relative=max(r['heat_chemical_relative'] for r in rows+kr+eq),max_common_mu_probability_relative=max(r['common_mu_probability_relative'] or 0 for r in eq),min_probability=min(r['min_probability'] for r in rows+kr+eq))
(ROOT/'results'/'followup_validation.json').write_text(json.dumps(validation,indent=2));print(validation)

(ROOT/'results'/'followup_provenance.json').write_text(json.dumps(dict(python=sys.version,code_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ['spatial_model.py','run_followup_controls.py','run_audit.py']}),indent=2))
