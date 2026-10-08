from pathlib import Path
import json,csv,platform,numpy,scipy
ROOT=Path(__file__).resolve().parents[1]
rows=[];lookup={};grids=[]
for d in [0,1]:
 for n in [81,161,321]:
  a=json.loads((ROOT/'data'/f'd{d}_N{n}.json').read_text());rows+=a['rows'];lookup[d,n]=a['rows']
 for a,b in zip(lookup[d,161][1:],lookup[d,321][1:]):grids.append(dict(distance_nm=d,V=a['V'],relative_change=abs(a['J_Acm2']/b['J_Acm2']-1)))
for r in rows:
 assert r['gate_passed']
 assert r['charge_relative_error']<1e-8
 assert r['energy_relative_error']<1e-8
 assert r['max_LDB_cycle_error']<1e-10
 assert r['occupancy_min']>=0 and r['occupancy_max']<=1
 if r['V']==0: assert r['J_Acm2']==r['pair_source_Acm2']==0
 assert abs(r['trap_areal_cm2']/9.6e9-1)<1e-12
assert max(r['relative_change'] for r in grids)<.001
with (ROOT/'data/point_summary.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
a=dict(solved_points=len(rows),continuation_points=sum(len(json.loads((ROOT/'data'/f'd{d}_N{n}.json').read_text())['trace']) for d in [0,1] for n in [81,161,321]),equilibrium_exact_zero=True,max_scaled_residual=max(r['scaled_residual'] for r in rows),max_relative_current_ledger_error=max(r['relative_error'] for r in rows),max_relative_charge_ledger_error=max(r['charge_relative_error'] for r in rows),max_relative_energy_ledger_error=max(r['energy_relative_error'] for r in rows),max_LDB_cycle_error=max(r['max_LDB_cycle_error'] for r in rows),grid_161_to_321=grids,parameters=dict(temperature_K=300,d_nm=100,Eg_eV=1.3,Nc_cm3=1e19,Nt_cm3=1e15,t0_eV=.001,lambda_eV=.3,xi_nm=.5,distance_per_leg_nm=[0,1],trap_support_nm=[2,98],integrated_traps_cm2=9.6e9),environment=dict(python=platform.python_version(),numpy=numpy.__version__,scipy=scipy.__version__))
(ROOT/'data/validation_summary.json').write_text(json.dumps(a,indent=2));print(json.dumps(a,indent=2))
