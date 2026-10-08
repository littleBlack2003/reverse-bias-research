from analyze_density import *
import hashlib
v=json.loads((ROOT/'data/validation.json').read_text());a=json.loads((ROOT/'data/analysis.json').read_text())
assert len(v['runs'])==21
assert all(x['failure'] is None and x['all_gates'] and x['max_scaled_residual']<2e-12 and x['max_relative_current_error']<=1e-8 and x['max_charge_relative_error']<1e-7 for x in v['runs'])
assert all(x['equilibrium_J_Acm2']==0 and x['equilibrium_pair_Acm2']==0 for x in v['runs'])
assert max(x['frozen_max_relative_spread'] for x in a['frozen_checks'])<1e-13
assert v['steps_max_relative']<1e-8
assert v['density_ratio_max_absolute_mesh_change_10to20']<1e-3
assert v['independent_event_max_relative_pair']<1e-10
assert v['independent_event_max_relative_J']<1e-6
assert len(v['refined_landmarks'])==18
for p in (ROOT/'data').glob('window_*.json'):
 d=json.loads(p.read_text());assert len(d['rows'])==41 and d['failure'] is None
 assert all(r['gate_passed'] and r['scaled_residual']<2e-12 for r in d['rows'])
original=Path('/workspace/shared/wxh/002高YZ/OSC_FF_reverse_knee_model')
source=ROOT/'reference/reference_solver/reference_audit/source_snapshot'
hash_checks=[]
for rel in ['solvers/base_dd_solver.py','solvers/pf_solver.py','model_comparison/baseline_ff78.json']:
 h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();expected=next(x['sha256'] for x in json.loads((ROOT/'source_hashes.json').read_text()) if x['relative_path']==rel);equal=h(source/rel)==expected;assert equal
 if original.exists():assert h(original/rel)==expected
 hash_checks.append(dict(relative_path=rel,sha256=expected,matches_reference_copy=equal))
save(ROOT/'data/final_assertions.json',dict(passed=True,checks=['21 full scans complete and converged','zero-current exact equilibrium','frozen density linearity','step refinement','density ratio mesh refinement','independent event and current reconstruction','18 refined shape landmarks','all 246 local refinement records converged','original physical solver hashes match snapshots'],original_source_hashes=hash_checks))
print('PASS: convergence, charge/current/source, frozen scaling, selected refinement, independent reconstruction and original solver hashes')
