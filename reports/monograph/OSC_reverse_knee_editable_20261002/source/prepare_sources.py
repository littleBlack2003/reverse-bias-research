from pathlib import Path
import shutil,json,hashlib
P=Path(__file__).resolve().parents[1]; W=Path('/workspace/shared/wxh'); M=W/'002高YZ/OSC_FF_reverse_knee_model/model_comparison'; R=W/'research'
items={
 'overview':M/'results_ff78/comparison_overview.png',
 'local_pf':M/'results_ff78/local_classical_pf/JV_comparison.png',
 'abcd':M/'results_ff78/four_rate_diagnostic/JV_ABCD.png',
 'abcd_profiles':M/'results_ff78/four_rate_diagnostic/profiles_ABCD.png',
 'density':M/'results_ff78/pf_density_scan/density_FF_V50.png',
 'thickness':M/'results_ff78/pf_thickness_scan/thickness_K_FF.png',
 'thickness_curves':M/'results_ff78/pf_thickness_scan/thickness_JV.png',
 'tat':M/'results_ff78/tat_results/tat_comparison.png',
 'thermal':M/'results_ff78/thermal_rate_dependence.png',
 'branch_mesh':R/'revalidated_rate_audit_20261002/figures/branch_and_mesh.png',
 'spatial_balance':R/'revalidated_rate_audit_20261002/figures/spatial_balance.png',
 'newton_failure':R/'asymmetric_solver_diagnostics_20261002/figures/newton_failure.png',
 'qf_recovery':R/'asymmetric_solver_diagnostics_20261002/figures/recovered_branch_and_gate.png',
 'quantum_kernel':R/'temperature_mechanism_audit_20261002/figures/quantum_and_balance.png',
 'pair_capacity':R/'temperature_mechanism_audit_20261002/figures/pair_capacity.png',
 'dos_extract':R/'temperature_mechanism_audit_20261002/figures/dos_and_extraction.png',
 'graph_field':R/'quantum_cycle_audit_20261002/figures/cycle_vs_field.png',
 'graph_temperature':R/'quantum_cycle_audit_20261002/figures/temperature_and_quantum_ratio.png',
 'graph_blocking':R/'quantum_cycle_audit_20261002/figures/extraction_blocking.png',
 'graph_feedback':R/'quantum_cycle_audit_20261002/figures/coulomb_extraction_feedback.png',
 'graph_wrong':R/'quantum_cycle_audit_20261002/figures/failed_naive_reverse.png',
}
manifest=[]
for name,src in items.items():
 if src.exists():
  dst=P/'figures'/f'{name}.png';shutil.copy2(src,dst)
  manifest.append({'id':name,'original':str(src),'copy':str(dst.relative_to(P)),'sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'kind':'existing_computation_figure','recomputed_here':False})
 else:print('MISSING',src)
# Freeze complete short reports and machine-readable support, not huge caches.
for directory in [M/'results_ff78',R/'revalidated_rate_audit_20261002',R/'asymmetric_solver_diagnostics_20261002',R/'temperature_mechanism_audit_20261002',R/'quantum_cycle_audit_20261002',R/'d18_l8bo_constraints_20261002']:
 for src in directory.rglob('*'):
  if src.is_file() and src.suffix.lower() in ['.md','.csv','.json'] and src.stat().st_size<10_000_000 and not any(s in src.parts for s in ['source_snapshot','__pycache__','cache']):
   rel=src.relative_to(W);dst=P/'data'/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
   manifest.append({'original':str(src),'copy':str(dst.relative_to(P)),'sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'kind':'research_source'})
(P/'data/source_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
print(len(items),'figure selections;',len(manifest),'provenance entries')
