from pathlib import Path
import json,hashlib,zipfile,shutil
ROOT=Path(__file__).resolve().parents[1]
needed=[
 'reference/code/stable_cycle_solver.py','reference/code/temperature_probe.py',
 'reference/reference_solver/code/qf_solver.py','reference/reference_solver/code/solver_probe.py',
 'reference/reference_solver/reference_audit/code/rate_audit.py',
 'reference/reference_solver/reference_audit/source_snapshot/solvers/__init__.py',
 'reference/reference_solver/reference_audit/source_snapshot/solvers/base_dd_solver.py',
 'reference/reference_solver/reference_audit/source_snapshot/solvers/pf_solver.py',
 'reference/reference_solver/reference_audit/source_snapshot/model_comparison/baseline_ff78.json']
paths=sorted([p for folder in ['code','data','figures','logs'] for p in (ROOT/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts]+[p for p in ROOT.iterdir() if p.is_file() and p.name!='MANIFEST_SHA256.json']+[ROOT/p for p in needed])
manifest=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths]
(ROOT/'MANIFEST_SHA256.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False),encoding='utf-8');paths.append(ROOT/'MANIFEST_SHA256.json')
out=Path('/workspace/scratch/d9cc890e37cd/density_knee_delivery');out.mkdir(exist_ok=True)
zpath=out/'density_knee_disentangling_300K_20261002.zip'
with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in paths:z.write(p,Path(ROOT.name)/p.relative_to(ROOT))
with zipfile.ZipFile(zpath) as z:
 assert z.testzip() is None
 for e in manifest:assert hashlib.sha256(z.read(ROOT.name+'/'+e['path'])).hexdigest()==e['sha256']
for p in [ROOT/'REPORT.md',*sorted((ROOT/'figures').glob('*.png'))]:shutil.copy2(p,out/p.name)
print(json.dumps(dict(archive=str(zpath),bytes=zpath.stat().st_size,sha256=hashlib.sha256(zpath.read_bytes()).hexdigest(),included_files=len(paths),all_manifest_hashes_verified=True),indent=2))
