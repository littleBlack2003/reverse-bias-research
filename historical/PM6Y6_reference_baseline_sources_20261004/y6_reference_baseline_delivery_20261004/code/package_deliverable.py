from pathlib import Path
import shutil,json,hashlib,zipfile
ROOT=Path(__file__).resolve().parents[1];STAGE=ROOT.parent/'y6_reference_baseline_delivery_20261004';STAGE.mkdir(exist_ok=True)
# Rebuild this dedicated generated staging directory only.
for p in STAGE.iterdir():
 if p.is_dir():shutil.rmtree(p)
 else:p.unlink()
for name in ['code','reference','figures','independent_review']:
 shutil.copytree(ROOT/name,STAGE/name,ignore=shutil.ignore_patterns('__pycache__','*.log','saved_state_audit.json','quadrature_probe.py'))
for name in ['README.md','requirements.lock.txt','runtime.json','reproduce_smoke.sh','reproduce_sweeps.sh']:
 shutil.copy2(ROOT/name,STAGE/name)
(STAGE/'data').mkdir()
for p in (ROOT/'data').iterdir():
 if not p.is_file():continue
 keep=(p.name in ['DELIVERY_STATUS.json','calibration_audit.json','observable_validation.json','metrics_summary.json','metrics_summary.csv'] or ('N321' in p.name and p.name!='T100_N321_I1_blocking_bn0_bp0.json' and p.name!='T100_N321_I1_blocking_bn0_bp0.npz'))
 if keep:shutil.copy2(p,STAGE/'data'/p.name)
(STAGE/'data/EXPLORATORY_EXCLUSIONS.txt').write_text('N161 and absolute-quasi-Fermi exploratory attempts are retained in the original workspace but excluded from this final bundle. The live100Kdark continuation is represented only by its frozen partial snapshot. No results beyond snapshot coverage are claimed. The snapshot contains one final-gate failure at0.025V, explicitly marked and excluded from certified plotting.\n')
files={str(p.relative_to(STAGE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in STAGE.rglob('*') if p.is_file()}
(STAGE/'SHA256.json').write_text(json.dumps(files,indent=2))
out=ROOT.parent/'PM6Y6_reference_baseline_sources_20261004.zip'
with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED,6) as z:
 for p in STAGE.rglob('*'):
  if p.is_file():z.write(p,p.relative_to(STAGE.parent))
with zipfile.ZipFile(out) as z:assert z.testzip() is None
print(json.dumps(dict(path=str(out),bytes=out.stat().st_size,files=len(files),sha256=hashlib.sha256(out.read_bytes()).hexdigest())))
