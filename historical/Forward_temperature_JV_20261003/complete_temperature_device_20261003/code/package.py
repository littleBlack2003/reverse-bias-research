from pathlib import Path
import shutil,json,hashlib,zipfile
R=Path(__file__).resolve().parents[1];P=R/'package';P.mkdir(exist_ok=True)
# Preserve relative import layout, with only scientific dependencies and deliverables.
for f in R.rglob('*'):
 if not f.is_file() or any(x in f.parts for x in ['package','__pycache__']):continue
 rel=f.relative_to(R)
 if rel.name in ['library_file_transfer.py','summary_progress.json','closure_review.md','LIBRARY_UPLOAD.json']:continue
 if rel.suffix not in ['.py','.json','.txt','.csv','.npz','.png','.pdf','.log','.sh']:continue
 out=P/R.name/rel;out.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(f,out)
for name in ['gaussian_transport_20261003','coupled_nonlocal_300K_20261003','localselfconsistent_srh_300K_20261003']:
 D=R.parent/name
 for f in D.rglob('*'):
  if not f.is_file() or '__pycache__' in f.parts:continue
  rel=f.relative_to(D)
  if f.suffix!='.py' and f.name!='baseline_ff78.json':continue
  out=P/name/rel;out.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(f,out)
D=P/'S10_arrhenius_digitization_20261003';D.mkdir(exist_ok=True)
for name in ['digitized_currents.csv','analysis.json','curvature_check.json','S10温度电流读图报告.txt']:
 shutil.copy2(R.parent/'S10_arrhenius_digitization_20261003'/name,D/name)
(P/'gaussian_transport_20261003/data').mkdir(exist_ok=True)
(P/'README.txt').write_text('Conditional complete-device temperature predictions, 2026-10-03\n\nPython 3 + numpy/scipy/matplotlib. No installation hooks. Run from any directory:\npython complete_temperature_device_20261003/code/run.py 81 300 .4\nFor the full grid: python complete_temperature_device_20261003/code/batch.py\nThe separate originally run probe was: python complete_temperature_device_20261003/code/run.py 161 420 .4\npython complete_temperature_device_20261003/code/controls.py\npython complete_temperature_device_20261003/code/mu_controls.py\npython complete_temperature_device_20261003/code/mu_refine.py\npython complete_temperature_device_20261003/code/refine641.py\nSet OPENBLAS_NUM_THREADS=1 and MPLCONFIGDIR to a writable directory.\n\nNo calibrated prediction of a real sample is claimed. Read the Chinese report and parameter/context limitations before using figures. All temperatures share material/contact parameters. 140 final saved states include equilibria, controls and 641-node refinements. S10 values are pixel digitization, not original measurement arrays.\n',encoding='utf-8')
manifest={str(f.relative_to(P)):hashlib.sha256(f.read_bytes()).hexdigest() for f in P.rglob('*') if f.is_file() and f.name!='SHA256.json'}
(P/'SHA256.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False))
out=R.parent/'Complete_device_temperature_predictions_20261003.zip'
with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED) as z:
 for f in P.rglob('*'):
  if f.is_file():z.write(f,f.relative_to(P))
with zipfile.ZipFile(out) as z:assert z.testzip() is None
print(out,out.stat().st_size,len(manifest))
