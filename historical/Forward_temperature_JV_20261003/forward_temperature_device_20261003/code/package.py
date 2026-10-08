from pathlib import Path
import zipfile,hashlib,json
R=Path(__file__).resolve().parents[1];files=[]
for f in R.rglob('*'):
 if f.is_file() and not any(x in f.parts for x in ['helpers','__pycache__']) and f.suffix in ['.py','.json','.npz','.csv','.txt','.png','.pdf','.log'] and not f.name.startswith('LIBRARY_'):files.append(f)
for dirname in ['complete_temperature_device_20261003','gaussian_transport_20261003','coupled_nonlocal_300K_20261003','localselfconsistent_srh_300K_20261003']:
 for f in (R.parent/dirname).rglob('*'):
  if f.is_file() and 'package' not in f.parts and '__pycache__' not in f.parts and (f.suffix=='.py' or f.name=='baseline_ff78.json') and 'library_' not in f.name:files.append(f)
manifest={str(f.relative_to(R.parent)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files}
out=R.parent/'Forward_temperature_JV_20261003.zip'
with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED) as z:
 for f in files:z.write(f,str(f.relative_to(R.parent)))
 z.writestr('SHA256.json',json.dumps(manifest,indent=2))
 z.writestr('README.txt','Same-parameter dark forward JV, 300/340/380/420 K, b=0.4/0.5 eV.\nRead Chinese report for assumptions and high-current limitations.\nRequires Python numpy scipy matplotlib. OPENBLAS_NUM_THREADS=1 python forward_temperature_device_20261003/code/run.py 161 300 .4\nNo illumination, solar FF or experimental validation is claimed.\n')
with zipfile.ZipFile(out) as z:assert z.testzip() is None
print(out,out.stat().st_size,len(files))
