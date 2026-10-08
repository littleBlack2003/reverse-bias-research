from pathlib import Path
import json,hashlib,zipfile
ROOT=Path(__file__).resolve().parents[1]
ALLOW={'code','reference','data','figures','validation','archive_record','input_with_shunt_saved_states'}
ROOT_FILES={'README.txt','requirements.lock.txt','runtime.json','PM6Y6_去除外部分流_重算说明.txt'}
files=[]
for f in ROOT.rglob('*'):
 if not f.is_file():continue
 rel=f.relative_to(ROOT)
 if '__pycache__' in rel.parts or f.suffix=='.pyc' or f.name=='reprocheck_path.txt' or f.suffix=='.log':continue
 if rel.parts[0] in ALLOW or str(rel) in ROOT_FILES:files.append(f)
files=sorted(files)
manifest={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files}
(ROOT/'SHA256.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
zip_path=ROOT.parent/'PM6Y6_no_shunt_100-300K_reproducible_20261004.zip'
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for f in files+[ROOT/'SHA256.json']:z.write(f,ROOT.name+'/'+str(f.relative_to(ROOT)))
with zipfile.ZipFile(zip_path) as z:
 assert z.testzip() is None
 assert set(z.namelist())=={ROOT.name+'/'+k for k in manifest}|{ROOT.name+'/SHA256.json'}
 for k,h in manifest.items():assert hashlib.sha256(z.read(ROOT.name+'/'+k)).hexdigest()==h
summary={'archive':zip_path.name,'bytes':zip_path.stat().st_size,'sha256':hashlib.sha256(zip_path.read_bytes()).hexdigest(),'files':len(files)+1,'CRC_passed':True,'all_file_hashes_passed':True}
(ROOT/'archive_validation.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
