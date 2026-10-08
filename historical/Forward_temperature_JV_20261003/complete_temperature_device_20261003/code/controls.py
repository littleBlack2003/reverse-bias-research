from pathlib import Path
import subprocess,concurrent.futures,os
root=Path(__file__).resolve().parents[1]
def run(T):
 with open(root/'logs'/f'constmu_T{T}_N321.log','w') as f:
  p=subprocess.run(['python',str(root/'code/run.py'),'321',str(T),'.4','gaussian_constmu'],stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'))
 print(T,p.returncode,flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:list(ex.map(run,[300,340,380,420]))
