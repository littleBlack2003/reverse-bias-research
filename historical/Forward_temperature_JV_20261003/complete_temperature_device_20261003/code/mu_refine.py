from pathlib import Path
import subprocess,os
root=Path(__file__).resolve().parents[1]
for T in [340,380,420]:
 with open(root/'logs'/f'mu300_T{T}_N161.log','w') as f:
  p=subprocess.run(['python',str(root/'code/run.py'),'161',str(T),'.4','gaussian_mu300'],stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'))
 print(T,p.returncode,flush=True)
