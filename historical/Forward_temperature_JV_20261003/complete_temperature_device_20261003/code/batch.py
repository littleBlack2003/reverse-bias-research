from pathlib import Path
import subprocess,concurrent.futures,os
root=Path(__file__).resolve().parents[1]
cases=[(N,T,b) for b in [.4,.5] for T in [300,340,380,420] for N in [81,161,321] if (N,T,b)!=(161,420,.4)]
def run(c):
 N,T,b=c; tag=f'b{b:g}_T{T}_N{N}'
 with open(root/'logs'/f'{tag}.log','w') as f:
  p=subprocess.run(['python',str(root/'code/run.py'),str(N),str(T),str(b)],stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'))
 print(tag,p.returncode,flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:list(ex.map(run,cases))
