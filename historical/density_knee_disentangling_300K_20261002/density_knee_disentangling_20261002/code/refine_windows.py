from run_density import *
def run(Nt,L):
 N=161;old=f'PF_Nt{Nt:.0e}_L{L}_N161_h0.5';a=np.load(ROOT/'data'/f'{old}.npz');idx=np.argmin(abs(a['V']+21));z=a['z'][idx].copy();p=dict(BASE,Nt=Nt);d=DensityDevice(PFParams(**p),N,light=L);trace=[];rows=[];states=[];prev=-21.
 for U in np.arange(21,26.001,.125):
  z,info=advance(d,z,prev,-float(U),trace);prev=-float(U);row,_=details(d,z,-float(U),info);rows.append(row);states.append(z.copy())
 name=f'window_PF_Nt{Nt:.0e}_L{L}_N161_h0.125';save(ROOT/'data'/f'{name}.json',dict(rows=rows,trace=trace,failure=None));np.savez_compressed(ROOT/'data'/f'{name}.npz',z=states,V=[r['V'] for r in rows]);print(name,'done',flush=True)
if __name__=='__main__':run(float(sys.argv[1]),int(sys.argv[2]))
