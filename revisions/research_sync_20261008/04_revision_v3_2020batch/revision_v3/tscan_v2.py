import sys,time,json;sys.path.insert(0,'.')
from run_v2 import *
from extend_device import temperature_continue
Ts=[300,275,250,225,200,175,150,125,100]
z=saved(300)['sc'][2];trace=[];zo={300:z};prev=300
for T in Ts[1:]:
    z=temperature_continue(z,prev,T,0.,1.,321,trace);zo[T]=z;prev=T
print('orig-model states ok',flush=True)
zeta={a:calibrate_zeta(a) for a in ('cal','one')}
res={}
for T in Ts:
    d0=make_device(T,321);o=d0.evaluate(zo[T],0.);J0=-float(np.mean(o['Jn']+o['Jp']))*1e3
    row={'orig':J0}
    for amp in ('cal','one'):
        for name,(ld,lb) in {'V1':(1.,0.),'V2':(1.,1.)}.items():
            try:
                d,zz=solve_state(T,amp,ld,lb,zeta[amp],start=zo[T]);m=metrics_at(d,zz);row[f'{name}_{amp}']=m['Jsc_mAcm2'];row[f'{name}_{amp}_nbulk']=m['n_bulk_median']
            except Exception as e:row[f'{name}_{amp}']=None
    res[T]=row;print(T,{k:(None if v is None else float(f'{v:.4g}')) for k,v in row.items()},flush=True)
json.dump(res,open('results/jsc_vs_T.json','w'),indent=1)
