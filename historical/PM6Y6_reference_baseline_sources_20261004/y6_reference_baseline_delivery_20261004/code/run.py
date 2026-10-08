from pathlib import Path
from dataclasses import asdict
import json,sys,time,traceback,os
import numpy as np
from device import *
from calibration import beta,G0,bulk_voc
ROOT=Path(__file__).resolve().parents[1]
def default(obj):
    if isinstance(obj,np.generic):return obj.item()
    raise TypeError(type(obj).__name__)
def run(T,N=161,I=1.,contact='blocking',bn=0.,bp=0.,Eg=1.42,step=.025):
    tag=f'T{T:g}_N{N}_I{I:g}_{contact}_bn{bn:g}_bp{bp:g}'+('_lightonly' if os.environ.get('ONLY_LIGHT')=='1' else '')
    p=Parameters(T=T,allow_mobility_extrapolation=(T<223 or T>328),beta_cm3s=float(beta(T)),G_cm3s=G0,minority_contact=contact,barrier_n_eV=bn,barrier_p_eV=bp,Eg_eV=Eg)
    d=Device(p,N);start=time.monotonic();trace=[];curves={};failure=None;states=[];labels=[]
    def persist():
        result=dict(tag=tag,parameters=asdict(p),state_coordinates='contact_relative_qf_v2',status='conditional simulation, not calibrated experimental JV',curves=curves,failure=failure,seconds=time.monotonic()-start,trace=trace)
        (ROOT/'data'/f'{tag}.json').write_text(json.dumps(result,indent=2,default=default))
        if states:np.savez_compressed(ROOT/'data'/f'{tag}.npz',z=np.array(states),labels=np.array(labels),x_cm=d.x*p.d_cm)
    failures=[]
    try:
        zeq,info=solve(d,0,d.initial());trace.append({k:v for k,v in info.items() if k!='history'})
        if not info['success']:raise RuntimeError('Equilibrium: '+info['reason'])
        for L in ([I] if os.environ.get("ONLY_LIGHT")=="1" else [0.,I]):
            arr=[];curves['dark' if not L else 'light']=arr
            try:
                if L:
                    vs=float(bulk_voc(T,I,Eg));z=poisson_seed(d,vs)
                    z,info=solve(d,vs,z,L,maxiter=240)
                    if not info['success']:raise RuntimeError('Light seed: '+info['reason'])
                    zs=z.copy();old=vs
                    targets=list(np.round(np.arange(0,vs,step),9)[::-1])
                    branches=[(zs,vs,targets),(zs,vs,list(np.round(np.arange(np.ceil(vs/step)*step,1.100001,step),9)))]
                else:branches=[(zeq,0.,list(np.round(np.arange(0,1.100001,step),9)))]
                for seed,old,targets in branches:
                    z=seed.copy()
                    for V in targets:
                        if V!=old:z,info=advance(d,z,old,float(V),L,L,trace)
                        old=float(V);out=d.ledger(z,V,L);out['newton_norm']=info.get('norm');arr.append(out)
                        states.append(z.copy());labels.append([L,V]);persist()
                        print(tag,L,V,out['J_Acm2'],out['continuity_relative'],out['energy_relative'],flush=True)
                arr.sort(key=lambda x:x['V'])
            except Exception as e:
                failures.append(dict(branch='dark' if not L else 'light',message=str(e),traceback=traceback.format_exc()))
                print('BRANCH FAILED',tag,L,str(e),flush=True)
    except Exception as e:failures.append(dict(branch='equilibrium',message=str(e),traceback=traceback.format_exc()))
    failure=failures or None
    persist();return tag
if __name__=='__main__':run(float(sys.argv[1]),int(sys.argv[2]) if len(sys.argv)>2 else 161,float(sys.argv[3]) if len(sys.argv)>3 else 1.,sys.argv[4] if len(sys.argv)>4 else'blocking',float(sys.argv[5]) if len(sys.argv)>5 else 0.,float(sys.argv[6]) if len(sys.argv)>6 else 0.,step=float(sys.argv[7]) if len(sys.argv)>7 else .025)
