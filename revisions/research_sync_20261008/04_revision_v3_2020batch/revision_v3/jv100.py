import sys,json,time;sys.path.insert(0,'.')
from lib_v2 import *
zo=orig_states([275,250,225,200,175,150,125,100])
out={}
for tag in ('V1_cal','V1_one','V2_cal','V2_one'):
    for T in (125,100):
        t=time.time()
        try:
            d,z=to_final(tag,T,zo[T]);vc,FF,Jsc=jv(d,z);out[f'{tag}_{T}']=dict(Jsc=Jsc,Voc=vc,FF=FF);print(tag,T,f'Jsc={Jsc:.3f} Voc={vc:.4f} FF={100*FF:.1f}%',f'{time.time()-t:.0f}s',flush=True)
        except Exception as e:print(tag,T,'FAIL',repr(e)[:150],flush=True)
json.dump(out,open('results/jv_100_125.json','w'),indent=1)
