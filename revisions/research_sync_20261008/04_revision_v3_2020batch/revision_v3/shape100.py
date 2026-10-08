import sys;sys.path.insert(0,'.')
from lib_v2 import *
zo=orig_states([275,250,225,200,175,150,125,100])
for tag in ('V1_cal','V1_one'):
    d,z=to_final(tag,100,zo[100]);z0=z;V=0.;rows=[(0.,Jmean(d,z,0.))]
    while V<1.13:
        V1=round(V+.04,6);z=adv(d,z,V,V1);V=V1;o=d.evaluate(z,V);rows.append((V,Jmean(d,z,V),float(np.median(o['n'][(d.x>.2)&(d.x<.8)]))))
    print(tag);print(' V      J(mA/cm2)  n_bulk   V*J')
    for r in rows:print(f' {r[0]:.2f}  {-r[1]*1e3:8.3f}  {(r[2] if len(r)>2 else float("nan")):.2e}  {-r[0]*r[1]*1e3:7.3f}')
