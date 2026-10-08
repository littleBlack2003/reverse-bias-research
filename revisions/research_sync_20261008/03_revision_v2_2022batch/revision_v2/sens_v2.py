import sys,os,json;sys.path.insert(0,'.')
import run_v2
from run_v2 import *
out={}
for amp in ('cal','one'):
    z0=calibrate_zeta(amp)
    for mult in (.3,1.,3.):
        d,z=solve_state(100,amp,1.,1.,z0*mult);m=metrics_at(d,z);out[f'V2_{amp}_zeta x{mult}']=m['Jsc_mAcm2'];print('V2',amp,'zeta x',mult,round(m['Jsc_mAcm2'],3),'R/G',round(m['R_over_G'],3),flush=True)
for cc in ('0.05','0.2'):
    os.environ['CLOSURE_CCUT']=cc
    for amp in ('cal','one'):
        d,z=solve_state(100,amp,1.,0.,None);m=metrics_at(d,z);out[f'V1_{amp}_ccut{cc}']=m['Jsc_mAcm2'];print('V1',amp,'c_cut',cc,round(m['Jsc_mAcm2'],3),flush=True)
json.dump(out,open('results/sens_100K.json','w'),indent=1)
