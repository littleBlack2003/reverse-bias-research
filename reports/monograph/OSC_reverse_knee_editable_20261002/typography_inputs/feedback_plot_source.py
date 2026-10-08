#!/usr/bin/env python3
"""Next justified diagnostic: isolate Coulomb-conditioned extraction feedback.
Only comparison switches Coulomb on/off; baseline H, geometry, lambda, T unchanged.
This control is NOT a new material parameter fit.
"""
from pathlib import Path
import os
os.environ.setdefault('MPLBACKEND','Agg');os.environ.setdefault('MPLCONFIGDIR','/tmp/quantum_cycle_mpl')
import csv,json
import numpy as np
import matplotlib.pyplot as plt
import model as m
ROOT=Path(__file__).parent;R=ROOT/'results';FIG=ROOT/'figures'
def write(name,rr):
 with (R/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
rows=[];edgerows=[];A0=m.A0
for scale in [0,1]:
 m.A0=A0*scale
 for kind in ['classical','quantum']:
  for gamma in np.logspace(2,10,161):
   row,p,rat,Q,ed,*_=m.solve(2,80,kind,gamma=float(gamma),detail=True)
   row['coulomb_scale']=scale;rows.append(row)
  for gamma in [1e2,1e4,1e6,1e8,1e14]:
   row,p,rat,Q,ed,*_=m.solve(2,80,kind,gamma=gamma,detail=True)
   for (i,j,name),(kf,kr,dg,tdg,dn),e in zip(m.EDGES,rat,ed):
    edgerows.append(dict(coulomb_scale=scale,T_K=80,F_MVcm=2,kernel=kind,Gamma_s=gamma,source=''.join(map(str,m.STATES[i])),destination=''.join(map(str,m.STATES[j])),process=name,source_probability=float(p[i]),destination_probability=float(p[j]),delta_energy_eV=float(dg),forward_rate_s=float(kf),reverse_rate_s=float(kr),net_edge_flux_s=float(e[5])))
m.A0=A0
write('extraction_coulomb_feedback_scan.csv',rows);write('extraction_coulomb_edge_fluxes.csv',edgerows)
summary=[]
for scale in [0,1]:
 for kind in ['classical','quantum']:
  cr=[r for r in rows if r['coulomb_scale']==scale and r['kernel']==kind];rates=np.array([r['R_s'] for r in cr]);pk=cr[np.argmax(rates)]
  tail=cr[-1];summary.append(dict(coulomb_scale=scale,kernel=kind,peak_Gamma_s=pk['GammaC_s'],peak_R_s=pk['R_s'],Gamma1e10_R_s=tail['R_s'],peak_over_fast_extraction=pk['R_s']/tail['R_s'],negative_slope_intervals=int(sum(np.diff(rates)<-1e-8*max(rates)))))
(R/'extraction_coulomb_summary.json').write_text(json.dumps(summary,indent=2))
fig,axes=plt.subplots(1,2,figsize=(10.8,4.3))
for scale in [0,1]:
 for kind,ls in [('classical','--'),('quantum','-')]:
  cr=[r for r in rows if r['coulomb_scale']==scale and r['kernel']==kind];ax=axes[scale]
  ax.loglog([r['GammaC_s'] for r in cr],[r['R_s'] for r in cr],ls,label=kind)
 axes[scale].set(xlabel='Equal endpoint exchange Gamma (s$^{-1}$)',ylabel='Net cycle rate (s$^{-1}$)',title='Coulomb '+('included' if scale else 'off control'));axes[scale].grid(alpha=.2);axes[scale].legend()
fig.suptitle('80 K, 2 MV/cm; fixed H, geometry and reorganization parameters')
plt.tight_layout();plt.savefig(FIG/'coulomb_extraction_feedback.png',dpi=170);plt.close()
print(json.dumps(summary,indent=2))
