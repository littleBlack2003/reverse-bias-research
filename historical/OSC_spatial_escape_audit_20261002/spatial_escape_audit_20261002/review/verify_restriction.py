"""Verify the added exclusion model by reducing the original full generator."""
from pathlib import Path
import sys,json,hashlib
import mpmath as mp
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'snapshot'))
import spatial_model as sm

rows=[]
for F in [1.5,2.,3.]:
    full=sm.solve(F,80,detail=True,dps=250)
    row,p,Q,rates,details,E,mu,g,states,edges=full
    keep=[i for i,s in enumerate(states) if 2-s[0]-s[1]<=1 and s[2]+s[3]+s[4]<=1]
    loc={full_i:i for i,full_i in enumerate(keep)}
    Qr=mp.zeros(len(keep))
    for i,ii in enumerate(keep):
        for j,jj in enumerate(keep):
            if i!=j:Qr[i,j]=Q[ii,jj]
    # Remove outgoing events as well as their reverse; rebuild diagonal losses.
    for j in range(len(keep)):Qr[j,j]=-sum(Qr[i,j] for i in range(len(keep)) if i!=j)
    A=Qr.copy();bb=mp.matrix([0]*len(keep));bb[len(keep)-1]=1
    for j in range(len(keep)):A[len(keep)-1,j]=1
    pr=mp.lu_solve(A,bb)
    left=mp.mpf(0);right=mp.mpf(0);heat=mp.mpf(0)
    for (i,j,ty,k),rate in zip(edges,rates):
        if i not in loc or j not in loc:continue
        current=pr[loc[i]]*rate['kf']-pr[loc[j]]*rate['kr']
        heat-=current*rate['tdg']
        if ty=='bath':
            if k==0:left+=current
            else:right-=current
    direct=sm.solve(F,80,variant='direct',dps=250)
    rows.append({'F_MVcm':F,'T_K':80,'full32_rate_s':row['R_s'],'direct_rate_s':direct['R_s'],
                 'restricted12_rate_s':float(left),'number_restricted_states':len(keep),
                 'normalization_error':float(abs(sum(pr)-1)), 'minimum_probability':float(min(pr)),
                 'continuity_relative':float(abs(left-right)/abs(left)),
                 'heat_chemical_relative':float(abs(heat-mp.mpf(str(F))*left)/(mp.mpf(str(F))*abs(left)))})
for row in rows:
    assert row['number_restricted_states']==12
    assert row['minimum_probability']>0
    for key in ['normalization_error','continuity_relative','heat_chemical_relative']:assert row[key]<1e-100,(key,row[key])
(HERE/'restricted_graph_independent.json').write_text(json.dumps(rows,indent=2)+'\n')
print(json.dumps(rows,indent=2))
