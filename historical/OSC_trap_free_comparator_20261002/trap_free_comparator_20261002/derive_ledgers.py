from pathlib import Path
import csv,json,hashlib
import comparator_model as m
root=Path(__file__).parent

def write(name,rows):
 with (root/'results'/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
rows=[]
for r in [1.,2.,5.,10.]:
 bind=m.COULOMB/3.5/r
 rows.append(dict(span_nm=r,binding_eV=bind,gap_eV=1.3,field_for_deltaG_zero_MVcm=(1.3-bind)/(.1*r),field_for_l0_activationless_MVcm=(1.3-bind+.05)/(.1*r),warning='point-charge and uncalibrated energy parameters; energy alignment only, not upturn or breakdown prediction'))
write('geometry_energy_screen.csv',rows)
row,p,Q,rates,details,E,states,edges=m.solve(1.5,80,'direct',detail=True)
idx={''.join(map(str,s)):i for i,s in enumerate(states)};led=[]
for a,b,process,work in [('10','01','V_C',0),('01','00','bathR',.75),('00','10','bathL',.75)]:
 de=float(E[idx[b]]-E[idx[a]])
 led.append(dict(from_state=a,to_state=b,process=process,delta_state_energy_eV=de,chemical_work_into_system_eV=work,total_heat_to_baths_eV=work-de))
write('direct_cycle_energy_ledger.csv',led)
assert abs(sum(x['delta_state_energy_eV'] for x in led))<1e-14
assert abs(sum(x['total_heat_to_baths_eV'] for x in led)-1.5)<1e-14
print('Cycle energy ledger and geometry screen passed')
