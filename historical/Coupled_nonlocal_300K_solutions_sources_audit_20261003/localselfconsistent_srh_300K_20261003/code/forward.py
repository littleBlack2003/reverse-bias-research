from integrate import *
d=make(161);z=np.load(ROOT/'data/baseline_N161_V0.npz')['z'];trace=[]
for L in np.linspace(0,1,21)[1:]:
 d.light=L;z,info=solve_stable(d,0,z);assert not info.get('status'),info
rows=[];old=0.
for target in [0.,.8297514975955914,.9512297784607306]:
 for V in np.linspace(old,target,max(1,int(np.ceil((target-old)/.02)))+1)[1:]:z,info=advance(d,z,old,float(V),trace);old=float(V)
 r,a=details(d,z,target,info);rows.append(r);print(target,r['J_Acm2'],r['relative_error'],flush=True)
ref=json.loads((ROOT/'baseline_ff78.json').read_text())['achieved']
Jsc=-rows[0]['J_Acm2']*1000;Jmpp=-rows[1]['J_Acm2']*1000
out=dict(rows=rows,reference=ref,Jsc_mAcm2=Jsc,J_at_reference_mpp_mAcm2=Jmpp,FF_at_reference_points=.8297514975955914*Jmpp/(.9512297784607306*Jsc),Voc_point_current_Acm2=rows[2]['J_Acm2'],trace=trace)
save(ROOT/'data/forward_check.json',out)
assert abs(Jsc/ref['Jsc_mAcm2']-1)<1e-7
assert abs(Jmpp/ref['Jmpp_mAcm2']-1)<1e-7
assert abs(rows[2]['J_Acm2'])<1e-8
