from analyze_density import *
from rate_audit import audit
import hashlib,platform,scipy

def main():
 runs=[];mesh=[];steps=[];equivs=[];ratios=[];windows=[];frozen=[]
 for p in sorted((ROOT/'data').glob('PF_Nt*.json'))+sorted((ROOT/'data').glob('DD_Nt*.json')):
  d=json.loads(p.read_text());rows=d['rows']+d.get('forward_rows',[]);tr=d['trace'];runs.append(dict(file=p.name,rows=len(rows),failure=d['failure'],max_scaled_residual=max(r['scaled_residual'] for r in rows),max_relative_current_error=max(r['relative_error'] for r in rows),max_charge_relative_error=max(r['charge_relative_error'] for r in rows),all_gates=all(r['gate_passed'] for r in rows),failed_continuation_attempts=sum('status' in r for r in tr),equilibrium_J_Acm2=d['equilibrium']['J_Acm2'],equilibrium_pair_Acm2=d['equilibrium']['pair_source_Acm2']))
 for L in [0,1]:
  for Nt in [1e14,1e15,1e16]:
   a=read('PF',Nt,L)['rows'];b=read('PF',Nt,L,161,.5)['rows'];ra=read('PF',1e15,L)['rows'];rb=read('PF',1e15,L,161,.5)['rows']
   for U in [5,10,15,20,23.5,25,30]:
    ar=min(a,key=lambda r:abs(r['U']-U));br=min(b,key=lambda r:abs(r['U']-U));mesh.append(dict(Nt_cm3=Nt,light=L,U=U,J_N81_Acm2=ar['J_Acm2'],J_N161_Acm2=br['J_Acm2'],relative_change=br['J_Acm2']/ar['J_Acm2']-1,field_center_relative_change=br['field_center_Vcm']/ar['field_center_Vcm']-1))
    rar=min(ra,key=lambda r:abs(r['U']-U));rbr=min(rb,key=lambda r:abs(r['U']-U));v1=(ar['pair_source_Acm2']/Nt)/(rar['pair_source_Acm2']/1e15);v2=(br['pair_source_Acm2']/Nt)/(rbr['pair_source_Acm2']/1e15);ratios.append(dict(Nt_cm3=Nt,light=L,U=U,ratio_N81=v1,ratio_N161=v2,ratio_absolute_mesh_change=v2-v1))
   arr=np.load(ROOT/'data'/f'PF_Nt{Nt:.0e}_L{L}_N81_h0.25.npz');d=DensityDevice(PFParams(**dict(BASE,Nt=Nt)),81,light=L)
   for U in [5,20,25]:
    k=np.argmin(abs(arr['V']+U));z=arr['z'][k];info=d.physical(z,-U);ev,_=audit(d,d.y_from_z(z,-U),-U,L);equivs.append(dict(Nt_cm3=Nt,light=L,U=U,J_direct_Acm2=ev['J_mAcm2']/1000,J_stable_Acm2=info['J_Acm2'],J_relative_difference=(ev['J_mAcm2']/1000-info['J_Acm2'])/max(abs(info['J_Acm2']),1e-300),pair_direct_Acm2=ev['pair_net_mAcm2']/1000,pair_stable_Acm2=info['pair_source_Acm2'],pair_relative_difference=(ev['pair_net_mAcm2']/1000-info['pair_source_Acm2'])/max(abs(info['pair_source_Acm2']),1e-300),independent_displacement_error_Acm2=ev['displacement_current_max_error_Acm2']))
   if L==0:
    h=read('PF',Nt,0,81,.125)['rows']
    for ar in a[1:]:
     br=min(h,key=lambda r:abs(r['U']-ar['U']));steps.append(dict(Nt_cm3=Nt,U=ar['U'],relative_J_change=br['J_Acm2']/ar['J_Acm2']-1))
   fn=ROOT/'data'/f'window_PF_Nt{Nt:.0e}_L{L}_N161_h0.125.json'
   if fn.exists():
    ww=json.loads(fn.read_text())['rows'];x=np.array([r['U'] for r in ww]);y=-np.array([r['J_Acm2'] for r in ww])
    for w in [5,9,13]:
     dd=savgol_filter(y,w,3,deriv=2,delta=.125);k=np.argmax(dd[5:-5])+5;windows.append(dict(Nt_cm3=Nt,light=L,nodes=161,step=.125,window_samples=w,support_V=(w-1)*.125,largest_curvature_U=float(x[k])))
 for name,data in [('mesh_checks.csv',mesh),('step_checks.csv',steps),('density_ratio_mesh_checks.csv',ratios),('independent_event_checks.csv',equivs),('refined_landmarks.csv',windows)]:csvout(name,data)
 summary=dict(runs=runs,mesh=mesh,steps_max_relative=max(abs(r['relative_J_change']) for r in steps),density_ratio_max_absolute_mesh_change_10to20=max(abs(r['ratio_absolute_mesh_change']) for r in ratios if 10<=r['U']<=20),independent_event_max_relative_J=max(abs(r['J_relative_difference']) for r in equivs),independent_event_max_relative_pair=max(abs(r['pair_relative_difference']) for r in equivs),refined_landmarks=windows)
 save(ROOT/'data/validation.json',summary);print(json.dumps({k:v for k,v in summary.items() if k not in ['runs','mesh']},indent=2))
 save(ROOT/'provenance.json',dict(date_UTC='2026-10-02',python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,source_project='/workspace/shared/wxh/002高YZ/OSC_FF_reverse_knee_model',source_qf='/workspace/shared/wxh/research/qf_temperature_confidence_20261002',design_sha256=hashlib.sha256((ROOT/'DESIGN_BEFORE_COMPUTING.md').read_bytes()).hexdigest(),baseline_sha256=hashlib.sha256((ROOT/'baseline_ff78.json').read_bytes()).hexdigest(),changed_original_sources=False))
if __name__=='__main__':main()
