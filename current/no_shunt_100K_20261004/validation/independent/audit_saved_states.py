"""Independent no-solve audit of saved states, EOS, precision and changed circuit.
Never edits numerical-core files or runs a Poisson/drift-diffusion solve.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import argparse,csv,hashlib,json,sys
from pathlib import Path
from functools import lru_cache
import numpy as np
from scipy.special import roots_legendre,expit
ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--prefix',default='completed');ap.add_argument('--source',type=Path);ap.add_argument('--temperatures',nargs='*',type=float);args=ap.parse_args()
ROOT=args.root.resolve();OUT=args.out.resolve();OUT.mkdir(parents=True,exist_ok=True)
# Portable archived inputs take priority by default; no legacy path is assumed.
if args.source is None:
 args.source=ROOT/'input_with_shunt_saved_states'
if not args.source.is_dir():
 raise FileNotFoundError(f'Prior input tree is required for frozen-input verification: {args.source}. Supply --source or bundle input_with_shunt_saved_states.')
args.source=args.source.resolve()
sys.path.insert(0,str(ROOT/'code'))
from device import Parameters,Q
from field_device import FieldDevice

def write(name,rows):
 if rows:
  with (OUT/(args.prefix+'_'+name+'.csv')).open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
@lru_cache(None)
def rule(n):
 z,w=roots_legendre(n);z*=18;w*=18*np.exp(-z*z/2)/np.sqrt(2*np.pi)
 # Exact Gaussian mass over this window differs from 1 by <2e-72.
 # Normalize only quadrature-rule roundoff, report its original error separately.
 norm=float(w.sum());return z,w/norm,norm

def ref(eta,s,order=4097):
 shape=np.shape(eta);eta=np.asarray(eta).ravel();z,w,_=rule(order);c=np.empty_like(eta);dc=np.empty_like(eta)
 for k in range(0,len(eta),128):
  sl=slice(k,k+128);a=-abs(eta[sl,None])-s*z;f=expit(a);cn=np.sum(f*w,axis=1)
  c[sl]=np.where(eta[sl]>0,1-cn,cn);dc[sl]=np.sum(f*expit(-a)*w,axis=1)
 return c.reshape(shape),dc.reshape(shape)

def independent_flux(eta,psi,qf,eos,mu,vt,dx,carrier,factor):
 # Reference EOS is quadrature, without spline interpolation.
 c,dc=ref(eta,eos.s);lc=np.log(c.astype(np.longdouble));de=np.diff(np.asarray(eta,np.longdouble));dl=np.diff(lc)
 cm,dm=ref((eta[:-1]+eta[1:])/2,eos.s)
 g=np.asarray(cm/dm,np.longdouble);large=abs(de)>1e-7;g[large]=de[large]/dl[large]
 a=np.diff(np.asarray(qf,np.longdouble))/g;dp=np.diff(np.asarray(psi,np.longdouble));x=(-dp if carrier=='n' else dp)/g
 b=np.empty_like(x);zero=abs(x)<1e-9;positive=x>50;negative=x<-50;mid=~(zero|positive|negative)
 b[zero]=np.log1p(-x[zero]/2+x[zero]**2/12)
 b[positive]=np.log(x[positive])-x[positive]-np.log1p(-np.exp(-x[positive]))
 b[negative]=np.log(-x[negative])-np.log1p(-np.exp(x[negative]))
 b[mid]=np.log(x[mid]/np.expm1(x[mid]))
 out=np.zeros_like(a);nz=abs(a)>0
 la=np.log(np.longdouble(Q*mu*vt))+np.log(g)-np.log(np.asarray(dx,np.longdouble))+lc[:-1]+b
 out[nz]=(1 if carrier=='n' else -1)*np.sign(a[nz])*np.exp(la[nz]+np.maximum(a[nz],0)+np.log(-np.expm1(-abs(a[nz]))))
 return np.asarray(out*factor,float),np.asarray(g,float),c

state_rows=[];eos_rows=[];metric_rows=[];frozen_rows=[];pending=[]
for path in sorted((ROOT/'data').glob('*.json')):
 try:data=json.loads(path.read_text())
 except Exception:continue
 if not isinstance(data,dict) or 'parameters' not in data or 'critical_states' not in data:continue
 T=data['parameters']['T'];N=data.get('nodes',data['parameters'].get('nodes'))
 if args.temperatures and T not in args.temperatures:continue
 if not path.with_suffix('.npz').exists() or not data.get('success'):
  pending.append(path.name);continue
 p=Parameters(**data['parameters']);d=FieldDevice(p,N,gamma_n=data['gamma_n'],gamma_p=data['gamma_p'],field_cap=data['field_cap'])
 ss=np.load(path.with_suffix('.npz'));labels=ss['labels'];states=ss['z'];critical=data['critical_states'];tag=data.get('tag',path.stem)
 selected=[(k,v) for k,v in critical.items() if isinstance(v,dict) and 'V' in v]
 selected += [(f'bracket_{i}',r) for i,r in enumerate(critical.get('bracket',[]))]
 i=np.where(labels[:,0]==1)[0];ih=i[np.argmax(labels[i,1])];selected.append(('max_saved_bias',{'V':float(labels[ih,1]),'light':1}))
 for name,row in selected:
  V=float(row['V']);L=float(row.get('light',0 if name=='equilibrium' else 1));idx=np.where((labels[:,0]==L)&(abs(labels[:,1]-V)<1e-12))[0]
  if not len(idx):raise RuntimeError(f'Missing state {tag} {name}')
  z=states[idx[-1]];o=d.evaluate(z,V);ledger=d.ledger(z,V,L);J=o['Jn']+o['Jp'];mean=float(J.mean());unc=ledger['current_closure_uncertainty_Acm2']
  rn=np.diff(o['Jn'])-Q*d.vol*(o['R'][1:-1]-L*p.G_cm3s);rp=np.diff(o['Jp'])+Q*d.vol*(o['R'][1:-1]-L*p.G_cm3s)
  E=abs(np.diff(z[:,0])*d.vt/d.dx);Ec=np.minimum(E,d.field_cap);sq=(Ec**2+1)**.25-1
  fn=np.exp(d.gamma_n*sq);fp=np.exp(d.gamma_p*sq)
  jnr,gn,cn=independent_flux(o['eta_n'],z[:,0],z[:,1],d.en,d.mn*p.N0_cm3,d.vt,d.dx,'n',fn)
  jpr,gp,cp=independent_flux(o['eta_p'],z[:,0],z[:,2],d.ep,d.mp*p.N0_cm3,d.vt,d.dx,'p',fp)
  A=V/d.vt+z[:,1]+z[:,2];rr=np.zeros_like(A);nz=abs(A)>0
  rr[nz]=np.sign(A[nz])*np.exp(np.log(p.beta_cm3s)+2*np.log(p.N0_cm3)+np.log(cn[nz])+np.log(cp[nz])+np.maximum(-A[nz],0)+np.log(-np.expm1(-abs(A[nz]))))
  rnr=np.diff(jnr)-Q*d.vol*(rr[1:-1]-L*p.G_cm3s);rpr=np.diff(jpr)+Q*d.vol*(rr[1:-1]-L*p.G_cm3s)
  ld=np.longdouble
  j_continuity_ld=ld(o['Jn'][0])+ld(o['Jp'][-1])+ld(Q)*np.sum(np.asarray(d.vol,ld)*(np.asarray(o['R'][1:-1],ld)-ld(L)*ld(p.G_cm3s)))
  Al=ld(V)/ld(d.vt)+np.asarray(z[:,1],ld)+np.asarray(z[:,2],ld);Rl=np.zeros_like(Al);az=abs(Al)>0
  logRl=np.log(ld(p.beta_cm3s))+2*np.log(ld(p.N0_cm3))+np.asarray(d.en.logc(o['eta_n']),ld)+np.asarray(d.ep.logc(o['eta_p']),ld)
  logRl[az]+=np.maximum(-Al[az],ld(0))+np.log(-np.expm1(-abs(Al[az])))
  Rl[az]=np.sign(Al[az])*np.exp(logRl[az]);j_continuity_ld_reaction=ld(o['Jn'][0])+ld(o['Jp'][-1])+ld(Q)*np.sum(np.asarray(d.vol,ld)*(Rl[1:-1]-ld(L)*ld(p.G_cm3s)))
  eta_max=max(o['eta_n'].max(),o['eta_p'].max());cmax=max(cn.max(),cp.max());min_reaction=np.min(Q*o['R']*o['A']*d.vt)
  res=d.residual(z,V,L);jac=d.jacobian(z,V,L);scale=np.maximum(np.asarray(abs(jac).max(axis=1).toarray()).ravel(),1e-250)
  state_rows.append(dict(tag=tag,T_K=T,N=N,state=name,V_V=V,J_intrinsic_Acm2=mean,closure_uncertainty_Acm2=unc,sign_margin=abs(mean)/max(unc,1e-300),current_spread_Acm2=float(np.ptp(J)),max_individual_carrier_Acm2=max(abs(o['Jn']).max(),abs(o['Jp']).max()),max_carrier_sum_magnitude_Acm2=float((abs(o['Jn'])+abs(o['Jp'])).max()),local_continuity_L1_Acm2=float(max(abs(rn).sum(),abs(rp).sum())),highorder_eos_current_Acm2=float(np.mean(jnr+jpr)),highorder_eos_current_shift_Acm2=float(abs(np.mean(jnr+jpr)-mean)),highorder_eos_max_edge_current_shift_Acm2=float(np.max(abs((jnr+jpr)-J))),highorder_eos_local_continuity_L1_Acm2=float(max(abs(rnr).sum(),abs(rpr).sum())),highorder_eos_integrated_recombination_shift_Acm2=float(abs(Q*np.dot(d.vol,(rr-o['R'])[1:-1]))),max_eta=float(eta_max),max_occupation=float(cmax),max_core_occupation=float(max(cn[(d.x>=.05)&(d.x<=.95)].max(),cp[(d.x>=.05)&(d.x<=.95)].max())),max_generalized_Einstein_factor=float(max(gn.max(),gp.max())),scaled_nonlinear_residual=float(max(abs(res)/scale)),charge_relative=ledger['charge_relative'],energy_residual_Wcm2=ledger['transport_dissipation_Wcm2']+ledger['recombination_dissipation_Wcm2']-ledger['light_chemical_work_Wcm2']-V*mean,min_reaction_dissipation_Wcm3=float(min_reaction),transport_dissipation_Wcm2=ledger['transport_dissipation_Wcm2'],gate_passed=ledger['gate_passed'],saved_intrinsic_error_Acm2=abs(row.get('J_intrinsic_Acm2',mean)-mean),longdouble_continuity_current_Acm2=float(j_continuity_ld),longdouble_continuity_current_difference_Acm2=float(abs(j_continuity_ld-ld(mean))),longdouble_reaction_continuity_current_Acm2=float(j_continuity_ld_reaction),longdouble_reaction_current_difference_Acm2=float(abs(j_continuity_ld_reaction-ld(mean)))))
  for car,eos in [('n',d.en),('p',d.ep)]:
   et=o['eta_'+car];c,dc,g=eos.evaluate(et);cr,dr=ref(et,eos.s);c2,d2=ref(et,eos.s,2049)
   eos_rows.append(dict(tag=tag,T_K=T,N=N,state=name,carrier=car,sigma_over_kT=eos.s,min_eta=float(et.min()),max_eta=float(et.max()),max_c_relative_error=float(max(abs(c/cr-1))),max_dc_relative_error=float(max(abs(dc/dr-1))),max_g_relative_error=float(max(abs(g/(cr/dr)-1))),max_2049_4097_c_relative_error=float(max(abs(c2/cr-1))),max_2049_4097_dc_relative_error=float(max(abs(d2/dr-1))),occupation_valid=bool(np.all(cr>0)&np.all(cr<1)),einstein_bound=bool(np.all(g>=1-1e-12))))
 m=data.get('metrics',{});br=critical.get('bracket',[]);curve=data.get('curves',{}).get('light',[])
 if all(k in critical for k in ('sc','oc','mpp')) and len(br)==2:
  a,b=br;ja=a['J_terminal_Acm2'];jb=b['J_terminal_Acm2'];slope=(jb-ja)/(b['V']-a['V']);sc=critical['sc'];mp=critical['mpp'];v=np.array([r['V'] for r in curve]);j=np.array([r['J_terminal_Acm2'] for r in curve]);pow=-v*j
  metric_rows.append(dict(tag=tag,T_K=T,N=N,Jsc_mAcm2=-sc['J_terminal_Acm2']*1000,Voc_V=critical['oc']['V'],FF=mp['external_terminal_Wcm2']/(sc['J_terminal_Acm2']*critical['oc']['V']),Pmax_mWcm2=-mp['external_terminal_Wcm2']*1000,lower_bracket_V=a['V'],upper_bracket_V=b['V'],lower_J_Acm2=ja,upper_J_Acm2=jb,lower_closure_Acm2=a['current_closure_uncertainty_Acm2'],upper_closure_Acm2=b['current_closure_uncertainty_Acm2'],lower_sign_margin=abs(ja)/max(a['current_closure_uncertainty_Acm2'],1e-300),upper_sign_margin=abs(jb)/max(b['current_closure_uncertainty_Acm2'],1e-300),signed_resolved_bracket=bool(ja < -a['current_closure_uncertainty_Acm2'] and jb>b['current_closure_uncertainty_Acm2']),local_slope_Acm2V=slope,root_closure_V_estimate=max(critical['oc']['current_closure_uncertainty_Acm2'],a['current_closure_uncertainty_Acm2'],b['current_closure_uncertainty_Acm2'])/abs(slope),number_curve_samples=len(curve),strictly_increasing_voltage=bool(np.all(np.diff(v)>0)),number_sampled_sign_crossings=int(np.sum(j[:-1]*j[1:]<0)),mpp_no_sampled_underperformance=bool(-mp['external_terminal_Wcm2']>=max(pow)*(1-1e-10)),max_shunt_current_Acm2=max(abs(r.get('J_shunt_Acm2',0)) for r in curve),max_terminal_intrinsic_difference_Acm2=max(abs(r['J_terminal_Acm2']-r['J_intrinsic_Acm2']) for r in curve)))
 if args.source:
  oldpath=args.source/'data'/f'latest_T{T:g}_N{N}.json'
  if oldpath.exists():
   old=json.loads(oldpath.read_text());frozen_rows.append(dict(tag=tag,parameters_exact=old['parameters']==data['parameters'],gamma_n_exact=old['gamma_n']==data['gamma_n'],gamma_p_exact=old['gamma_p']==data['gamma_p'],field_cap_exact=old['field_cap']==data['field_cap'],Rs_zero=data.get('Rs_ohmcm2')==0,Jsc_difference_Acm2=critical['sc']['J_intrinsic_Acm2']-old['critical_states']['sc']['J_intrinsic_Acm2'],Voc_shift_V=m['Voc_V']-old['metrics']['Voc_V']))
write('states',state_rows);write('eos',eos_rows);write('metrics',metric_rows);write('frozen_inputs',frozen_rows)
summary=dict(scope='Read-only independent saved-state audit; direct 2049/4097-node normalized Gauss-Legendre FD quadrature on [-18,18], no PDE rerun. Direct-EOS residual is a constitutive perturbation check, not a re-solved solution.',completed_tags=sorted(set(r['tag'] for r in state_rows)),pending=pending,states=len(state_rows),eos_tests=len(eos_rows),rule_original_gaussian_mass={k:rule(k)[2] for k in (2049,4097)},all_saved_state_gates=all(r['gate_passed'] for r in state_rows),all_brackets_resolved=all(r['signed_resolved_bracket'] for r in metric_rows),all_mpp_sampled_bounds=all(r['mpp_no_sampled_underperformance'] for r in metric_rows),all_occupations_valid=all(r['occupation_valid'] for r in eos_rows),max_eos_c_relative=max((r['max_c_relative_error'] for r in eos_rows),default=None),max_eos_dc_relative=max((r['max_dc_relative_error'] for r in eos_rows),default=None),max_eos_g_relative=max((r['max_g_relative_error'] for r in eos_rows),default=None),max_highorder_crosscheck_relative=max((max(r['max_2049_4097_c_relative_error'],r['max_2049_4097_dc_relative_error']) for r in eos_rows),default=None),max_direct_eos_current_shift_Acm2=max((r['highorder_eos_current_shift_Acm2'] for r in state_rows),default=None),max_direct_eos_local_continuity_L1_Acm2=max((r['highorder_eos_local_continuity_L1_Acm2'] for r in state_rows),default=None),frozen_input_comparison_count=len(frozen_rows),source_tree=str(args.source.relative_to(ROOT)) if args.source.is_relative_to(ROOT) else str(args.source),all_frozen_input_comparisons_exact=bool(frozen_rows) and all(all(r[k] for k in ['parameters_exact','gamma_n_exact','gamma_p_exact','field_cap_exact','Rs_zero']) for r in frozen_rows),caveats=['Closure-based current estimates exclude mesh, constitutive, contact, calibration and extrapolation uncertainties.','Near intrinsic open circuit, opposite carrier currents may exceed terminal current by many orders of magnitude; endpoint signs are tested against closure uncertainty.','Root voltage digits must also respect independent mesh shifts and direct-EOS perturbations; optimizer tolerance alone is not an uncertainty bound.'])
(OUT/(args.prefix+'_summary.json')).write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
