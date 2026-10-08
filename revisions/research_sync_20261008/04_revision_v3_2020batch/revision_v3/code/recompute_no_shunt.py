"""Recompute no-shunt JV, reusing exact same-bias PDE states and solving all new biases.
No constitutive or physical parameters are fitted or modified. Root calls solve PDE.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import json,time,csv,traceback,hashlib
from pathlib import Path
from dataclasses import asdict
import numpy as np
from scipy.optimize import brentq,minimize_scalar
from device import solve,advance
from extend_device import make_device,BETA_EA,GP_SLOPE
from warm_device import observable
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'input_with_shunt_saved_states'
if not SOURCE.exists():SOURCE=ROOT.parent/'latest_calibrated_100K_20261004'

def run(T,N=321,metrics_only=False):
    tag=f'no_shunt_T{T:g}_N{N}';out=ROOT/'data'/f'{tag}.json'
    if out.exists():
        old=json.loads(out.read_text())
        if old.get('success'):
            d=make_device(T,N)
            assert old['parameters']==asdict(d.p) and old['gamma_n']==d.gamma_n and old['gamma_p']==d.gamma_p and old['field_cap']==d.field_cap and old['external_shunt_conductance_Scm2']==0
            print('CACHED',tag,flush=True);return old
    d=make_device(T,N);p=d.p;start=time.monotonic();trace=[];cache={};critical={};metrics={};counts={'reused_same_grid':0,'interpolated_initial_guess':0,'new_solved_states':0}
    result=dict(tag=tag,model='latest_frozen_parameters_no_external_shunt',nodes=N,parameters=asdict(p),gamma_n=d.gamma_n,gamma_p=d.gamma_p,field_cap=d.field_cap,Rsh_ohmcm2=None,external_shunt_enabled=False,external_shunt_conductance_Scm2=0.,Rs_ohmcm2=0.,beta_extension_Ea_eV=BETA_EA,gamma_p_extension_slope=GP_SLOPE,curves={'light':[]},metrics=metrics,critical_states=critical,counts=counts,success=False,status='Conditional temperature extrapolation, external shunt removed, no parameter refit')
    def persist():
        result['seconds']=time.monotonic()-start;result['trace']=trace
        result['curves']['light']=sorted([r for (v,l),(r,z) in cache.items() if l==1],key=lambda r:r['V'])
        out.write_text(json.dumps(result,indent=2))
        if cache:np.savez_compressed(out.with_suffix('.npz'),z=np.array([z for r,z in cache.values()]),labels=np.array([[l,v] for v,l in cache]),x_cm=d.x*p.d_cm)
    def save(z,V,L,category):
        row=observable(d,z,float(V),float(L),zeq);row['category']=category
        r=d.residual(z,V,L);jac=d.jacobian(z,V,L);scale=np.maximum(np.asarray(abs(jac).max(axis=1).toarray()).ravel(),1e-250)
        row['scaled_nonlinear_residual']=float(max(abs(r)/scale))
        o=d.evaluate(z,V);row['min_density_fraction']=float(min(o['n'].min(),o['p'].min())/p.N0_cm3)
        row['max_generalized_Einstein_factor']=float(max(o['g_n'].max(),o['g_p'].max()))
        row['min_eta_n']=float(o['eta_n'].min());row['max_eta_n']=float(o['eta_n'].max());row['min_eta_p']=float(o['eta_p'].min());row['max_eta_p']=float(o['eta_p'].max())
        if not row['gate_passed'] or row['scaled_nonlinear_residual']>=3e-10:raise RuntimeError(f'State fails gates at {V}: {row["gate_passed"]}, {row["scaled_nonlinear_residual"]}')
        cache[(float(V),float(L))]=(row,z.copy());return row
    try:
        oldfile=SOURCE/'data'/f'latest_T{T:g}_N{N}.npz'
        oldjson=oldfile.with_suffix('.json')
        if oldfile.exists():
            old=json.loads(oldjson.read_text());assert old['parameters']==asdict(p) and old['gamma_n']==d.gamma_n and old['gamma_p']==d.gamma_p and old['field_cap']==d.field_cap and old['Rs_ohmcm2']==0
            ss=np.load(oldfile);ii=np.where(ss['labels'][:,0]==0)[0];assert len(ii)
            zeq=ss['z'][ii[np.argmin(abs(ss['labels'][ii,1]))]].copy()
            for z,(L,V) in zip(ss['z'],ss['labels']):save(z,V,L,'reused_same_bias_same_PDE');counts['reused_same_grid']+=1
            result['reuse_source']=str(oldfile.relative_to(ROOT.parent));result['reuse_source_sha256']=hashlib.sha256(oldfile.read_bytes()).hexdigest()
        else:
            # Mesh refinement starts from solved coarse profiles, never declares interpolation a solution.
            base=np.load(ROOT/'data'/f'no_shunt_T{T:g}_N321.npz');xb=base['x_cm']/p.d_cm
            eqi=np.where(base['labels'][:,0]==0)[0][0];guess=np.column_stack([np.interp(d.x,xb,base['z'][eqi,:,c]) for c in range(3)])
            zeq,info=solve(d,0,guess,0,maxiter=240);trace.append({k:v for k,v in info.items() if k!='history'})
            if not info['success']:raise RuntimeError('Refined equilibrium failed '+info.get('reason',''))
            save(zeq,0,0,'new_refined_equilibrium');counts['new_solved_states']+=1
            row0=json.loads((ROOT/'data'/f'no_shunt_T{T:g}_N321.json').read_text());target=row0['metrics']['Voc_V']
            ii=np.where(base['labels'][:,0]==1)[0];i=ii[np.argmin(abs(base['labels'][ii,1]-target))];V=float(base['labels'][i,1]);guess=np.column_stack([np.interp(d.x,xb,base['z'][i,:,c]) for c in range(3)])
            z,info=solve(d,V,guess,1,maxiter=240);trace.append({k:v for k,v in info.items() if k!='history'})
            if not info['success']:raise RuntimeError('Refined illuminated seed failed '+info.get('reason',''))
            save(z,V,1,'new_refined_illuminated_seed');counts['new_solved_states']+=1;counts['interpolated_initial_guess']+=2
        critical['equilibrium']=cache[(0.,0.)][0]
        def evaluate(V):
            V=float(V);key=(V,1.)
            if key in cache:return cache[key]
            candidates=[(abs(v-V),v,z) for (v,L),(r,z) in cache.items() if L==1.];_,v0,z0=min(candidates,key=lambda a:a[0])
            z,info=advance(d,z0,v0,V,1,1,trace);r=save(z,V,1,'new_selfconsistent_solve');counts['new_solved_states']+=1
            if counts['new_solved_states']%10==0:persist()
            return r,z
        sc,zsc=evaluate(0.)
        if sc['J_terminal_Acm2']>=0:raise RuntimeError('Non-photovoltaic short circuit')
        arr=sorted([r for (v,L),(r,z) in cache.items() if L==1.],key=lambda r:r['V'])
        vhi=arr[-1]['V'];hi=arr[-1]
        while hi['J_terminal_Acm2']<=0:
            vhi+=.02;hi,_=evaluate(vhi)
            if vhi>2:raise RuntimeError('No signed root by2V; inspect mechanism before widening')
        arr=sorted([r for (v,L),(r,z) in cache.items() if L==1.],key=lambda r:r['V'])
        a,b=next((a,b) for a,b in zip(arr[:-1],arr[1:]) if a['J_terminal_Acm2']<0<b['J_terminal_Acm2'])
        vc=float(brentq(lambda v:evaluate(v)[0]['J_terminal_Acm2'],a['V'],b['V'],xtol=2e-10))
        oc,zoc=evaluate(vc);eps=5e-6;br=[evaluate(vc-eps)[0],evaluate(vc+eps)[0]]
        if not (br[0]['J_terminal_Acm2']<0<br[1]['J_terminal_Acm2'] and all(r['terminal_sign_resolved'] for r in br)):raise RuntimeError('Root bracket not sign resolved')
        if not metrics_only:
            # Uniform linear voltage coverage through the true intrinsic Voc and 40mV beyond.
            for v in np.linspace(0,vc+.04,121):evaluate(float(v))
        else:evaluate(vc+.04)
        opt=minimize_scalar(lambda v:evaluate(v)[0]['external_terminal_Wcm2'],bounds=(0,vc),method='bounded',options={'xatol':2e-8})
        mp,zmp=evaluate(float(opt.x));critical.update(sc=sc,oc=oc,mpp=mp,bracket=br,beyond_oc=evaluate(vc+.04)[0])
        metrics.update(T_K=T,N=N,Jsc_mAcm2=-sc['J_terminal_Acm2']*1000,Voc_V=vc,FF=mp['external_terminal_Wcm2']/(sc['J_terminal_Acm2']*vc),Pmax_mWcm2=-mp['external_terminal_Wcm2']*1000,Vmp_V=float(opt.x),Jmp_mAcm2=-mp['J_terminal_Acm2']*1000,beta_cm3s=p.beta_cm3s,gamma_n=d.gamma_n,gamma_p=d.gamma_p,mu_n_zero_cm2Vs=d.mn,mu_p_zero_cm2Vs=d.mp,n_oc_center_cm3=oc['n_center_cm3'],p_oc_center_cm3=oc['p_center_cm3'],Voc_resolved_bracket=True,Voc_bracket_half_width_V=eps,Voc_bracket_low_V=vc-eps,Voc_bracket_high_V=vc+eps,minimum_bracket_signal_to_closure=min(abs(r['J_terminal_Acm2'])/r['current_closure_uncertainty_Acm2'] for r in br))
        allrows=[r for r,z in cache.values()]
        result['audit_summary']={f'max_{key}':max(r[key] for r in allrows) for key in ['scaled_nonlinear_residual','continuity_abs_Acm2','current_closure_uncertainty_Acm2','charge_relative','max_density_fraction','max_core_density_fraction','field_capped_length_fraction','max_generalized_Einstein_factor']}
        result['audit_summary'].update(all_intrinsic_gates=all(r['gate_passed'] for r in allrows),max_abs_whole_circuit_energy_residual_Wcm2=max(abs(r['whole_circuit_energy_residual_Wcm2']) for r in allrows),external_shunt_identically_zero=all(r['J_shunt_Acm2']==0 and r['external_shunt_heat_Wcm2']==0 for r in allrows))
        profiles=[]
        for name,v,z in [('sc',0,zsc),('oc',vc,zoc),('mpp',opt.x,zmp)]:
            o=d.evaluate(z,v)
            for i,x in enumerate(d.x):profiles.append(dict(condition=name,T_K=T,V_V=v,x_cm=x*p.d_cm,x_fraction=x,n_cm3=o['n'][i],p_cm3=o['p'][i],G_cm3s=p.G_cm3s,R_cm3s=o['R'][i],QFLS_eV=o['A'][i]*d.vt,phi_V=z[i,0]*d.vt))
        with out.with_name(tag+'_profiles.csv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(profiles[0]));w.writeheader();w.writerows(profiles)
        result['success']=True
    except Exception as e:result.update(error=str(e),traceback=traceback.format_exc())
    persist();print('DONE',tag,'success',result['success'],'seconds',result['seconds'],'metrics',metrics,'error',result.get('error'),flush=True);return result

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--temperatures',nargs='+',type=float,default=list(range(100,301,25)));ap.add_argument('--N',type=int,default=321);ap.add_argument('--metrics-only',action='store_true');a=ap.parse_args()
    for T in a.temperatures:
        r=run(T,a.N,a.metrics_only)
        if not r['success']:raise SystemExit(1)
