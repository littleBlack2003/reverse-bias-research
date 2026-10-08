"""Latest frozen RT calibration, conditional 100–300 K extension.
No per-temperature inverse calibration. All old core constitutive code unchanged.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
from dataclasses import asdict
import csv,json,time,traceback
import numpy as np
from scipy.optimize import brentq,minimize_scalar
from device import Parameters,solve,advance,poisson_seed,KB
from field_device import FieldDevice
from warm_device import ANCHOR,RSH,AUTHOR,source_gammas,observable
ROOT=Path(__file__).resolve().parents[1]
WARM=ROOT.parent/'latest_warm_source_20261004/extracted/calibrated_warm_exploration_20261004'
BETA_EA=KB*np.log(AUTHOR[225.]/AUTHOR[200.])/(1/200.-1/225.)
GP225=ANCHOR['gamma_p']*source_gammas(225.)[1]
GP250=ANCHOR['gamma_p']*source_gammas(250.)[1]
GP_SLOPE=(GP225-GP250)/(1/225.**2-1/250.**2)

def parameters(T):
    if not 100<=T<=300:raise ValueError('Declared conditional extension covers 100–300 K only')
    tt=np.array(sorted(AUTHOR))
    if T>=200: beta=ANCHOR['beta']*np.exp(np.interp(1/T,1/tt[::-1],np.log([AUTHOR[x] for x in tt[::-1]])))/AUTHOR[300.]
    else: beta=ANCHOR['beta']*AUTHOR[200.]/AUTHOR[300.]*np.exp(-BETA_EA/KB*(1/T-1/200.))
    if T>=225:
        t=np.array([225.,250.,275.,300.]); gp=np.interp(T,t,[ANCHOR['gamma_p']*source_gammas(x)[1] for x in t])
    else: gp=GP225+GP_SLOPE*(1/T**2-1/225.**2)
    # Preserve original warm source values at exact nodes to floating-point precision.
    if T in [225.,250.,275.,300.]: beta=ANCHOR['beta']*AUTHOR[T]/AUTHOR[300.]
    p=Parameters(T=T,G_cm3s=ANCHOR['G'],Eg_eV=ANCHOR['Eg'],beta_cm3s=float(beta),barrier_n_eV=ANCHOR['bn'],barrier_p_eV=ANCHOR['bp'],allow_mobility_extrapolation=T<223)
    return p,ANCHOR['gamma_n'],float(gp)

def make_device(T,N=321,cap=None):
    p,gn,gp=parameters(float(T)); return FieldDevice(p,N,gamma_n=gn,gamma_p=gp,field_cap=ANCHOR['field_cap'] if cap is None else cap)

def temperature_continue(z,T0,T1,V,L,N,trace,cap=None,depth=0):
    d=make_device(T1,N,cap);guess=z*(T0/T1)
    zz,info=solve(d,V,guess,L,maxiter=220)
    trace.append(dict(stage='temperature_continuation',from_T=T0,to_T=T1,depth=depth,**{k:v for k,v in info.items() if k!='history'}))
    if info['success']:return zz
    if depth>=8:raise RuntimeError(f'Temperature continuation failed {T0}->{T1} K at V={V}: {info.get("reason")}')
    tm=(T0+T1)/2;mid=temperature_continue(z,T0,tm,V,L,N,trace,cap,depth+1)
    return temperature_continue(mid,tm,T1,V,L,N,trace,cap,depth+1)

def run(T,N=321,previous=None,step=.02,cap=None):
    d=make_device(T,N,cap);p=d.p;tag=f'latest_T{T:g}_N{N}'+('' if cap is None else f'_cap{cap:g}');out=ROOT/'data'/f'{tag}.json'
    if out.exists():
        old=json.loads(out.read_text())
        if old.get('success'):
            if old['parameters']!=asdict(p) or old['gamma_n']!=d.gamma_n or old['gamma_p']!=d.gamma_p or old['field_cap']!=d.field_cap or old['Rsh_ohmcm2']!=RSH or old['beta_extension_Ea_eV']!=BETA_EA or old['gamma_p_extension_slope']!=GP_SLOPE:
                raise RuntimeError('Cached result does not match declared frozen inputs; use a fresh data directory')
            print('CACHED',tag,flush=True);return old
    start=time.monotonic();trace=[];states=[];labels=[];cache={};curves={'light':[]};metrics={};critical={}
    result=dict(tag=tag,model='latest_rt_author_hole_conditional_lowT',nodes=N,parameters=asdict(p),gamma_n=d.gamma_n,gamma_p=d.gamma_p,field_cap=d.field_cap,Rsh_ohmcm2=RSH,Rs_ohmcm2=0.,beta_extension_Ea_eV=BETA_EA,gamma_p_extension_slope=GP_SLOPE,curves=curves,metrics=metrics,critical_states=critical,success=False,status='conditional extrapolation below mobility/gamma source range; fixed latest 300 K calibration; no per-T refit')
    def persist():
        result['seconds']=time.monotonic()-start;result['trace']=trace;out.write_text(json.dumps(result,indent=2))
        if states:np.savez_compressed(out.with_suffix('.npz'),z=np.array(states),labels=np.array(labels),x_cm=d.x*p.d_cm)
    def save(z,V,L,category):
        row=observable(d,z,float(V),float(L),zeq);row['category']=category
        r=d.residual(z,V,L);jac=d.jacobian(z,V,L);scale=np.maximum(np.asarray(abs(jac).max(axis=1).toarray()).ravel(),1e-250)
        row['scaled_nonlinear_residual']=float(max(abs(r)/scale))
        oo=d.evaluate(z,V);row['min_density_fraction']=float(min(oo['n'].min(),oo['p'].min())/p.N0_cm3)
        row['max_generalized_Einstein_factor']=float(max(oo['g_n'].max(),oo['g_p'].max()))
        states.append(z.copy());labels.append([L,V]);cache[(float(V),float(L))]=(row,z.copy());return row
    try:
        zeq,info=solve(d,0,d.initial());trace.append(dict(stage='equilibrium',**{k:v for k,v in info.items() if k!='history'}))
        if not info['success']:raise RuntimeError('Equilibrium '+info.get('reason',''))
        critical['equilibrium']=save(zeq,0,0,'equilibrium')
        nc=np.sqrt(p.G_cm3s/p.beta_cm3s);vseed=float(p.Eg_eV+d.vt*(d.en.eta_from_c(nc/p.N0_cm3)+d.ep.eta_from_c(nc/p.N0_cm3)))
        if previous:
            fn=ROOT/'data'/f"{previous['tag']}.npz";ss=np.load(fn);ii=np.where(ss['labels'][:,0]==1)[0];i=ii[np.argmin(abs(ss['labels'][ii,1]-vseed))];vold=float(ss['labels'][i,1])
            zlight=temperature_continue(ss['z'][i],previous['parameters']['T'],T,vold,1,N,trace,cap)
            zlight,info=advance(d,zlight,vold,vseed,1,1,trace)
        else:
            zlight=poisson_seed(d,vseed);zlight,info=solve(d,vseed,zlight,1,maxiter=240)
            trace.append(dict(stage='light_seed',**{k:v for k,v in info.items() if k!='history'}))
            if not info['success']:raise RuntimeError('Light seed '+info.get('reason',''))
        save(zlight,vseed,1,'seed');persist()
        def evaluate(V):
            k=(float(V),1.)
            if k in cache:return cache[k]
            cand=[(abs(v-V),v,z) for (v,L),(r,z) in cache.items() if L==1.];_,v0,z0=min(cand,key=lambda a:a[0])
            z,info=advance(d,z0,v0,float(V),1,1,trace);r=save(z,V,1,'evaluation');return r,z
        # Follow a solved illuminated seed down to short circuit, retaining every requested bias.
        for V in sorted(np.arange(0,max(vseed,step),step),reverse=True):
            row,z=evaluate(float(round(V,9)));curves['light'].append(row);persist()
        sc,zsc=evaluate(0.)
        # Establish a signed positive bracket even if the physical open circuit moves.
        seedrow,_=evaluate(vseed);vhi=vseed
        while seedrow['J_terminal_Acm2']<=0 and vhi<1.5:
            vhi+=step;seedrow,_=evaluate(vhi)
        if sc['J_terminal_Acm2']>=0 or seedrow['J_terminal_Acm2']<=0:raise RuntimeError('No positive-power signed terminal root')
        arr=sorted([r for (v,L),(r,z) in cache.items() if L==1],key=lambda r:r['V'])
        a,b=next((a,b) for a,b in zip(arr[:-1],arr[1:]) if a['J_terminal_Acm2']<0<b['J_terminal_Acm2'])
        vc=float(brentq(lambda v:evaluate(v)[0]['J_terminal_Acm2'],a['V'],b['V'],xtol=2e-11))
        oc,zoc=evaluate(vc)
        # Add a linear voltage grid resolving the shrinking low-T power quadrant.
        for v in np.linspace(0,vc+min(.04,.2*vc+.002),61):evaluate(float(v))
        opt=minimize_scalar(lambda v:evaluate(v)[0]['external_terminal_Wcm2'],bounds=(0,vc),method='bounded',options={'xatol':max(2e-10,vc*1e-7)})
        mp,zmp=evaluate(float(opt.x));eps=min(5e-6,vc*1e-3);br=[evaluate(vc-eps)[0],evaluate(vc+eps)[0]]
        critical.update(sc=sc,oc=oc,mpp=mp,bracket=br)
        metrics.update(T_K=T,N=N,Jsc_mAcm2=-sc['J_terminal_Acm2']*1000,Voc_V=vc,FF=mp['external_terminal_Wcm2']/(sc['J_terminal_Acm2']*vc),Pmax_mWcm2=-mp['external_terminal_Wcm2']*1000,Vmp_V=float(opt.x),Jmp_mAcm2=-mp['J_terminal_Acm2']*1000,beta_cm3s=p.beta_cm3s,gamma_n=d.gamma_n,gamma_p=d.gamma_p,mu_n_zero_cm2Vs=d.mn,mu_p_zero_cm2Vs=d.mp,n_oc_center_cm3=oc['n_center_cm3'],p_oc_center_cm3=oc['p_center_cm3'],Voc_resolved_bracket=bool(br[0]['J_terminal_Acm2']<0<br[1]['J_terminal_Acm2'] and all(r['terminal_sign_resolved'] and r['gate_passed'] for r in br)),Voc_bracket_half_width_V=eps)
        maxV=vc+min(.04,.2*vc+.002)
        curves['light']=sorted([r for (v,L),(r,z) in cache.items() if L==1 and 0<=v<=maxV+1e-10],key=lambda r:r['V'])
        allrows=[r for r,z in cache.values()]
        result['audit_summary']=dict(all_intrinsic_gates=all(r['gate_passed'] for r in allrows),max_scaled_nonlinear_residual=max(r['scaled_nonlinear_residual'] for r in allrows),max_continuity_Acm2=max(r['continuity_abs_Acm2'] for r in allrows),max_current_closure_uncertainty_Acm2=max(r['current_closure_uncertainty_Acm2'] for r in allrows),max_charge_relative=max(r['charge_relative'] for r in allrows),max_whole_circuit_energy_residual_Wcm2=max(abs(r['whole_circuit_energy_residual_Wcm2']) for r in allrows),max_core_density_fraction=max(r['max_core_density_fraction'] for r in allrows),max_density_fraction=max(r['max_density_fraction'] for r in allrows),max_field_capped_length_fraction=max(r['field_capped_length_fraction'] for r in allrows),max_field_Vcm=max(r['field_max_Vcm'] for r in allrows),max_mobility_factor_p=max(r['mobility_factor_p_max'] for r in allrows),max_mobility_factor_n=max(r['mobility_factor_n_max'] for r in allrows))
        profiles=[]
        for name,v,z in [('sc',0,zsc),('oc',vc,zoc),('mpp',opt.x,zmp)]:
            o=d.evaluate(z,v)
            for i,x in enumerate(d.x):profiles.append(dict(condition=name,T_K=T,V_V=v,x_cm=x*p.d_cm,x_fraction=x,n_cm3=o['n'][i],p_cm3=o['p'][i],G_cm3s=p.G_cm3s,R_cm3s=o['R'][i],QFLS_eV=o['A'][i]*d.vt,phi_V=z[i,0]*d.vt))
        with out.with_name(tag+'_profiles.csv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(profiles[0]));w.writeheader();w.writerows(profiles)
        result['success']=True
    except Exception as e:result.update(error=str(e),traceback=traceback.format_exc())
    persist();print('DONE',tag,'success',result['success'],metrics,result.get('error'),flush=True);return result

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--temperatures',nargs='+',type=float,default=[300,275,250,225,200,175,150,125,100]);ap.add_argument('--N',type=int,default=321);ap.add_argument('--cap',type=float);ap.add_argument('--independent',action='store_true');a=ap.parse_args()
    prev=None
    for T in a.temperatures:
        res=run(T,N=a.N,previous=None if a.independent else prev,cap=a.cap)
        if not res['success']:break
        prev=res
