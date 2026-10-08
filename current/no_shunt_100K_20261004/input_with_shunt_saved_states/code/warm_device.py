"""Fixed-RT-anchor warm exploration. No per-temperature optimization.
Temperature nodes restricted to 225–300 K. Signed terminal J includes the frozen
room-temperature parallel leakage; contacts, gap and generation are unchanged.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
from dataclasses import asdict
import csv,json,time,traceback,hashlib
import numpy as np
from scipy.optimize import brentq,minimize_scalar
from device import Parameters,Device,solve,advance,poisson_seed,KB,Q
from field_device import FieldDevice
ROOT=Path(__file__).resolve().parents[1]
FIT=json.loads((ROOT/'reference/room_fit_frozen.json').read_text())
ANCHOR=FIT['parameters'];RSH=FIT['Rsh_ohmcm2']
AUTHOR={float(r['T_K']):float(r['k2_author_plot_rounded_cm3_s']) for r in csv.DictReader((ROOT/'reference/author_fit_lines.csv').open())}
INVERSE={float(r['T_K']):float(r['beta_local_cm3s']) for r in csv.DictReader((ROOT/'reference/observable_local_beta.csv').open())}
# Source-vector interpolation at target temperatures. Exact source table copied
# by evidence audit; values below are replaced from that CSV when present.
TP=np.array([225.,245.,267.,284.,306.]);GP=np.array([.00369813,.002993,.002125,.001629,.000964])
TN=np.array([225.,250.,275.,300.]);GN=np.array([.00157325,.00195068,.00207492,.00137675])

def source_gammas(T):
    f=ROOT/'reference/transport_at_target_temperatures.csv'
    if f.exists():
        rr=list(csv.DictReader(f.open()));r=next(r for r in rr if float(r['T_K'])==T)
        return float(r['gamma_n_relative_to_300']),float(r['gamma_p_relative_to_300'])
    raise FileNotFoundError('Exact audited transport source table is required')

def parameters(T,model):
    if T not in AUTHOR or not 225<=T<=300:raise ValueError('Only supported warm 225/250/275/300 K nodes are permitted')
    if model not in ['rt_author_hole','rt_author_frozen','rt_author_conditional','rt_inverse_hole','absolute_author_hole','absolute_author_frozen']:raise ValueError(model)
    beta=AUTHOR[T] if model.startswith('absolute') else ANCHOR['beta']*((INVERSE[T]/INVERSE[300.]) if model=='rt_inverse_hole' else AUTHOR[T]/AUTHOR[300.])
    rn,rp=source_gammas(T)
    gn=ANCHOR['gamma_n']*(rn if model=='rt_author_conditional' else 1.)
    gp=ANCHOR['gamma_p']*(1. if model.endswith('frozen') else rp)
    p=Parameters(T=T,G_cm3s=ANCHOR['G'],Eg_eV=ANCHOR['Eg'],beta_cm3s=beta,barrier_n_eV=ANCHOR['bn'],barrier_p_eV=ANCHOR['bp'])
    return p,gn,gp

def make_device(T,N,model,cap=None):
    p,gn,gp=parameters(T,model)
    return FieldDevice(p,N,gamma_n=gn,gamma_p=gp,field_cap=ANCHOR['field_cap'] if cap is None else cap)

def observable(d,z,V,L,zeq):
    r=d.ledger(z,V,L);o=d.evaluate(z,V);ref=d.evaluate(zeq,0);x=d.x
    integ=lambda y:float(np.trapezoid(y,x))
    dn=o['n']-ref['n'];dp=o['p']-ref['p'];js=V/RSH;j=r['J_Acm2']+js
    row=dict(**r,V_intrinsic_V=float(V),V_terminal_V=float(V),J_intrinsic_Acm2=r['J_Acm2'],J_shunt_Acm2=float(js),J_terminal_Acm2=float(j),external_shunt_heat_Wcm2=float(V*js),external_terminal_Wcm2=float(V*j))
    row['whole_circuit_energy_residual_Wcm2']=r['transport_dissipation_Wcm2']+r['recombination_dissipation_Wcm2']-r['light_chemical_work_Wcm2']+V*js-V*j
    row['terminal_sign_resolved']=bool(abs(j)>r['current_closure_uncertainty_Acm2'])
    core=(x>=.05)&(x<=.95);cent=(x>=.45)&(x<=.55)
    R=o['R'][1:-1];near=(x[1:-1]<.05)|(x[1:-1]>.95)
    row.update(n_bace_weighted_excess_cm3=integ((1-x)*dn+x*dp),n_pair_excess_mean_cm3=.5*integ(dn+dp),
        n_center_cm3=float(np.interp(.5,x,o['n'])),p_center_cm3=float(np.interp(.5,x,o['p'])),
        center_QFLS_eV=float(np.interp(.5,x,o['A'])*d.vt),volume_mean_QFLS_eV=float(np.dot(d.vol,o['A'][1:-1])*d.vt/sum(d.vol)),
        max_core_density_fraction=float(max(o['n'][core].max(),o['p'][core].max())/d.p.N0_cm3),
        min_local_reaction_dissipation_Wcm3=float(np.min(Q*o['R']*o['A']*d.vt)),
        near_contact_recombination_Acm2=float(Q*np.dot(d.vol[near],R[near])),
        source_minus_recombination_Acm2=float(r['generation_Acm2']-r['recombination_Acm2']),
        electron_J_left_Acm2=float(o['Jn'][0]),electron_J_right_Acm2=float(o['Jn'][-1]),hole_J_left_Acm2=float(o['Jp'][0]),hole_J_right_Acm2=float(o['Jp'][-1]),
        BACE_proxy_status='conditional static excess extraction-weighted proxy, not transient BACE',
        frozen_shunt_status='300 K fitted shunt held constant; unvalidated temperature law')
    return row

def run(T,N=321,model='rt_author_hole',step=.02,I=1.,metrics_only=False,cap=None):
    tag=f'{model}_T{T:g}_N{N}'+(f'_I{I:g}' if I!=1 else '')+('_metrics' if metrics_only else '')+(f'_cap{cap:g}' if cap else '')
    out=ROOT/'data'/f'{tag}.json'
    if out.exists():
        old=json.loads(out.read_text())
        if old.get('success'):
            want,gn,gp=parameters(T,model);expected_cap=ANCHOR['field_cap'] if cap is None else cap
            if old['parameters']!=asdict(want) or old['gamma_n']!=gn or old['gamma_p']!=gp or old['field_cap']!=expected_cap or old['Rsh_ohmcm2']!=RSH:
                raise RuntimeError('Cached result does not match frozen inputs; preserve it and recompute in a fresh data directory')
            print('CACHED',tag,flush=True);return old
    d=make_device(T,N,model,cap);p=d.p;start=time.monotonic();trace=[];curves={};states=[];labels=[];metrics={};critical={};cache={}
    result=dict(tag=tag,model=model,intensity=I,nodes=N,parameters=asdict(p),gamma_n=d.gamma_n,gamma_p=d.gamma_p,field_cap=d.field_cap,Rsh_ohmcm2=RSH,Rs_ohmcm2=0.,status='fixed-300K-anchor conditional temperature extrapolation within source warm range; no per-T retuning',curves=curves,metrics=metrics,critical_states=critical,success=False)
    def persist():
        result['seconds']=time.monotonic()-start;result['trace']=trace
        out.write_text(json.dumps(result,indent=2))
        if states:np.savez_compressed(out.with_suffix('.npz'),z=np.array(states),labels=np.array(labels),x_cm=d.x*p.d_cm)
    def save(z,V,L,category):
        row=observable(d,z,V,L,zeq);row['category']=category
        states.append(z.copy());labels.append([L,V]);cache[(float(V),float(L))]=(row,z.copy());return row
    try:
        zeq,info=solve(d,0,d.initial());trace.append({k:v for k,v in info.items() if k!='history'})
        if not info['success']:raise RuntimeError('Equilibrium '+info.get('reason',''))
        critical['equilibrium']=save(zeq,0,0,'equilibrium')
        nc=np.sqrt(I*p.G_cm3s/p.beta_cm3s);vseed=float(p.Eg_eV+d.vt*(d.en.eta_from_c(nc/p.N0_cm3)+d.ep.eta_from_c(nc/p.N0_cm3)))
        zlight=poisson_seed(d,vseed);zlight,info=solve(d,vseed,zlight,I,maxiter=240)
        if not info['success']:raise RuntimeError('Light seed '+info.get('reason',''))
        save(zlight,vseed,I,'seed')
        for L in ([I] if metrics_only else [0.,I]):
            arr=[];curves['dark' if L==0 else 'light']=arr
            vv=np.array([0.,.4,.6,.7,.75,.8,.85,.9,.95,1.,1.1]) if metrics_only else np.round(np.arange(0,1.100001,step),9)
            if L:branches=[(zlight.copy(),vseed,sorted([v for v in vv if v<vseed],reverse=True)),(zlight.copy(),vseed,sorted([v for v in vv if v>=vseed]))]
            else:branches=[(zeq.copy(),0.,list(vv))]
            for z,old,targets in branches:
                for V in targets:
                    if V!=old:z,info=advance(d,z,float(old),float(V),L,L,trace)
                    old=float(V);row=save(z,float(V),L,'sweep');arr.append(row)
            arr.sort(key=lambda r:r['V'])
        def evaluate(V,L=I):
            k=(float(V),float(L))
            if k in cache:return cache[k]
            candidates=[(abs(v-V),v,z) for (v,ell),(r,z) in cache.items() if ell==L];_,v0,z0=min(candidates,key=lambda x:x[0])
            z,info=advance(d,z0,v0,float(V),L,L,trace);r=save(z,float(V),L,'metric');return r,z
        light=curves['light'];cross=[(a,b) for a,b in zip(light[:-1],light[1:]) if a['J_terminal_Acm2']<0<b['J_terminal_Acm2']]
        if not cross:raise RuntimeError('No signed terminal light root within 0–1.1 V')
        a,b=cross[0];vc=float(brentq(lambda v:evaluate(v)[0]['J_terminal_Acm2'],a['V'],b['V'],xtol=2e-9))
        oc,zoc=evaluate(vc);sc,zsc=evaluate(0.)
        opt=minimize_scalar(lambda v:evaluate(v)[0]['external_terminal_Wcm2'],bounds=(0,vc),method='bounded',options={'xatol':2e-7})
        mp,zmp=evaluate(float(opt.x));br=[evaluate(vc-5e-6)[0],evaluate(vc+5e-6)[0]]
        critical.update(sc=sc,oc=oc,mpp=mp,bracket=br)
        metrics.update(T_K=T,model=model,N=N,intensity=I,Jsc_mAcm2=-sc['J_terminal_Acm2']*1000,Voc_V=vc,FF=mp['external_terminal_Wcm2']/(sc['J_terminal_Acm2']*vc),Pmax_mWcm2=-mp['external_terminal_Wcm2']*1000,Vmp_V=float(opt.x),Jmp_mAcm2=-mp['J_terminal_Acm2']*1000,
            beta_cm3s=p.beta_cm3s,gamma_n=d.gamma_n,gamma_p=d.gamma_p,mu_n_zero_cm2Vs=d.mn,mu_p_zero_cm2Vs=d.mp,
            n_oc_bace_proxy_cm3=oc['n_bace_weighted_excess_cm3'],n_oc_center_cm3=oc['n_center_cm3'],p_oc_center_cm3=oc['p_center_cm3'],QFLS_center_oc_eV=oc['center_QFLS_eV'],
            Voc_resolved_10uV_bracket=bool(br[0]['J_terminal_Acm2']<0<br[1]['J_terminal_Acm2'] and all(r['terminal_sign_resolved'] and r['gate_passed'] for r in br)))
        profiles=[];ref=d.evaluate(zeq,0.)
        for name,vv,zz in [('sc',0.,zsc),('oc',vc,zoc),('mpp',float(opt.x),zmp)]:
            o=d.evaluate(zz,vv)
            for k,x in enumerate(d.x):profiles.append(dict(condition=name,T_K=T,V_V=vv,x_cm=x*p.d_cm,x_fraction=x,n_cm3=o['n'][k],p_cm3=o['p'][k],n_dark0_cm3=ref['n'][k],p_dark0_cm3=ref['p'][k],G_cm3s=p.G_cm3s*I,R_cm3s=o['R'][k],QFLS_eV=o['A'][k]*d.vt,phi_V=zz[k,0]*d.vt))
        with out.with_name(tag+'_profiles.csv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(profiles[0]));w.writeheader();w.writerows(profiles)
        allrows=[r for arr in curves.values() for r in arr]+[r for r,z in cache.values()]
        result['audit_summary']=dict(all_intrinsic_gates=all(r['gate_passed'] for r in allrows),max_continuity_Acm2=max(r['continuity_abs_Acm2'] for r in allrows),max_charge_relative=max(r['charge_relative'] for r in allrows),max_whole_circuit_energy_residual_Wcm2=max(abs(r['whole_circuit_energy_residual_Wcm2']) for r in allrows),min_local_reaction_dissipation_Wcm3=min(r['min_local_reaction_dissipation_Wcm3'] for r in allrows),max_core_density_fraction=max(r['max_core_density_fraction'] for r in allrows))
        result['success']=True
    except Exception as e:result.update(error=str(e),traceback=traceback.format_exc())
    persist();print('DONE',tag,'success',result['success'],metrics,flush=True);return result

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--T',type=float,required=True);ap.add_argument('--N',type=int,default=321);ap.add_argument('--model',default='rt_author_hole');ap.add_argument('--step',type=float,default=.02);ap.add_argument('--I',type=float,default=1.);ap.add_argument('--metrics-only',action='store_true');ap.add_argument('--cap',type=float);a=ap.parse_args();run(**vars(a))
