"""Independent, bounded scientific review; read frozen source, write review only."""
from pathlib import Path
import sys, json, hashlib, itertools, platform
import numpy as np
import scipy
import mpmath as mp
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'snapshot'))
import spatial_model as sm

def rel(a,b):
    return float(abs(a-b)/max(abs(a),abs(b),mp.mpf('1e-200')))

def energy_independent(st,g,F,gauge=0,epsr=3.5):
    # Addition energies integrated from the neutral reference; pair interaction
    # expanded into core-core, core-electron and electron-electron terms.
    m=lambda x:mp.mpf(str(x))
    n=len(st);value=mp.mpf(0)
    for i in range(n):
        value+=m(g['eps'][i])*(st[i]-g['b'][i])
        value+=(g['b'][i]-st[i])*(m('.1')*m(F)*m(g['x'][i])+m(gauge))
    if epsr is not None:
        for i in range(n):
            for j in range(i+1,n):
                product=g['b'][i]*g['b'][j]-g['b'][i]*st[j]-st[i]*g['b'][j]+st[i]*st[j]
                value+=m(sm.A0)/m(epsr)*product/abs(m(g['x'][i])-m(g['x'][j]))
    return value

def check_state(result,F,T,delta_mu,gauge=0,epsr=3.5):
    row,p,Q,rates,details,E,mu,g,states,edges=result
    beta=1/(mp.mpf(str(sm.KB))*T)
    independent_E=[energy_independent(st,g,F,gauge,epsr) for st in states]
    ledger={}
    # Check detailed balance, number and electric dipole from state differences,
    # not the implementation's stored dx, tdg or mean charge fields.
    max_ldberr=mp.mpf(0);dN=mp.mpf(0);dp=mp.mpf(0);qtotal=mp.mpf(0)
    left=mp.mpf(0);right=mp.mpf(0);reservoir_shift=mp.mpf(0)
    for edge,r in zip(edges,rates):
        i,j,ty,k=edge
        dn=[states[j][z]-states[i][z] for z in range(len(g['x']))]
        nchange=sum(dn)
        a=independent_E[j]-independent_E[i]
        ldb=-beta*(a-(mu[k]*nchange if ty=='bath' else 0))
        max_ldberr=max(max_ldberr,abs(mp.log(r['kf']/r['kr'])-ldb))
        flux=p[i]*r['kf']-p[j]*r['kr']
        dN+=flux*nchange
        dp+=flux*sum(mp.mpf(str(g['x'][z]))*dn[z] for z in range(len(dn)))
        if ty=='bath':
            reservoir_shift-=flux*mp.mpf(str(g['xb'][k]))*nchange
            if k==0:left+=flux*nchange
            else:right+=flux*nchange
        qtotal-=flux*(a-(mu[k]*nchange if ty=='bath' else 0))
    R=left
    return {
        'max_energy_error_eV':float(max(abs(a-b) for a,b in zip(E,independent_E))),
        'max_edge_LDB_log_error':float(max_ldberr),
        'number_conservation_abs_s':float(abs(dN)),
        'dipole_stationarity_abs_nm_s':float(abs(dp)),
        'left_right_current_relative':rel(left,-right),
        'reported_current_relative':rel(R,mp.mpf(str(row['R_s']))),
        'heat_chemical_relative':rel(qtotal,mp.mpf(str(delta_mu))*R) if delta_mu else float(abs(qtotal)),
        'full_displacement_relative':rel(dp+reservoir_shift,mp.mpf(str(g['xb'][1]-g['xb'][0]))*R),
        'actual_drop_V':float(mp.mpf('.1')*mp.mpf(str(F))*mp.mpf(str(g['xb'][1]-g['xb'][0]))),
        'reported_drop_V':row['total_domain_potential_drop_V'],
        'R_s':float(R),'minimum_probability':float(min(p)),
        'charge_from_population_e':float(sum(p[i]*(sum(g['b'])-sum(st)) for i,st in enumerate(states))),
    }

def schur_and_lumping(result):
    row,p,Q,rates,details,E,mu,g,states,edges=result
    # Retain reference outer tag configuration Vo=1, Co=0, eight microstates.
    keep=[i for i,st in enumerate(states) if st[0]==1 and st[-1]==0]
    eliminate=[i for i in range(len(states)) if i not in keep]
    block=lambda ii,jj:mp.matrix([[Q[i,j] for j in jj] for i in ii])
    AA=block(keep,keep);AB=block(keep,eliminate);BA=block(eliminate,keep);BB=block(eliminate,eliminate)
    lift=-(BB**-1)*BA
    eff=AA+AB*lift
    w=mp.matrix([1+sum(lift[j,i] for j in range(lift.rows)) for i in range(lift.cols)])
    M=eff.copy();rhs=mp.matrix([0]*len(keep));rhs[len(keep)-1]=1
    for i in range(len(keep)):M[len(keep)-1,i]=w[i]
    pa=mp.lu_solve(M,rhs);pb=lift*pa
    restored=mp.matrix([0]*len(states))
    for i,ii in enumerate(keep):restored[ii]=pa[i]
    for i,ii in enumerate(eliminate):restored[ii]=pb[i]
    flux=sum(restored[i]*r['kf']-restored[j]*r['kr'] for (i,j,ty,k),r in zip(edges,rates) if ty=='bath' and k==0)
    R=sum(p[i]*r['kf']-p[j]*r['kr'] for (i,j,ty,k),r in zip(edges,rates) if ty=='bath' and k==0)
    # Explicit strong-lumpability counterexample: same VDC=100, Vo empty/full.
    state_index={s:i for i,s in enumerate(states)}
    a=state_index[(0,1,0,0,0)];b=state_index[(1,1,0,0,0)]
    dest=[i for i,st in enumerate(states) if st[1:4]==(0,1,0)]
    wa=sum(Q[i,a] for i in dest);wb=sum(Q[i,b] for i in dest)
    return {'retained_probability':float(sum(pa)),
            'naively_renormalized_current_multiplier':float(1/sum(pa)),
            'full_probability_max_relative':float(max(abs(restored[i]/p[i]-1) for i in range(len(p)))),
            'restored_current_relative':rel(flux,R),
            'schur_column_sum_abs':float(max(abs(sum(eff[:,j])) for j in range(eff.cols))),
            'strong_lumpability_example':{'source_microstates':['01000','11000'],'source_inner':'100','destination_inner':'010','rates_s':[float(wa),float(wb)],'ratio':float(wa/wb)}}

def main():
    results={};solves={}
    cases=[('colocated',{}),('direct',{}),('explicit',{}),('explicit_gauge',{'gauge':.137}),('explicit_common_mu',{'delta_mu':0.,'mean_mu':.023}),('direct_no_Coulomb',{'epsr':None}),('explicit_no_Coulomb',{'epsr':None})]
    for label,kwargs in cases:
        variant=label.split('_')[0]
        sol=sm.solve(1.5,80,variant=variant,detail=True,**kwargs)
        solves[label]=sol
        results[label]=check_state(sol,1.5,80,kwargs.get('delta_mu',1.5),kwargs.get('gauge',0),kwargs.get('epsr',3.5))
        print(label,results[label]['R_s'],flush=True)
    a=solves['explicit'];b=solves['explicit_gauge']
    results['gauge_invariance']={'population_max_relative':float(max(abs(x/y-1) for x,y in zip(a[1],b[1]))),'current_relative':rel(a[0]['R_s'],b[0]['R_s']),'phonon_heat_relative':rel(a[0]['phonon_heat_eV_s'],b[0]['phonon_heat_eV_s'])}
    results['schur_and_lumping']=schur_and_lumping(a)
    # Analytic final-exchange mismatch for isolated outer pair (Vo hole, Co electron).
    A=sm.A0/3.5;F=1.5
    results['terminal_pair_barrier']={str(s):{'binding_eV':A/(2*s),'remaining_field_gain_eV':.1*F*(6-s),'outgoing_min_absorption_eV':A/(2*s)-.1*F*(6-s)} for s in [5.1,5.5,5.8,5.9,6.]}
    results['last_tag_zero_mismatch_s_nm']=(6+np.sqrt(36-2*A/(.1*F)))/2
    results['residual_binding_at_12nm_in_kBT']={str(T):A/12/(sm.KB*T) for T in [80,300]}
    results['distance_for_binding_below_point1_kBT_nm']={str(T):A/(.1*sm.KB*T) for T in [80,300]}
    results['environment']={'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'mpmath':mp.__version__}
    results['snapshot_sha256']={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((HERE/'snapshot').rglob('*')) if p.is_file()}
    (HERE/'independent_checks.json').write_text(json.dumps(results,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(results['schur_and_lumping'],indent=2),flush=True)

if __name__=='__main__':main()
