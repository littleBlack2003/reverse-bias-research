"""Spatial escape audit. Energies eV, positions nm, rates s^-1. Exploratory only."""
from pathlib import Path
import sys,itertools,json
from functools import lru_cache
import numpy as np
import mpmath as mp
sys.path.insert(0,str(Path(__file__).parent/'inputs'))
import reservoir as old
core=old.core
KB=core.KB; HBAR=core.HBAR; B=float(old.HRES**2); A0=core.A0

def geometry(variant='explicit',s=6.,eta=.5,decay_nm=None):
    if variant in ('colocated','direct'):
        labels=['V','D','C'];x=[-5.,0.,5.];eps=[-.65,0.,.65];b=[1,0,0]
        bondH2=[1e-10,1e-10];bathH2=[B,B]
    elif variant=='explicit':
        labels=['Vo','V','D','C','Co'];x=[-s,-5.,0.,5.,s];eps=[-.65,-.65,0.,.65,.65];b=[1,1,0,0,0]
        bondH2=[eta*B,1e-10,1e-10,eta*B];bathH2=[(1-eta)*B]*2
    elif variant in ('left_only','right_only'):
        if variant=='left_only':
            labels=['Vo','V','D','C'];x=[-s,-5.,0.,5.];eps=[-.65,-.65,0.,.65];b=[1,1,0,0]
            bondH2=[eta*B,1e-10,1e-10];bathH2=[(1-eta)*B,B]
        else:
            labels=['V','D','C','Co'];x=[-5.,0.,5.,s];eps=[-.65,0.,.65,.65];b=[1,0,0,0]
            bondH2=[1e-10,1e-10,eta*B];bathH2=[B,(1-eta)*B]
    else:raise ValueError(variant)
    xb=[-5.,5.] if variant=='colocated' else [-6.,6.]
    if decay_nm is not None:
        # Additional distance suppression at fixed chosen prefactors; NEVER renormalize it away.
        # Internal two 5 nm bonds already have their specified coupling and are untouched.
        for k in range(len(bondH2)):
            if labels[k] in ('Vo','C') and labels[k+1] in ('V','Co'):
                bondH2[k]*=np.exp(-2*(x[k+1]-x[k])/decay_nm)
        bathH2=[h*np.exp(-2*abs(xx-yy)/decay_nm) for h,xx,yy in zip(bathH2,[x[0],x[-1]],xb)]
    return dict(variant=variant,labels=labels,x=x,eps=eps,b=b,bondH2=bondH2,bathH2=bathH2,xb=xb,eta=eta,s=s,decay_nm=decay_nm)

@lru_cache(None)
def graph(n,nval=0,single_pair=False):
    states=[s for s in itertools.product([0,1],repeat=n) if not single_pair or (nval-sum(s[:nval])<=1 and sum(s[nval:])<=1)];idx={s:i for i,s in enumerate(states)};edges=[]
    for i,s in enumerate(states):
        for k in range(n-1):
            if s[k] and not s[k+1]:
                t=list(s);t[k]=0;t[k+1]=1
                if tuple(t) in idx:edges.append((i,idx[tuple(t)],'hop',k))
        for k,side in [(0,0),(n-1,1)]:
            if not s[k]:
                t=list(s);t[k]=1
                if tuple(t) in idx:edges.append((i,idx[tuple(t)],'bath',side))
    return states,edges

def solve(F,T,variant='explicit',s=6.,eta=.5,delta_mu=None,mean_mu=0.,gauge=0.,epsr=3.5,decay_nm=None,kind='quantum',order=24,dps=250,detail=False,single_pair=False):
    mp.mp.dps=dps;m=lambda z:mp.mpf(str(z));beta=1/(m(KB)*m(T));g=geometry(variant,s,eta,decay_nm);n=len(g['x']);states,edges=graph(n,sum(g['b']),single_pair)
    if delta_mu is None:delta_mu=F
    mu=[m(mean_mu)+m(delta_mu)/2-m(gauge),m(mean_mu)-m(delta_mu)/2-m(gauge)]
    band=[m(-.65)-m('.1')*m(F)*m(g['xb'][0])-m(gauge),m(.65)-m('.1')*m(F)*m(g['xb'][1])-m(gauge)]
    E=[]
    for state in states:
        q=[b-ni for b,ni in zip(g['b'],state)]
        en=sum(m(e)*(ni-bi)+qi*(m('.1')*m(F)*m(x)+m(gauge)) for e,ni,bi,qi,x in zip(g['eps'],state,g['b'],q,g['x']))
        if epsr is not None:en+=sum(m(A0)/m(epsr)*q[i]*q[j]/abs(m(g['x'][i])-m(g['x'][j])) for i in range(n) for j in range(i+1,n))
        E.append(en)
    Q=mp.zeros(len(states));rates=[];dbmax=0.;meanmax=0.
    for i,j,ty,k in edges:
        a=E[j]-E[i];diag=None
        if ty=='hop':
            H=np.sqrt(g['bondH2'][k]);dg=float(a)
            # Evaluate downhill/larger direction for numerical fidelity, enforce exact ratio.
            lf=float(core.logkernel(dg,T,kind,H=H));lr=float(core.logkernel(-dg,T,kind,H=H));dbmax=max(dbmax,abs(lf-lr+float(beta*a)))
            if lf>=lr:kf=mp.exp(m(lf));kr=kf*mp.exp(beta*a)
            else:kr=mp.exp(m(lr));kf=kr*mp.exp(-beta*a)
            tdg=a;mean=None;dx=m(g['x'][k+1])-m(g['x'][k]);name=g['labels'][k]+'_'+g['labels'][k+1]
        else:
            tdg=a-mu[k];da=round(float(a-band[k]),15);mur=round(float(mu[k]-band[k]),15)
            diag=old.exchange(T,da,mur,'V' if k==0 else 'C',kind,1.,order)
            factor=m(g['bathH2'][k])/m(B)
            if diag['log_add']>=diag['log_out']:
                kf=mp.exp(m(diag['log_add']))*factor;kr=kf*mp.exp(beta*tdg);mean=band[k]+m(diag['mean_relative_energy_add_eV'])
            else:
                kr=mp.exp(m(diag['log_out']))*factor;kf=kr*mp.exp(-beta*tdg);mean=band[k]+m(diag['mean_relative_energy_out_eV'])
            dbmax=max(dbmax,abs(diag['direct_DB_log_error']));meanmax=max(meanmax,abs(diag['conditional_mean_energy_error_eV']))
            dx=m(g['x'][0 if k==0 else -1])-m(g['xb'][k]);name='bathL' if k==0 else 'bathR'
        rates.append(dict(kf=kf,kr=kr,a=a,tdg=tdg,mean=mean,dx=dx,name=name,diag=diag))
        Q[j,i]+=kf;Q[i,j]+=kr;Q[i,i]-=kf;Q[j,j]-=kr
    p=core.stationary(Q)
    flux={name:mp.mpf(0) for name in [r['name'] for r in rates]};gross={name:mp.mpf(0) for name in flux};qph={name:mp.mpf(0) for name in flux};qel=[mp.mpf(0)]*2
    node=mp.matrix([0]*n);heat=mp.mpf(0);entropy=mp.mpf(0);energy=mp.mpf(0);disp=mp.mpf(0);details=[]
    for (i,j,ty,k),r in zip(edges,rates):
        f=p[i]*r['kf'];rev=p[j]*r['kr'];net=f-rev;name=r['name'];flux[name]+=net;gross[name]+=f+rev
        for z in range(n):node[z]+=net*(states[j][z]-states[i][z])
        heat-=net*r['tdg'];energy+=net*r['a'];entropy+=net*mp.log(f/rev);disp+=net*r['dx']
        if ty=='hop':qph[name]-=net*r['a']
        else:qph[name]+=net*(r['mean']-r['a']);qel[k]+=net*(mu[k]-r['mean'])
        details.append(dict(source=''.join(map(str,states[i])),destination=''.join(map(str,states[j])),process=name,forward_rate_s=float(r['kf']),reverse_rate_s=float(r['kr']),net_flux_s=float(net),gross_flux_s=float(f+rev),energy_change_eV=float(r['a']),phonon_heat_eV_s=float(-net*r['a'] if ty=='hop' else net*(r['mean']-r['a'])),electron_bath_heat_eV_s=float(0 if ty=='hop' else net*(mu[k]-r['mean'])),electron_displacement_nm_s=float(net*r['dx'])))
    R=flux['V_D'];maxgross=max(gross.values());scale=maxgross if delta_mu==0 else max(abs(R),mp.mpf('1e-250'));chem=mu[0]*flux['bathL']+mu[1]*flux['bathR'];powScale=max(abs(heat),abs(chem),maxgross*mp.mpf('1e-90'))
    z=[-beta*(en-(m(mean_mu)-m(gauge))*sum(st)) for en,st in zip(E,states)];zmax=max(z);peq=mp.matrix([mp.exp(a-zmax) for a in z]);peq/=sum(peq)
    row=dict(variant=variant,T_K=T,F_MVcm=F,delta_mu_eV=delta_mu,s_nm=s,eta=eta,decay_nm=decay_nm,epsr=epsr,kernel=kind,number_states=len(states),single_pair_constraint=single_pair,R_s=float(R),bath_left_in_s=float(flux['bathL']),bath_right_out_s=float(-flux['bathR']),total_coupling_H2_eV2=sum(g['bondH2'])+sum(g['bathH2']),branch_budget_nominal_H2_eV2=B,fixed_outer_window_potential_drop_V=1.2*F,represented_bath_potential_drop_V=.1*F*(g['xb'][1]-g['xb'][0]),represented_bath_span_nm=g['xb'][1]-g['xb'][0],max_node_continuity_relative=float(max(abs(z) for z in node)/scale),generator_residual_scaled=float(max(abs(z) for z in Q*p)/maxgross),min_probability=float(min(p)),normalization_error=float(abs(sum(p)-1)),total_bath_heat_eV_s=float(heat),chemical_work_eV_s=float(chem),phonon_heat_eV_s=float(sum(qph.values())),electron_bath_heat_eV_s=float(sum(qel)),left_electron_heat_eV_s=float(qel[0]),right_electron_heat_eV_s=float(qel[1]),heat_chemical_relative=float(abs(heat-chem)/powScale),heat_partition_relative=float(abs(sum(qph.values())+sum(qel)-heat)/powScale),state_energy_derivative_eV_s=float(energy),entropy_kB_s=float(entropy),entropy_heat_relative=float(abs(entropy-beta*heat)/max(abs(entropy),maxgross*mp.mpf('1e-85'))),electron_displacement_nm_s=float(disp),displacement_relative=float(abs(disp-(m(g['xb'][1])-m(g['xb'][0]))*R)/(scale*(m(g['xb'][1])-m(g['xb'][0])))),direct_DB_max_log_error=dbmax,direct_mean_energy_error_eV=meanmax,common_mu_probability_relative=float(max(abs(p[i]/peq[i]-1) for i in range(len(p)))) if delta_mu==0 else None,common_mu_edge_flux_over_gross=float(max(abs(mp.mpf(str(d['net_flux_s']))) for d in details)/maxgross) if delta_mu==0 else None,precision_dps=dps)
    for k,label in enumerate(g['labels']):row['occupation_'+label]=float(sum(p[i]*st[k] for i,st in enumerate(states)))
    vi=g['labels'].index('V');ci=g['labels'].index('C');di=g['labels'].index('D')
    row['inner_pair_probability']=float(sum(p[i] for i,st in enumerate(states) if st[vi]==0 and st[di]==0 and st[ci]==1))
    row['core_formation_gross_s']=float(gross['V_D']);row['core_generation_gross_s']=float(gross['D_C']);row['net_over_generation_gross']=float(R/gross['D_C'])
    row['mean_total_charge_e']=float(sum(p[i]*sum(b-ni for b,ni in zip(g['b'],st)) for i,st in enumerate(states)))
    row['p_most_likely']=float(max(p));row['most_likely_state']=''.join(map(str,states[max(range(len(p)),key=lambda i:p[i])]))
    if detail:return row,p,Q,rates,details,E,mu,g,states,edges
    return row
