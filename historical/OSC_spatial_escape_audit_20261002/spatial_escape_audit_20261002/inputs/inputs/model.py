"""Reversible 3-site / 8-state benchmark reconstructed 2026-10-02.
Electron energy eV; distance nm; time s; F=MV/cm (positive favors electrons +x).
Single-occupancy tagged orbitals, grand-canonical baths behind finite exchange.
"""
from pathlib import Path
from functools import lru_cache
import json, itertools
import numpy as np
from scipy.special import logsumexp,gammaln,eval_genlaguerre
import mpmath as mm
P=json.loads((Path(__file__).parent/'parameters.json').read_text())
KB=8.617333262145e-5; HBAR=6.582119569e-16; QE=1.602176634e-19; A0=1.4399645478
STATES=list(itertools.product([0,1],repeat=3)); INDEX={s:i for i,s in enumerate(STATES)}
# Every undirected edge is oriented by a physical process, not index order.
EDGES=[]
for i,s in enumerate(STATES):
    for a,b,name in [(0,1,'VD'),(1,2,'DC')]:
        if s[a] and not s[b]:
            t=list(s);t[a]=0;t[b]=1
            EDGES.append((i,INDEX[tuple(t)],name))
    for a,name in [(0,'V'),(2,'C')]:
        if not s[a]:
            t=list(s);t[a]=1
            EDGES.append((i,INDEX[tuple(t)],name))
assert len(EDGES)==12
@lru_cache(None)
def sidebands(T,L=60):
    l=np.arange(-L,L+1);x=P['quantum_mode_eV']/(KB*T);nb=1/np.expm1(x);S=P['Huang_Rhys_S']
    z=2*S*np.sqrt(nb*(nb+1));j=np.arange(180);nu=np.abs(l)[:,None]
    logI=logsumexp((2*j+nu)*np.log(z/2)-gammaln(j+1)-gammaln(j+nu+1),axis=1)
    lp=-S*(2*nb+1)+l*x/2+logI
    # Intentionally NOT normalized after truncation.
    return l,lp

def logkernel(dg,T,kind='quantum',H=1e-5,L=60):
    dg=np.asarray(dg)
    if kind=='classical':
        lam=P['classical_lambda_eV']
        return np.log(2*np.pi/HBAR)+2*np.log(H)-.5*np.log(4*np.pi*lam*KB*T)-(dg+lam)**2/(4*lam*KB*T)
    lam=P['quantum_low_frequency_lambda_eV'];hw=P['quantum_mode_eV']
    if kind=='naive_ground_mlj':
        l=np.arange(L+1);S=P['Huang_Rhys_S'];lp=-S+l*np.log(S)-gammaln(l+1)
    else:l,lp=sidebands(T,L)
    return np.log(2*np.pi/HBAR)+2*np.log(H)-.5*np.log(4*np.pi*lam*KB*T)+logsumexp(lp-(dg[...,None]+lam+l*hw)**2/(4*lam*KB*T),axis=-1)

def franck_condon_double_lograte(dg,T,N=64):
    """Independent thermal-initial-state m,n sum, no Bessel formula."""
    m,n=np.meshgrid(np.arange(N),np.arange(N),indexing='ij');mi=np.minimum(m,n);ma=np.maximum(m,n)
    S=P['Huang_Rhys_S']; hw=P['quantum_mode_eV'];lam=P['quantum_low_frequency_lambda_eV']
    lag=eval_genlaguerre(mi,ma-mi,S)
    with np.errstate(divide='ignore'):
        lf=-S+gammaln(mi+1)-gammaln(ma+1)+(ma-mi)*np.log(S)+2*np.log(np.abs(lag))
    lt=np.log(-np.expm1(-hw/(KB*T)))-m*hw/(KB*T)
    return np.log(2*np.pi/HBAR)+2*np.log(1e-5)-.5*np.log(4*np.pi*lam*KB*T)+logsumexp(lt+lf-(dg+lam+(n-m)*hw)**2/(4*lam*KB*T))

def energy(F,mp=False):
    conv=mm.mpf if mp else float
    cv=lambda x:conv(str(x));x=list(map(cv,P['positions_nm']));eps=list(map(cv,P['bare_electron_energies_eV']));b=P['core_charge_e']
    phi=[cv('.1')*cv(F)*xx for xx in x];a=cv(A0)/cv(P['relative_dielectric'])
    E=[]
    for s in STATES:
        q=[b[i]-s[i] for i in range(3)]
        ee=sum(eps[i]*(s[i]-b[i])+q[i]*phi[i] for i in range(3))
        ee+=sum(a*q[i]*q[j]/abs(x[i]-x[j]) for i in range(3) for j in range(i+1,3))
        E.append(ee)
    return E

def build(F,T,kind='quantum',gamma=1e6,equilibrium=False,dps=110,naive_irreversible_contacts=False):
    mm.mp.dps=dps;m=mm.mpf;beta=1/(m(str(KB))*m(str(T)));E=energy(F,True)
    mu={'V':m('.5')*m(str(F)),'C':-m('.5')*m(str(F))}
    if equilibrium:mu={'V':m('.023'),'C':m('.023')}
    rates=[];Q=mm.zeros(8)
    gammas=gamma if isinstance(gamma,dict) else {'V':gamma,'C':gamma}
    for i,j,name in EDGES:
        dg=E[j]-E[i]
        if name in ('VD','DC'):
            kf=mm.exp(m(str(float(logkernel(float(dg),T,kind,H=P['couplings_eV'][name])))))
            if kind=='naive_ground_mlj':kr=mm.exp(m(str(float(logkernel(float(-dg),T,kind,H=P['couplings_eV'][name])))))
            else:
                # The FULL thermal kernel is independently checked in both directions.
                # Exact ratio assembly prevents float roundoff leaking into rare currents.
                kr=kf*mm.exp(beta*dg)
            thermodg=dg;dn=0
        else:
            thermodg=dg-mu[name];dn=1;g=m(str(gammas[name]))
            kf=g/(1+mm.exp(beta*thermodg));kr=g/(1+mm.exp(-beta*thermodg))
            if naive_irreversible_contacts:
                if name=='V':kf=g;kr=m(0)
                else:kf=m(0);kr=g
        rates.append((kf,kr,dg,thermodg,dn))
        Q[j,i]+=kf;Q[i,j]+=kr;Q[i,i]-=kf;Q[j,j]-=kr
    return E,mu,rates,Q

def stationary(Q):
    n=Q.rows;A=Q.copy();scale=max(abs(z) for z in Q)
    A=A/scale
    for j in range(n):A[n-1,j]=1
    b=mm.matrix([0]*(n-1)+[1]);p=mm.lu_solve(A,b)
    assert min(p)>0, 'Never clip negative stationary probabilities'
    assert abs(sum(p)-1)<mm.mpf('1e-80')
    return p

def solve(F,T,kind='quantum',gamma=1e6,equilibrium=False,dps=110,detail=False,naive_irreversible_contacts=False):
    E,mu,rates,Q=build(F,T,kind,gamma,equilibrium,dps,naive_irreversible_contacts)
    p=stationary(Q);flux={name:mm.mpf(0) for name in ['VD','DC','V','C']};gross={k:mm.mpf(0) for k in flux}
    en=mm.mpf(0);heat=mm.mpf(0);ep=mm.mpf(0);edge=[]
    for (i,j,name),(kf,kr,dg,tdg,dn) in zip(EDGES,rates):
        f=p[i]*kf;r=p[j]*kr;net=f-r;flux[name]+=net;gross[name]+=f+r
        en+=net*dg;heat-=net*tdg
        if f>0 and r>0:ep+=net*mm.log(f/r)
        edge.append((i,j,name,f,r,net,dg,tdg))
    chem=mu['V']*flux['V']+mu['C']*flux['C'];rate=(flux['VD']+flux['DC'])/2
    maxgross=max(gross.values());scale=max(abs(rate),mm.mpf('1e-200')) if not (equilibrium or F==0) else maxgross
    continuity=max(abs(flux['VD']-flux['DC']),abs(flux['V']-flux['VD']),abs(flux['C']+flux['DC']))/scale
    resid=max(abs(z) for z in Q*p)/maxgross
    powerrel=abs(heat-chem)/max(abs(heat),abs(chem),maxgross*mm.mpf('1e-95'))
    row={'T_K':T,'F_MVcm':F,'kernel':kind,'GammaV_s':float(gamma['V'] if isinstance(gamma,dict) else gamma),'GammaC_s':float(gamma['C'] if isinstance(gamma,dict) else gamma),
      'R_s':float(rate),'V_bath_in_s':float(flux['V']),'C_bath_out_s':float(-flux['C']),
      'VD_gross_s':float(gross['VD']),'DC_gross_s':float(gross['DC']),
      'nV':float(sum(p[i]*s[0] for i,s in enumerate(STATES))),'nD':float(sum(p[i]*s[1] for i,s in enumerate(STATES))),'nC':float(sum(p[i]*s[2] for i,s in enumerate(STATES))),
      'total_charge_e':float(sum(p[i]*(1-sum(s)) for i,s in enumerate(STATES))),
      'continuity_relative_error':float(continuity),'generator_scaled_residual':float(resid),'energy_balance_eV_s':float(en),'heat_bath_eV_s':float(heat),'chemical_power_eV_s':float(chem),
      'heat_vs_chemical_relative_error':float(powerrel),'entropy_production_kB_s':float(ep),
      'minimum_probability':float(min(p)),'probability_normalization_error':float(abs(sum(p)-1)),
      'conditional_J_Acm2':float(rate)*QE*P['independent_graph_density_cm3_for_conditional_conversion']*P['device_thickness_cm_for_conditional_conversion'],
      'local_internal_displacement_nm_s':float(rate)*10,'conditional_internal_Ramo_fraction':.1,
      'local_power_W_per_graph':float(heat)*QE,'working_dps':dps}
    if detail:return row,p,rates,Q,edge,E,mu
    return row

def numpy_stationary(Q):
    a=np.array(Q.tolist(),float);scale=np.max(np.abs(a));a=a/scale;a[-1,:]=1;b=np.zeros(8);b[-1]=1
    return np.linalg.solve(a,b),np.linalg.cond(a)

def spanning_trees():
    result=[]
    for chosen in itertools.combinations(range(len(EDGES)),7):
        adjacency=[[] for _ in range(8)]
        for e in chosen:
            i,j,_=EDGES[e];adjacency[i].append((j,e));adjacency[j].append((i,e))
        reached={0};queue=[0]
        for v in queue:
            for w,e in adjacency[v]:
                if w not in reached:reached.add(w);queue.append(w)
        if len(reached)!=8:continue
        oriented=[]
        for root in range(8):
            done={root};queue=[root];arcs=[]
            for target in queue:
                for source,e in adjacency[target]:
                    if source in done:continue
                    done.add(source);queue.append(source)
                    arcs.append((e,0 if EDGES[e][0]==source else 1))
            oriented.append(arcs)
        result.append(oriented)
    return result

def tree_stationary(rates,trees):
    logs=np.array([[float(mm.log(r[0])),float(mm.log(r[1]))] for r in rates])
    weights=[logsumexp([sum(logs[e,d] for e,d in t[root]) for t in trees]) for root in range(8)]
    return np.exp(weights-logsumexp(weights))
