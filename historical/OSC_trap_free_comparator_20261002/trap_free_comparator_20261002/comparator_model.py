"""Trap-free native frontier-state pair versus one midpoint defect.
Exploratory nonadiabatic thermal-vibronic Markov graph; eV, nm, s.
Not coherent Zener/WKB and not material fit. All charge states retained.
"""
from pathlib import Path
from functools import lru_cache
import sys,itertools,math
import numpy as np
import mpmath as mp
sys.path.insert(0,str(Path(__file__).resolve().parent/'inputs/previous'))
import reservoir as old
core=old.core
KB=core.KB;HBAR=core.HBAR;QELEC=core.QE;COULOMB=core.A0
INTERNAL_BUDGET=2e-10

@lru_cache(None)
def graph(topology):
 n=2 if topology=='direct' else 3
 states=list(itertools.product([0,1],repeat=n));ind={s:i for i,s in enumerate(states)};edges=[]
 for i,s in enumerate(states):
  for k in range(n-1):
   if s[k] and not s[k+1]:
    t=list(s);t[k]=0;t[k+1]=1;edges.append((i,ind[tuple(t)],'hop',k))
  for k,side in [(0,0),(n-1,1)]:
   if not s[k]:
    t=list(s);t[k]=1;edges.append((i,ind[tuple(t)],'bath',side))
 return states,edges

def solve(F,T,topology='direct',span=10.,mode='fixed_budget',decay=None,delta_mu=None,mean_mu=0.,gauge=0.,epsr=3.5,bath_scale=(1.,1.),order=24,dps=220,detail=False):
 mp.mp.dps=dps;m=lambda v:mp.mpf(str(v));beta=1/(m(KB)*m(T))
 n=2 if topology=='direct' else 3
 x=[-span/2,span/2] if n==2 else [-span/2,0,span/2]
 bare=[-.65,.65] if n==2 else [-.65,0,.65];cores=[1]+[0]*(n-1)
 labels=['V','C'] if n==2 else ['V','D','C']
 if mode=='fixed_budget':h=[np.sqrt(INTERNAL_BUDGET/(n-1))]*(n-1)
 elif mode=='distance':
  assert decay is not None and decay>0
  # Same distance law calibrated to H(5 nm)=10 micro-eV. Never renormalized.
  h=[1e-5*np.exp(-((x[k+1]-x[k])-5)/decay) for k in range(n-1)]
 else:raise ValueError(mode)
 states,edges=graph(topology)
 dmu=.1*F*span if delta_mu is None else delta_mu
 mu=[m(mean_mu)+m(dmu)/2-m(gauge),m(mean_mu)-m(dmu)/2-m(gauge)]
 bands=[m(-.65)-m('.1')*m(F)*m(x[0])-m(gauge),m(.65)-m('.1')*m(F)*m(x[-1])-m(gauge)]
 E=[]
 for s in states:
  q=[c-ni for c,ni in zip(cores,s)]
  en=sum(m(e)*(ni-c)+m(qi)*(m('.1')*m(F)*m(xx)+m(gauge)) for e,ni,c,qi,xx in zip(bare,s,cores,q,x))
  if epsr is not None:en+=sum(m(COULOMB)/m(epsr)*m(q[i]*q[j])/abs(m(x[i])-m(x[j])) for i in range(n) for j in range(i+1,n))
  E.append(en)
 Q=mp.zeros(len(states));rates=[];db=0.;mer=0.
 for i,j,kind,k in edges:
  a=E[j]-E[i];diag=None;mean=None
  if kind=='hop':
   lf=float(core.logkernel(float(a),T,H=h[k]));lr=float(core.logkernel(float(-a),T,H=h[k]))
   db=max(db,abs(lf-lr+float(beta*a)))
   if lf>=lr:kf=mp.exp(m(lf));kr=kf*mp.exp(beta*a)
   else:kr=mp.exp(m(lr));kf=kr*mp.exp(-beta*a)
   tdg=a;name=labels[k]+'_'+labels[k+1];dx=m(x[k+1])-m(x[k])
  else:
   tdg=a-mu[k];dx=m(0);name=['bathL','bathR'][k]
   da=round(float(a-bands[k]),15);mur=round(float(mu[k]-bands[k]),15)
   diag=old.exchange(T,da,mur,['V','C'][k],order=order)
   db=max(db,abs(diag['direct_DB_log_error']));mer=max(mer,abs(diag['conditional_mean_energy_error_eV']))
   scale=m(bath_scale[k])
   if diag['log_add']>=diag['log_out']:
    kf=mp.exp(m(diag['log_add']))*scale;kr=kf*mp.exp(beta*tdg);mean=bands[k]+m(diag['mean_relative_energy_add_eV'])
   else:
    kr=mp.exp(m(diag['log_out']))*scale;kf=kr*mp.exp(-beta*tdg);mean=bands[k]+m(diag['mean_relative_energy_out_eV'])
  Q[j,i]+=kf;Q[i,j]+=kr;Q[i,i]-=kf;Q[j,j]-=kr
  rates.append(dict(kf=kf,kr=kr,a=a,tdg=tdg,mean=mean,kind=kind,side=k,name=name,dx=dx))
 p=core.stationary(Q)
 flux={r['name']:m(0) for r in rates};gross={name:m(0) for name in flux};node=mp.matrix([0]*n)
 heat=m(0);ph=m(0);el=m(0);entropy=m(0);stateen=m(0);disp=m(0);details=[]
 for (i,j,kind,k),r in zip(edges,rates):
  f=p[i]*r['kf'];rev=p[j]*r['kr'];net=f-rev;name=r['name'];flux[name]+=net;gross[name]+=f+rev
  for l in range(n):node[l]+=net*(states[j][l]-states[i][l])
  heat-=net*r['tdg'];stateen+=net*r['a'];disp+=net*r['dx']
  if f>0 and rev>0:entropy+=net*mp.log(f/rev)
  if kind=='hop':pq=-net*r['a'];eq=m(0)
  else:pq=net*(r['mean']-r['a']);eq=net*(mu[k]-r['mean'])
  ph+=pq;el+=eq
  details.append(dict(source=''.join(map(str,states[i])),destination=''.join(map(str,states[j])),process=name,delta_energy_eV=float(r['a']),thermodynamic_delta_eV=float(r['tdg']),source_probability=float(p[i]),destination_probability=float(p[j]),kf_s=float(r['kf']),kr_s=float(r['kr']),net_flux_s=float(net),gross_flux_s=float(f+rev),phonon_heat_eV_s=float(pq),electronic_heat_eV_s=float(eq)))
 R=flux['V_C' if n==2 else 'D_C'];g=max(gross.values());sc=g if dmu==0 or 0 in bath_scale else max(abs(R),m('1e-250'))
 power=mu[0]*flux['bathL']+mu[1]*flux['bathR'];scale=max(abs(power),abs(heat),g*m('1e-150'))
 probs=[-beta*(ee-m(mean_mu-gauge)*sum(s)) for ee,s in zip(E,states)];z=max(probs);peq=mp.matrix([mp.exp(v-z) for v in probs]);peq/=sum(peq)
 pair=[0]*n;pair[-1]=1;paired=states.index(tuple(pair));neutral=states.index(tuple(cores))
 row=dict(topology=topology,T_K=T,F_MVcm=F,span_nm=span,voltage_drop_V=.1*F*span,delta_mu_eV=dmu,mode=mode,decay_nm=decay,epsr=epsr,internal_total_H2_eV2=sum(z*z for z in h),internal_H_eV=h[0],bath_H2_eV2_each=old.HRES**2,bath_left_scale=bath_scale[0],bath_right_scale=bath_scale[1],number_states=len(states),R_pair_s=float(R),conditional_J_for_Npair_1e15_and_d100nm_Acm2=float(R)*QELEC*1e10,required_Npair_for_50mA_at_d100nm_cm3=.05/(float(R)*QELEC*1e-5) if R>0 and dmu>0 and all(bath_scale) else None,pair_probability=float(p[paired]),neutral_probability=float(p[neutral]),pair_minus_neutral_eV=float(E[paired]-E[neutral]),min_probability=float(min(p)),normalization_error=float(abs(sum(p)-1)),generator_residual=float(max(abs(z) for z in Q*p)/g),continuity_relative=float(max(abs(z) for z in node)/sc),heat_eV_s=float(heat),phonon_heat_eV_s=float(ph),electron_heat_eV_s=float(el),chemical_power_eV_s=float(power),heat_vs_power_relative=float(abs(heat-power)/scale),heat_partition_relative=float(abs(ph+el-heat)/scale),state_energy_derivative_eV_s=float(stateen),entropy_kB_s=float(entropy),entropy_balance_relative=float(abs(entropy-beta*heat)/max(abs(entropy),g*m('1e-150'))),electron_displacement_nm_s=float(disp),displacement_relative=float(abs(disp-m(span)*R)/(m(span)*sc)),direct_DB_max_log_error=db,conditional_spectrum_mean_error_eV=mer,equilibrium_probability_relative=float(max(abs(p[i]/peq[i]-1) for i in range(len(p)))) if dmu==0 else None,equilibrium_edge_net_over_gross=max(abs(r['net_flux_s'])/(r['gross_flux_s'] or 1) for r in details) if dmu==0 else None,precision_dps=dps,quadrature_order=order)
 if detail:return row,p,Q,rates,details,E,states,edges
 return row

def independent_tree_probabilities(Q):
 """Enumerate all directed in-trees for four states, independent of linear solve."""
 n=Q.rows;assert n==4
 weights=[]
 for root in range(n):
  others=[i for i in range(n) if i!=root];wt=mp.mpf(0)
  options=[[j for j in range(n) if j!=i and Q[j,i]>0] for i in others]
  for dest in itertools.product(*options):
   arrows=dict(zip(others,dest));ok=True
   for start in others:
    visited=set();i=start
    while i!=root:
     if i in visited:ok=False;break
     visited.add(i);i=arrows[i]
    if not ok:break
   if ok:wt+=mp.fprod(Q[j,i] for i,j in arrows.items())
  weights.append(wt)
 return mp.matrix([v/sum(weights) for v in weights])
