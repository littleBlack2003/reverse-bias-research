"""Finite normalized-DOS reversible reservoir exchange, new phase 2026-10-02.
Equations were recorded first in THEORY_BEFORE_CODE.md. All parameters exploratory.
"""
from pathlib import Path
import sys
from functools import lru_cache
import numpy as np
from scipy.special import logsumexp
from scipy.special import roots_legendre
import mpmath as mm
sys.path.insert(0,str(Path(__file__).parent/'inputs'))
import model as core
KB=core.KB;HBAR=core.HBAR;QE=core.QE
GAMMA0=1e6;W0=1.;HRES=np.sqrt(GAMMA0*HBAR*W0/(2*np.pi))
@lru_cache(None)
def quadrature_nodes(T,da,murel,alpha,kind,W,order,support):
 lo,hi=((-W,0) if alpha=='V' else (0,W)) if support=='one_sided' else (-W/2,W/2)
 lam=.2 if kind=='classical' else .05;hw=.15
 knots=[lo,hi]
 # Split at physical Gaussian centers and the Fermi edge plus thermal shoulders.
 for z in [murel,murel-12*KB*T,murel+12*KB*T,da]:
  if lo<z<hi:knots.append(z)
 for direction in [-1,1]:
  for l in ([0] if kind=='classical' else range(-40,41)):
   z=da+direction*(lam+l*hw)
   if lo<z<hi:knots.append(z)
 knots=np.unique(knots);x,w=roots_legendre(order);nodes=[];weights=[]
 for a,b in zip(knots[:-1],knots[1:]):
  nodes.extend((a+b)/2+(b-a)/2*x);weights.extend(w*(b-a)/2/W)
 return np.array(nodes),np.array(weights)
@lru_cache(None)
def exchange(T,da,murel,alpha,kind='quantum',W=1.,order=24,wrong_fermi=False,support='one_sided'):
 """All energies relative to this bath's band edge; da is charge-conditioned.
 Integral spectral weight=1. HRES is fixed even when W changes.
 Outputs independent direct add/out logs and conditional energy first moments.
 """
 u,w=quadrature_nodes(T,da,murel,alpha,kind,W,order,support)
 lf=-np.logaddexp(0,(u-murel)/(KB*T));lh=-np.logaddexp(0,-(u-murel)/(KB*T))
 ladd=core.logkernel(da-u,T,kind,H=HRES)+lf+np.log(w)
 lout=core.logkernel(u-da,T,kind,H=HRES)+(lf if wrong_fermi else lh)+np.log(w)
 a=float(logsumexp(ladd));b=float(logsumexp(lout));ma=float(np.sum(np.exp(ladd-a)*u));mb=float(np.sum(np.exp(lout-b)*u))
 return dict(log_add=a,log_out=b,mean_relative_energy_add_eV=ma,mean_relative_energy_out_eV=mb,direct_DB_log_error=a-b+(da-murel)/(KB*T),conditional_mean_energy_error_eV=ma-mb,node_count=len(u),DOS_mass=float(sum(w)),W_eV=W,H_res_eV=HRES,total_coupling_H2_eV2=HRES**2)

def solve(F,T,kind='quantum',delta_mu=None,mean_mu=0.,bath_model='finite_dos',bath_kind=None,W=1.,order=24,dps=130,detail=False,wrong_fermi=False,support='one_sided'):
 mm.mp.dps=dps;m=mm.mpf;beta=1/(m(str(KB))*m(str(T)));E=core.energy(F,True)
 if delta_mu is None:delta_mu=F
 mu={'V':m(str(mean_mu))+m(str(delta_mu))/2,'C':m(str(mean_mu))-m(str(delta_mu))/2}
 band={'V':m('-.65')+m(str(F))/2,'C':m('.65')-m(str(F))/2}
 rates=[];Q=mm.zeros(8);bath_kind=bath_kind or kind;details=[]
 for i,j,name in core.EDGES:
  a=E[j]-E[i]
  if name in ['VD','DC']:
   kf=mm.exp(m(str(float(core.logkernel(float(a),T,kind,H=core.P['couplings_eV'][name])))));kr=kf*mm.exp(beta*a);tdg=a;mean=None;diag=None
  else:
   tdg=a-mu[name]
   if bath_model=='legacy':
    kf=m(str(GAMMA0))/(1+mm.exp(beta*tdg));kr=m(str(GAMMA0))/(1+mm.exp(-beta*tdg));mean=a;diag=None
   else:
    da=round(float(a-band[name]),15);mur=round(float(mu[name]-band[name]),15)
    diag=exchange(T,da,mur,name,bath_kind,W,order,wrong_fermi,support)
    if wrong_fermi:
     kf=mm.exp(m(str(diag['log_add'])));kr=mm.exp(m(str(diag['log_out'])))
    elif diag['log_add']>=diag['log_out']:
     kf=mm.exp(m(str(diag['log_add'])));kr=kf*mm.exp(beta*tdg)
    else:
     kr=mm.exp(m(str(diag['log_out'])));kf=kr*mm.exp(-beta*tdg)
    # Same conditional spectrum follows analytically from pointwise detailed balance.
    me=diag['mean_relative_energy_add_eV'] if diag['log_add']>=diag['log_out'] else diag['mean_relative_energy_out_eV']
    mean=band[name]+m(str(me))
  rates.append((kf,kr,a,tdg,mean,diag));Q[j,i]+=kf;Q[i,j]+=kr;Q[i,i]-=kf;Q[j,j]-=kr
 p=core.stationary(Q);flux={z:m(0) for z in ['VD','DC','V','C']};gross={z:m(0) for z in flux};heat=m(0);stateen=m(0);entropy=m(0);phononheat=m(0);fermiheat=m(0);maxdb=0;maxmeanerr=0
 for (i,j,name),(kf,kr,a,tdg,mean,diag) in zip(core.EDGES,rates):
  forward=p[i]*kf;reverse=p[j]*kr;net=forward-reverse;flux[name]+=net;gross[name]+=forward+reverse;heat-=net*tdg;stateen+=net*a;entropy+=net*mm.log(forward/reverse)
  if name in ['V','C']:
   phononheat+=net*(mean-a);fermiheat+=net*(mu[name]-mean)
   if diag:maxdb=max(maxdb,abs(diag['direct_DB_log_error']));maxmeanerr=max(maxmeanerr,abs(diag['conditional_mean_energy_error_eV']))
  else:phononheat-=net*a
  details.append(dict(source=''.join(map(str,core.STATES[i])),destination=''.join(map(str,core.STATES[j])),process=name,source_probability=float(p[i]),destination_probability=float(p[j]),forward_rate_s=float(kf),reverse_rate_s=float(kr),effective_Gamma_s=float(kf+kr),net_flux_s=float(net),gross_flux_s=float(forward+reverse),addition_energy_eV=float(a),thermodynamic_delta_eV=float(tdg),mean_exchange_energy_eV=float(mean) if mean is not None else '',direct_integral_DB_error=diag['direct_DB_log_error'] if diag else '',direct_mean_energy_error_eV=diag['conditional_mean_energy_error_eV'] if diag else ''))
 rate=flux['DC'];chem=mu['V']*flux['V']+mu['C']*flux['C'];maxgross=max(gross.values());eq=(delta_mu==0);sc=maxgross if eq else max(abs(rate),m('1e-200'))
 cont=max(abs(flux['VD']-flux['DC']),abs(flux['V']-flux['VD']),abs(flux['C']+flux['DC']))/sc
 heaterr=abs(heat-chem)/max(abs(heat),abs(chem),maxgross*m('1e-90'))
 energyerr=abs(phononheat+fermiheat-heat)/max(abs(heat),maxgross*m('1e-90'))
 enterr=abs(entropy-heat*beta)/max(abs(entropy),maxgross*m('1e-85'))
 row=dict(T_K=T,F_MVcm=F,delta_mu_eV=delta_mu,mean_mu_eV=mean_mu,internal_kernel=kind,bath_kernel=bath_kind,bath_model=bath_model,DOS_support=support,W_eV=W,R_s=float(rate),nV=float(sum(p[i]*s[0] for i,s in enumerate(core.STATES))),nD=float(sum(p[i]*s[1] for i,s in enumerate(core.STATES))),nC=float(sum(p[i]*s[2] for i,s in enumerate(core.STATES))),total_charge_e=float(sum(p[i]*(1-sum(s)) for i,s in enumerate(core.STATES))),min_probability=float(min(p)),continuity_relative_error=float(cont),generator_scaled_residual=float(max(abs(z) for z in Q*p)/maxgross),chemical_power_eV_s=float(chem),total_bath_heat_eV_s=float(heat),phonon_bath_heat_eV_s=float(phononheat),electron_bath_heat_eV_s=float(fermiheat),state_energy_derivative_eV_s=float(stateen),heat_vs_chemical_relative_error=float(heaterr),heat_partition_relative_error=float(energyerr),entropy_production_kB_s=float(entropy),entropy_vs_heat_relative_error=float(enterr),direct_integral_DB_max_log_error=maxdb,direct_integral_mean_energy_max_error_eV=maxmeanerr,conditional_J_Acm2=float(rate)*QE*1e10,precision_dps=dps)
 if detail:return row,p,rates,Q,details,E,mu
 return row
