#!/usr/bin/env python3
"""Conditional equilibrium harmonic-bath benchmark; no fitted material spectrum.
Energy eV, physical time s; Fourier variable u=t/hbar and contour tau are eV^-1.
Exact quantum Ohmic exponential-cutoff bath plus discrete modes.
No empirical linewidth or posthoc detailed-balance enforcement in rate evaluation.
"""
from pathlib import Path
from functools import lru_cache
import json,csv,sys,hashlib,platform,os
os.environ.setdefault('MPLBACKEND','Agg')
import numpy as np
from scipy.special import loggamma,gammaln,digamma,polygamma,logsumexp
from scipy.integrate import quad
from scipy.optimize import brentq
import scipy
KB=8.617333262145e-5;HBAR=6.582119569e-16
ROOT=Path(__file__).resolve().parent
R=ROOT/'results';R.mkdir(exist_ok=True)

def savecsv(name,rows):
    with (R/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

@lru_cache(None)
def thermal_weights(T,energy,S,L=40):
    ls=np.arange(-L,L+1);x=energy/(KB*T);nb=1/np.expm1(x)
    z=2*S*np.sqrt(nb*(nb+1));j=np.arange(220)
    n=np.abs(ls)[:,None]
    logI=logsumexp((2*j+n)*np.log(z/2)-gammaln(j+1)-gammaln(j+n+1),axis=1)
    return ls,-S*(2*nb+1)+ls*x/2+logI

def gaussian_logp(E,T,lam=.05):
    return -.5*np.log(4*np.pi*lam*KB*T)-(E-lam)**2/(4*lam*KB*T)

@lru_cache(None)
def ohmic_logp(E,T,Ec,lam=.05,tol=2e-10):
    """P(E), bath energy absorption density. Reorg density λ/Ec exp(-ε/Ec).
    Fourier inversion on the steepest real saddle contour, not a linewidth.
    C(z) exact via thermal Gamma product. Both signs computed independently.
    """
    beta=1/(KB*T);alpha=lam/Ec;a=1/(beta*Ec)
    def mean(tau):
        return alpha/beta*(digamma(1+a-tau/beta)-digamma(a+tau/beta))
    lo=-1/Ec+1e-9;hi=beta+1/Ec-1e-9
    tau=brentq(lambda t:mean(t)-E,lo,hi,xtol=1e-11,rtol=1e-14)
    z1=a+tau/beta;z2=1+a-tau/beta
    logc0=alpha*(gammaln(z1)+gammaln(z2)-np.log(beta*Ec)-2*gammaln(1+a))
    var=alpha/beta**2*(polygamma(1,z1)+polygamma(1,z2));scale=1/np.sqrt(var)
    def fun(w):
        x=w*scale;y=x/beta
        z=alpha*(loggamma(z1+1j*y)+loggamma(z2-1j*y)-gammaln(z1)-gammaln(z2))+1j*E*x
        return float(np.exp(z).real)
    integral,error=quad(fun,0,np.inf,epsabs=tol,epsrel=tol,limit=350)
    if integral<=0 or error>max(5e-8,abs(integral)*1e-7):
        raise ArithmeticError((E,T,Ec,lam,integral,error,tau))
    return float(E*tau+logc0+np.log(scale*integral/np.pi))

def quantum_logk(dg,T,Ec=None,H=1e-5,L=30,modes=((.15,1.),),tol=2e-10):
    energies=np.array([0.]);weights=np.array([0.])
    for hw,S in modes:
        ls,lp=thermal_weights(T,hw,S,L)
        energies=(energies[:,None]+ls*hw).ravel();weights=(weights[:,None]+lp).ravel()
    # Group equal discrete energies exactly at 1 peV; does not change distributions.
    ee=np.round(energies,12);unique=np.unique(ee)
    pp=np.array([logsumexp(weights[ee==e]) for e in unique])
    if Ec is None: slow=gaussian_logp(-dg-unique,T)
    else:slow=np.array([ohmic_logp(round(float(-dg-e),12),T,Ec,.05,tol) for e in unique])
    return float(np.log(2*np.pi/HBAR)+2*np.log(H)+logsumexp(pp+slow))

def exact_ohmic_var(T,Ec,lam=.05):
    a=KB*T/Ec
    return lam*Ec*(1+2*a*a*polygamma(1,1+a))

def main():
    rows=[];checks={}
    for T in [80,150,300]:
        for eps in [.001,.002,.005,.01,.02,.05,.15]:
            x=eps/(KB*T)
            rows.append(dict(T_K=T,energy_meV=eps*1000,wavenumber_cm_inverse=eps*8065.54394,x_beta_hw=x,variance_ratio=x/2/np.tanh(x/2),excited_thermal_probability=np.exp(-x),period_fs=2*np.pi*HBAR/eps*1e15))
    savecsv('classical_mode_criteria.csv',rows)
    rows=[]
    for tolerance in [.01,.05,.1]:
        x=brentq(lambda x:x/2/np.tanh(x/2)-1-tolerance,.001,5)
        for T in [80,150,300]:rows.append(dict(variance_tolerance=tolerance,T_K=T,max_energy_meV=x*KB*T*1000,max_wavenumber_cm_inverse=x*KB*T*8065.54394,x_max=x,scope='variance criterion only; rate-tail error can be larger'))
    savecsv('classical_variance_thresholds.csv',rows)
    rows=[]
    for T in [80,150,300]:
      for H in [1e-5,1e-4,1e-3,.015,.02]:
        lam=.05;tc=HBAR/np.sqrt(lam*KB*T);kmax=2*np.pi/HBAR*H*H/np.sqrt(4*np.pi*lam*KB*T)
        for eps in [.001,.005,.02]:
            rows.append(dict(T_K=T,H_meV=H*1000,slow_energy_meV=eps*1000,gaussian_tc_fs=tc*1e15,H_tc_over_hbar=H*tc/HBAR,kmax_s=kmax,kmax_tc=kmax*tc,LZ_typical_exponent=2*np.pi*H*H/(eps*np.sqrt(2*lam*KB*T)),tau_relax_for_1percent_kmax_ns=.01/kmax*1e9,scope='sufficient-scale diagnostics not proof; kmax is kernel upper bound'))
    savecsv('weak_coupling_scales.csv',rows)
    rows=[];db=[];quaderr=[];checksign=[]
    for T in [80.,300.]:
      for Ec in [.001,.005,.02]:
        # Independent signs; no rate-ratio enforcement inside algorithm.
        for E in [0.,.01,.05,.15,.5,1.0]:
            lp=ohmic_logp(E,T,Ec);lm=ohmic_logp(-E,T,Ec)
            db.append(abs(lm-lp+E/(KB*T)))
            lq=ohmic_logp(E,T,Ec,tol=1e-12)
            quaderr.append(abs(lp-lq))
        for dg in [-.9,-.65,-.45,-.3,-.2,-.15,-.1,-.05,0.,.05,.15,.3,.65,.9]:
            lc=quantum_logk(dg,T);lq=quantum_logk(dg,T,Ec)
            rows.append(dict(T_K=T,Ec_meV=Ec*1000,delta_G_eV=dg,log_k_existing_MLJ_s=lc,log_k_quantum_continuum_s=lq,k_existing_MLJ_s=np.exp(lc),k_quantum_continuum_s=np.exp(lq),quantum_over_existing=np.exp(lq-lc),log10_ratio=(lq-lc)/np.log(10),low_bath_variance_ratio=exact_ohmic_var(T,Ec)/(.1*KB*T)))
        # finite-T vibronic balance independently in both directions
        for dg in [.05,.2,.65]:db.append(abs(quantum_logk(dg,T,Ec)-quantum_logk(-dg,T,Ec)+dg/(KB*T)))
    savecsv('fully_quantum_continuum_comparison.csv',rows)
    # Both one-mode and two-mode high-frequency baths have lambda=.15 eV and
    # zero-T second cumulant .0225 eV². Their detailed spectral support differs.
    rows=[]
    for T in [80.,300.]:
      for dg in [-.65,-.45,-.3,-.2,-.15,-.1,-.05,0.,.15,.3]:
        single=quantum_logk(dg,T,.005)
        multi=quantum_logk(dg,T,.005,L=18,modes=((.1,.75),(.2,.375)))
        rows.append(dict(T_K=T,delta_G_eV=dg,Ec_meV=5.,single_mode_k_s=np.exp(single),two_mode_k_s=np.exp(multi),two_over_single=np.exp(multi-single)))
    savecsv('matched_lambda_high_mode_comparison.csv',rows)
    # Integrate slow-bath density and moments independently of characteristic function.
    moments=[]
    for T,Ec in [(80.,.005),(80.,.02),(300.,.005)]:
      vals=[]
      for n in [0,1,2]:
        v,e=quad(lambda E:E**n*np.exp(ohmic_logp(round(float(E),13),T,Ec)), -1.,2.,epsabs=2e-8,epsrel=2e-8,limit=120)
        vals.append(v)
      target=exact_ohmic_var(T,Ec)
      moments.append(dict(T_K=T,Ec_meV=Ec*1000,mass=vals[0],mean_eV=vals[1],variance_eV2=vals[2]-vals[1]**2,exact_variance_eV2=target))
    savecsv('quantum_density_moment_checks.csv',moments)
    checks.update(detailed_balance_max_log_error=max(db),quadrature_tolerance_max_log_error=max(quaderr),density_mass_max_error=max(abs(x['mass']-1) for x in moments),density_mean_max_error_eV=max(abs(x['mean_eV']-.05) for x in moments),density_variance_max_error_eV2=max(abs(x['variance_eV2']-x['exact_variance_eV2']) for x in moments))
    assert checks['detailed_balance_max_log_error']<1e-8
    assert checks['quadrature_tolerance_max_log_error']<1e-7
    assert checks['density_mass_max_error']<2e-7
    assert checks['density_mean_max_error_eV']<2e-7
    # classical continuum limit and sideband truncation separately
    lim=[];trunc=[]
    for T in [80.,300.]:
      for dg in [-.3,-.1,0,.1]:
        lc=quantum_logk(dg,T);lq=quantum_logk(dg,T,.0001)
        lim.append(dict(T_K=T,delta_G_eV=dg,Ec_meV=.1,quantum_over_MLJ=np.exp(lq-lc)))
      for dg in [-.9,.65]:trunc.append(abs(quantum_logk(dg,T,.005,L=30)-quantum_logk(dg,T,.005,L=40)))
    savecsv('classical_limit_check.csv',lim)
    checks['sideband_L30_L40_max_log_error']=max(trunc)
    assert max(trunc)<1e-8
    checks['all_assertions_passed']=True
    checks['python']=platform.python_version();checks['numpy']=np.__version__;checks['scipy']=scipy.__version__
    (R/'validation.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
    print(json.dumps(checks,indent=2),flush=True)

if __name__=='__main__':main()
