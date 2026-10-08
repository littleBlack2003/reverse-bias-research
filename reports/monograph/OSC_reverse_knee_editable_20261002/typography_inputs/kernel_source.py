#!/usr/bin/env python3
"""New 2026-10-02 mechanism benchmarks, NOT sample fits or DD predictions.
Units: eV, nm, K, s, cm, V/cm. All model assumptions in parameters.json.
"""
from pathlib import Path
import os
os.environ.setdefault('MPLBACKEND','Agg')
os.environ.setdefault('MPLCONFIGDIR','/tmp/osc_temperature_matplotlib')
import json, csv, sys, hashlib, platform
import numpy as np
import scipy
from scipy.special import logsumexp, gammaln, expit, eval_genlaguerre
from scipy.integrate import quad
from scipy.optimize import brentq
import sympy as sp
import matplotlib
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
R=ROOT/'results'; R.mkdir(exist_ok=True)
FIG=ROOT/'figures'; FIG.mkdir(exist_ok=True)
KB=8.617333262145e-5
HBAR=6.582119569e-16
Q=1.602176634e-19
EPS0_CM=8.8541878128e-14
COULOMB=1.4399645478 # eV nm, vacuum
P=json.loads((ROOT/'parameters.json').read_text())

def writecsv(name, rows):
    with (R/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def classical_lograte(dg,T,H=P['elementary_H_eV'],lam=P['classical_lambda_eV']):
    return np.log(2*np.pi/HBAR)+2*np.log(H)-.5*np.log(4*np.pi*lam*KB*T)-(np.asarray(dg)+lam)**2/(4*lam*KB*T)

def ma_lograte(dg,T,nu=1e12,r=5.,a=.5):
    return np.log(nu)-2*r/a-np.maximum(np.asarray(dg),0)/(KB*T)

def fc_logweights(T,S=1.,hw=.15,L=100):
    """Thermal displaced-harmonic-oscillator distribution of bath quanta.
    I_-l=I_l; negative l are retained even when very small.
    Do not renormalize truncation: check omitted mass instead.
    """
    l=np.arange(-L,L+1); x=hw/(KB*T); nb=1/np.expm1(x)
    z=2*S*np.sqrt(nb*(nb+1))
    # Bessel I via convergent positive series, avoids underflow at large |l|.
    j=np.arange(0,180)
    order=np.abs(l)[:,None]
    terms=(2*j+order)*np.log(z/2)-gammaln(j+1)-gammaln(j+order+1)
    logI=logsumexp(terms,axis=1)
    lp=-S*(2*nb+1)+l*x/2+logI
    return l,lp

def quantum_lograte(dg,T,H=P['elementary_H_eV'],lams=P['quantum_classical_bath_lambda_eV'],S=P['quantum_Huang_Rhys_S'],hw=P['quantum_mode_hw_eV'],L=100):
    l,lp=fc_logweights(T,S,hw,L)
    dg=np.asarray(dg)
    return (np.log(2*np.pi/HBAR)+2*np.log(H)-.5*np.log(4*np.pi*lams*KB*T)
      +logsumexp(lp-(dg[...,None]+lams+l*hw)**2/(4*lams*KB*T),axis=-1))

def ground_mlj_lograte(dg,T,H=P['elementary_H_eV'],lams=P['quantum_classical_bath_lambda_eV'],S=P['quantum_Huang_Rhys_S'],hw=P['quantum_mode_hw_eV'],L=100):
    l=np.arange(L+1); lp=-S+l*np.log(S)-gammaln(l+1)
    dg=np.asarray(dg)
    return (np.log(2*np.pi/HBAR)+2*np.log(H)-.5*np.log(4*np.pi*lams*KB*T)
      +logsumexp(lp-(dg[...,None]+lams+l*hw)**2/(4*lams*KB*T),axis=-1))

def cycle(a,b,c,d):
    # 0->1 valence supply a; its reverse b. 1->0 conduction release c; reverse d.
    s=a+b+c+d; occ=(a+d)/s
    return (a*c-b*d)/s,occ

def log_depleted_cycle(la,lc):
    return la+lc-logsumexp(np.stack(np.broadcast_arrays(la,lc)),axis=0)

def gaussian_filling(mu,T,sigma=.0668):
    f=lambda z: np.exp(-z*z/2)/np.sqrt(2*np.pi)*expit((mu-sigma*z)/(KB*T))
    return quad(f,-12,12,epsabs=1e-14,epsrel=2e-11,limit=200)[0]

def symbolic_checks():
    a,b,c,d=sp.symbols('a b c d',positive=True)
    f=(a+d)/(a+b+c+d); rate=(a*c-b*d)/(a+b+c+d)
    dg,lam,kT=sp.symbols('dg lam kT',positive=True)
    checks={
      'trap_stationarity':sp.simplify((a+d)*(1-f)-(b+c)*f)==0,
      'valence_flux_equals_cycle':sp.simplify(a*(1-f)-b*f-rate)==0,
      'conduction_flux_equals_cycle':sp.simplify(c*f-d*(1-f)-rate)==0,
      'marcus_log_ratio':sp.simplify(-((dg+lam)**2-(-dg+lam)**2)/(4*lam*kT)+dg/kT)==0,
      'depleted_harmonic_mean':sp.simplify(rate.subs({b:0,d:0})-1/(1/a+1/c))==0,
      'single_branch_bound':sp.limit(rate.subs({b:0,d:0}),a,sp.oo)==c}
    assert all(checks.values())
    return checks

def main():
    tests={'symbolic':symbolic_checks()}
    Ts=P['temperature_grid_K']; Fs=np.linspace(0,3,601)
    # Actual recovered FF78 constants are loaded, never rebuilt or overwritten.
    base=json.loads((ROOT/P['baseline_file']).read_text())['full_local_parameters']
    Eg=base['Eg']; er=base['eps_r']; Nt=base['Nt']; N=base['Nc']; cap=base['capture']
    dcm=base['d']; pref=cap*N; fac=Q*Nt*dcm
    # PF: directional upper envelope, emission supply capacity only, not local closure.
    rows=[]
    for T in Ts:
      en0=pref*np.exp(-Eg/2/(KB*T))
      for F in [0,.1,.2,.5,1.,1.5,2.,3.]:
        A=COULOMB/er; fvn=F*.1
        dl=2*np.sqrt(A*fvn); rs=np.sqrt(A/fvn) if F else np.inf
        en=en0*np.exp(min(dl,Eg/2)/(KB*T))
        mid=en*en0/(en+en0)
        # Each point optimizes a different trap; NOT one trap's field trace.
        et=max((Eg-dl)/2,0)
        env=pref/2*np.exp(-max(Eg-dl,0)/2/(KB*T))
        rows.append(dict(T_K=T,F_MVcm=F,delta_pf_eV=dl,saddle_nm=rs,
            midgap_no_field_R_s=en0/2,midgap_singlebranch_R_s=mid,
            midgap_bound_J_Acm2=fac*en0,optimized_Et_from_VB_eV=et,
            optimized_upper_envelope_J_Acm2=fac*env,
            continuum_geometry_at_least_2nm=bool(rs>=2),
            warning='depletion supply upper envelope; continuum core validity unproven'))
    writecsv('pf_supply_bounds.csv',rows)
    tests['pf_midgap_gain_le_2']=all(x['midgap_singlebranch_R_s']<=2*x['midgap_no_field_R_s']*(1+1e-13) for x in rows)
    # Quantum vs classical kernels; same total reorganization .20 eV.
    rows=[]; curves={}; db_errors=[]; mass_errors=[]; trunc_errors=[]
    for T in Ts:
      l,lp=fc_logweights(T)
      mass_errors.append(abs(np.exp(logsumexp(lp))-1))
      mass=np.exp(lp)
      tests[f'fc_nonnegative_{T}']=bool(np.all(mass>=0))
      for dg in np.linspace(-1.2,1.2,97):
        lc=classical_lograte(dg,T); lq=quantum_lograte(dg,T); lg=ground_mlj_lograte(dg,T)
        db=quantum_lograte(dg,T)-quantum_lograte(-dg,T)+dg/(KB*T)
        db_errors.append(abs(float(db)))
        trunc_errors.append(abs(float(quantum_lograte(dg,T,L=80)-quantum_lograte(dg,T,L=100))))
        rows.append(dict(T_K=T,delta_G_eV=dg,classical_k_s=float(np.exp(lc)),
                        quantum_thermal_k_s=float(np.exp(lq)),ground_only_k_s=float(np.exp(lg)),
                        ground_only_DB_error_log=float(lg-ground_mlj_lograte(-dg,T)+dg/(KB*T)),
                        thermal_DB_error_log=float(db)))
      dg=Eg/2-Fs*.1*P['pair_span_nm']/2
      for label,fun in [('Classical Marcus',classical_lograte),('Thermal vibronic',quantum_lograte),('Miller-Abrahams',ma_lograte)]:
        lr=fun(dg,T)-np.log(2)
        curves[label,T]=fac*np.exp(lr)
    writecsv('rate_kernel_comparison.csv',rows)
    tests['fc_weight_mass_max_error']=max(mass_errors)
    tests['quantum_DB_logratio_max_error']=max(db_errors)
    tests['quantum_L80_vs_L100_max_lograte_error']=max(trunc_errors)
    assert max(mass_errors)<5e-13 and max(db_errors)<5e-12 and max(trunc_errors)<1e-12
    rows=[]
    for T in Ts:
      for i,F in enumerate(Fs):
        rows.append(dict(T_K=T,F_MVcm=float(F),delta_leg_eV=float(Eg/2-F*.1*P['pair_span_nm']/2),
          marcus_Jcapacity_Acm2=float(curves['Classical Marcus',T][i]),
          quantum_Jcapacity_Acm2=float(curves['Thermal vibronic',T][i]),
          ma_Jcapacity_Acm2=float(curves['Miller-Abrahams',T][i])))
    writecsv('ideal_pair_capacity.csv',rows)
    # Rate upper bounds, H(distance) and fixed total coupling capacity.
    rows=[]
    for T in Ts:
      for a in [.3,.5,1.]:
        r=P['pair_span_nm']/2; H=.020*np.exp(-(r-1)/a)
        for lam,label in [(.2,'classical'),(.05,'thermal_vibronic_bound')]:
          upper=2*np.pi/HBAR*H*H/np.sqrt(4*np.pi*lam*KB*T)
          rows.append(dict(T_K=T,amplitude_decay_length_nm=a,leg_nm=r,H_at_1nm_eV=.020,
                  H_at_leg_eV=H,kernel=label,leg_upper_rate_s=upper,
                  symmetric_cycle_upper_J_Acm2=fac*upper/2,
                  scope='mathematical kernel upper bound not necessarily attainable or nonadiabatic-valid'))
    writecsv('coupling_capacity_bounds.csv',rows)
    # Ground-only DB failures and Bose/classical validity diagnostics.
    rows=[]
    for T in Ts:
      for hw in [.002,.005,.010,.020,.050,.150,.200]:
        x=hw/(KB*T)
        rows.append(dict(T_K=T,hw_eV=hw,hw_over_kBT=x,nB=1/np.expm1(x),
          thermally_excited_fraction=np.exp(-x),quantum_to_classical_coordinate_variance=(x/2)/np.tanh(x/2)))
    writecsv('vibration_regimes.csv',rows)
    # Equal-chemical-potential equilibrium at nonzero spatial potential drop.
    # Stable log occupations: use exact log sigmoid; avoids 1-f rounded to zero.
    eqerr=[]; thermodynamic=[]
    for T in Ts:
      for F in [.2,1,2.]:
        EV=-Eg/2+F*.1*P['pair_span_nm']/2; EC=-EV; Et=.07; mu=.023
        lFV=-np.logaddexp(0,(EV-mu)/(KB*T)); lHV=-np.logaddexp(0,-(EV-mu)/(KB*T))
        lFC=-np.logaddexp(0,(EC-mu)/(KB*T)); lHC=-np.logaddexp(0,-(EC-mu)/(KB*T))
        la=quantum_lograte(Et-EV,T)+lFV; lb=quantum_lograte(EV-Et,T)+lHV
        lc=quantum_lograte(EC-Et,T)+lHC; ld=quantum_lograte(Et-EC,T)+lFC
        z=np.exp(np.array([la,lb,lc,ld])-max(la,lb,lc,ld)); rr,occ=cycle(*z)
        error=abs(float((la+lc)-(lb+ld))); eqerr.append(error)
        thermodynamic.append(dict(T_K=T,F_MVcm=F,log_cycle_affinity=la+lc-lb-ld,
          trap_occupancy=occ,fermi_expected=expit((mu-Et)/(KB*T)),scaled_net_cycle=rr))
    writecsv('equilibrium_nonzero_field.csv',thermodynamic)
    tests['common_mu_nonzero_field_max_log_affinity']=max(eqerr)
    assert max(eqerr)<5e-12
    # Fixed carrier number versus fixed chemical potential.
    sigma=.0668; frac=1e-4
    mu300=brentq(lambda mu:gaussian_filling(mu,300,sigma)-frac,-1.,.2,xtol=1e-14)
    rows=[]
    for T in Ts:
      mu=brentq(lambda m:gaussian_filling(m,T,sigma)-frac,-1.,.2,xtol=1e-14)
      filling=gaussian_filling(mu300,T,sigma)
      nboltz=np.exp(mu300/(KB*T)+sigma*sigma/(2*(KB*T)**2))
      rows.append(dict(T_K=T,sigma_eV=sigma,fixed_density_fraction=frac,
                      fixed_n_mu_from_DOS_center_eV=mu,fixed_mu_eV=mu300,
                      fixed_mu_density_fraction=filling,
                      naive_boltzmann_fixed_mu_fraction=nboltz,
                      sharp_DOS_ni_cm3=N*np.exp(-Eg/2/(KB*T)),
                      parabolic_DOS_ni_cm3=N*(T/300)**1.5*np.exp(-Eg/2/(KB*T))))
    writecsv('fixed_mu_vs_fixed_n.csv',rows)
    tests['fixed_n_integration_max_abs_error']=max(abs(gaussian_filling(x['fixed_n_mu_from_DOS_center_eV'],x['T_K'])-frac) for x in rows)
    # Transport consistency, not a self-consistent transport calculation.
    rows=[]; mu300base=base['mu_n']; F=1e6; J=.05
    for T in Ts:
      for Eact in [0,.05,.1,.15]:
        mu=mu300base*np.exp(-Eact/KB*(1/T-1/300))
        n=J/(Q*mu*F); tt=dcm/(mu*F)
        rows.append(dict(T_K=T,activation_eV=Eact,mu_cm2Vs=mu,
          transit_time_s=tt,extraction_rate_estimate_s=1/tt,
          carrier_density_needed_for_50mA_cm3=n,needed_occupancy_fraction=n/N,
          unipolar_spacecharge_to_imposed_field=Q*n*dcm/(EPS0_CM*er*F),
          field_Vcm=F,current_Acm2=J,
          warning='illustrative Arrhenius law; not D18:L8-BO mobility data; unipolar estimate'))
    writecsv('transport_capacity.csv',rows)
    # Elastic sequential TAT benchmark: contact-fed transport, NOT body generation.
    pth=(ROOT/P['tat_source_file']); sys.path.insert(0,str(pth.parent))
    from tat_model import tat_state,TATParams
    tatrows=[]
    for T in Ts:
      for F in [0,.05,.1,.2,.5,1,1.5,2,3]:
        V=F*1e6*dcm
        ans=tat_state(V,T,100,TATParams(energy_nodes=320))
        tatrows.append(dict(T_K=T,F_MVcm=F,**ans))
    writecsv('elastic_contact_tat.csv',tatrows)
    tests['elastic_tat_max_conservation_error_Acm2']=max(x['conservation_error_Acm2'] for x in tatrows)
    tests['elastic_tat_zero_bias_max_J_Acm2']=max(abs(x['J_Acm2']) for x in tatrows if x['F_MVcm']==0)
    tatconvergence=[]
    for T in [80,300]:
      for F in [.1,1,3]:
        vals=[tat_state(F*1e6*dcm,T,100,TATParams(energy_nodes=n))['J_Acm2'] for n in [160,320,640]]
        tatconvergence.append(abs(vals[-1]-vals[-2])/max(abs(vals[-1]),1e-300))
    tests['elastic_tat_320_vs_640_relative_max_error']=max(tatconvergence)
    # Nonanalytic WKB turning points slow global Gaussian quadrature convergence.
    assert max(tatconvergence)<1e-7
    # Continuum phonon-assisted escape candidate with phonon exchange at core,
    # then WKB propagation; same spectral coupling used in the reverse process.
    from tat_model import wkb_action
    phrows=[]; pherrs=[]; phdb=[]
    def ph_escape(F,T,rtol=2e-9):
      Et=0.; mu=-.4; ell=2.; drop=F*.1*ell; edge=.25-drop
      l,lp=fc_logweights(T)
      def logP(E):
        return -.5*np.log(4*np.pi*.05*KB*T)+logsumexp(lp-(E-.05-l*.15)**2/(4*.05*KB*T))
      def integrand(eps,reverse=False):
        gamma=1e12*np.exp(-float(wkb_action(.65-eps,.65-drop-eps,ell,.2)))
        return gamma*np.exp(logP(eps if reverse else -eps))*expit((mu-eps if reverse else eps-mu)/(KB*T))
      pts=[v for v in [.65,.65-drop] if edge<v<edge+1.]
      ko,err=quad(lambda e:integrand(e),edge,edge+1.,points=pts,epsabs=1e-50,epsrel=rtol,limit=250)
      ki,erri=quad(lambda e:integrand(e,True),edge,edge+1.,points=pts,epsabs=1e-50,epsrel=rtol,limit=250)
      return ko,ki,abs(np.log(ko)-np.log(ki)-(Et-mu)/(KB*T)),max(err/ko,erri/ki)
    for T in Ts:
      for F in [0,.2,.5,1.,1.5,2.,3.]:
        ko,ki,db,err=ph_escape(F,T)
        phrows.append(dict(T_K=T,F_MVcm=F,escape_rate_s=ko,capture_rate_per_empty_center_s=ki,
          log_detailed_balance_error=db,quadrature_relative_error_estimate=err,
          warning='single trap to one reservoir; no pair cycle; exploratory spectral coupling and band edge'))
        phdb.append(db)
    writecsv('phonon_assisted_continuum.csv',phrows)
    for T in [80,300]:
      for F in [0,1.,3.]:
        k1=ph_escape(F,T,2e-9)[0];k2=ph_escape(F,T,1e-11)[0]
        pherrs.append(abs(k1-k2)/k2)
    tests['phonon_continuum_max_log_DB_error']=max(phdb)
    tests['phonon_continuum_max_refinement_relative_error']=max(pherrs)
    assert max(phdb)<1e-8 and max(pherrs)<1e-7
    # Selected finite-T vs ground-only counterexample.
    cases=[]
    for T in [80,150,300]:
      for dg in [.05,.2,.35,-.6,-1.]:
        lq=quantum_lograte(dg,T); lc=classical_lograte(dg,T); lg=ground_mlj_lograte(dg,T)
        cases.append(dict(T_K=T,delta_eV=dg,k_quantum_s=float(np.exp(lq)),k_classical_s=float(np.exp(lc)),
           quantum_over_classical=float(np.exp(lq-lc)),ground_MLJ_DB_log_defect=float(lg-ground_mlj_lograte(-dg,T)+dg/(KB*T))))
    writecsv('selected_counterexamples.csv',cases)
    # Independent full m,n Franck-Condon sum using Laguerre overlap factors.
    fccheck=[]
    m,n=np.meshgrid(np.arange(64),np.arange(64),indexing='ij')
    small=np.minimum(m,n); large=np.maximum(m,n); diff=np.abs(m-n)
    poly=eval_genlaguerre(small,diff,1.)
    with np.errstate(divide='ignore'):
      lFC=-1.+gammaln(small+1)-gammaln(large+1)+2*np.log(np.abs(poly))
    for T in [80,150,330]:
      lpn=np.log(-np.expm1(-.15/(KB*T)))-n*.15/(KB*T)
      for dg in [-1.2,-.6,-.2,0,.2,.6,1.2]:
        full=np.log(2*np.pi/HBAR)+2*np.log(1e-5)-.5*np.log(4*np.pi*.05*KB*T)+logsumexp(lpn+lFC-(dg+.05+(m-n)*.15)**2/(4*.05*KB*T))
        fccheck.append(abs(float(full-quantum_lograte(dg,T))))
    tests['Bessel_vs_independent_Laguerre_FC_max_lograte_error']=max(fccheck)
    assert max(fccheck)<1e-11
    # Threshold crossing records include every crossing; no endpoint substitution.
    throws=[]
    for T in Ts:
      for label in ['Classical Marcus','Thermal vibronic','Miller-Abrahams']:
        curve=curves[label,T]
        for threshold in [.0001,.001,.005,.050]:
          crossings=[]
          for i in range(len(Fs)-1):
            if (curve[i]-threshold)*(curve[i+1]-threshold)<0:
              # Interpolate in log current; enough for discrete benchmark indicator.
              cross=Fs[i]+(Fs[i+1]-Fs[i])*(np.log(threshold)-np.log(curve[i]))/(np.log(curve[i+1])-np.log(curve[i]))
              crossings.append(float(cross))
          throws.append(dict(T_K=T,model=label,threshold_Acm2=threshold,
             first_crossing_MVcm=crossings[0] if crossings else '',
             all_crossings_MVcm=json.dumps(crossings),crossing_count=len(crossings),
             reason='' if crossings else 'not reached in scan',
             scope='source capacity threshold not measured knee; log-linear interpolation on 0.005 MV/cm grid'))
    writecsv('source_threshold_crossings.csv',throws)
    # Consistent plot units and explicit exploratory caption.
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    colors={80:'#244c8c',150:'#5aa6a8',300:'#d97a22',330:'#a74c65'}
    fig,axs=plt.subplots(1,3,figsize=(14,4.5),sharey=True)
    for ax,label in zip(axs,['Classical Marcus','Thermal vibronic','Miller-Abrahams']):
      for T in [80,150,300,330]: ax.semilogy(Fs,curves[label,T],label=f'{T} K',color=colors[T])
      ax.set(xlabel='Imposed local field (MV/cm)',title=label,ylim=(1e-22,1),xlim=(0,3)); ax.grid(alpha=.2); ax.legend()
    axs[0].set_ylabel('Ideal pair-source capacity (A/cm$^2$)')
    fig.suptitle('Exploratory absorbing-reservoir bound; 10 nm span; not device J-V or a material fit',y=1.02)
    fig.tight_layout();fig.savefig(FIG/'pair_capacity.png',dpi=180,bbox_inches='tight');plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(10,4))
    ds=np.linspace(-1.2,.7,501)
    for T in [80,300]:
      axs[0].semilogy(ds,np.exp(classical_lograte(ds,T)),ls='--',color=colors[T],label=f'classical {T} K')
      axs[0].semilogy(ds,np.exp(quantum_lograte(ds,T)),color=colors[T],label=f'vibronic {T} K')
      er=ground_mlj_lograte(ds,T)-ground_mlj_lograte(-ds,T)+ds/(KB*T)
      axs[1].plot(ds,er,color=colors[T],label=f'ground-only {T} K')
    axs[0].set(xlabel='Elementary free-energy change (eV)',ylabel='Rate (s$^{-1}$)',ylim=(1e-30,1e9));axs[0].legend()
    axs[1].axhline(0,color='k',lw=1);axs[1].set(xlabel='Elementary free-energy change (eV)',ylabel='ln(k forward/k reverse) + deltaG/kBT');axs[1].legend()
    for a in axs:a.grid(alpha=.2)
    fig.suptitle('Exploratory H=10 micro-eV, total lambda=0.20 eV; exact thermal sidebands retain reverse processes')
    fig.tight_layout();fig.savefig(FIG/'quantum_and_balance.png',dpi=180);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(10,4))
    data=list(csv.DictReader((R/'fixed_mu_vs_fixed_n.csv').open()))
    axs[0].semilogy([float(x['T_K']) for x in data],[float(x['fixed_mu_density_fraction']) for x in data],'o-',label='fixed chemical potential')
    axs[0].axhline(frac,ls='--',color='k',label='fixed carrier density');axs[0].legend()
    axs[0].set(xlabel='Temperature (K)',ylabel='Occupied Gaussian DOS fraction',title='Equilibrium DOS illustration, sigma=66.8 meV')
    for Ea in [.05,.1,.15]:
      xx=[x for x in rows if x['activation_eV']==Ea]
      axs[1].semilogy([x['T_K'] for x in xx],[x['unipolar_spacecharge_to_imposed_field'] for x in xx],'o-',label=f'{Ea*1e3:g} meV mobility activation')
    axs[1].axhline(1,ls='--',color='k');axs[1].legend();axs[1].set(xlabel='Temperature (K)',ylabel='Space-charge / imposed-field estimate',title='50 mA/cm2, 1 MV/cm, 100 nm')
    for a in axs:a.grid(alpha=.2)
    fig.suptitle('Illustrations of inconsistent fixed-density / fixed-mobility temperature extrapolations')
    fig.tight_layout();fig.savefig(FIG/'dos_and_extraction.png',dpi=180);plt.close(fig)
    tests['all_assertions_passed']=True
    (R/'validation.json').write_text(json.dumps(tests,indent=2),encoding='utf-8')
    provenance={'created_utc':'2026-10-02','description':'NEW derivations and calculations, not recovered historical calculations; copied inputs identified in inputs/SOURCE_PROVENANCE.json',
      'command':'python audit.py','python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'sympy':sp.__version__,'matplotlib':matplotlib.__version__,
      'inputs':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'parameters.json',(ROOT/P['baseline_file']),(ROOT/P['tat_source_file'])]}}
    (R/'provenance.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    print(json.dumps(tests,indent=2));print('Selected counterexamples:');print(json.dumps(cases,indent=2))

if __name__=='__main__':main()
