"""Gaussian-DOS transport closure; cm,eV,V,s units; isothermal local equilibrium.
Primary definitions: Pasveer et al. PRL 94,206601 Eqs3-5; Doan et al.
WIAS2493 (2018) Eq2.3,2.11 and inverse-activity flux Sec4.1.
Not a calibrated D18/L8-BO model; validity metadata is mandatory output.
"""
from dataclasses import dataclass
from functools import lru_cache
import numpy as np
from scipy.special import roots_hermitenorm, expit
KB=8.617333262145e-5
Q=1.602176634e-19
@lru_cache(None)
def normal_rule(order):
    x,w=roots_hermitenorm(order)
    return x,w/np.sqrt(2*np.pi)

def gaussian_eos(eta,s,order=256):
    """Return occupation c=n/N, dc/deta, Einstein factor c/(dc/deta).
    eta=(electron EF-local LUMO center)/kT; hole eta has opposite energy sign.
    s=sigma/kT. Independent Gaussian energies and Fermi occupations.
    """
    if not np.isfinite(s) or not 0<=s<=6: raise ValueError('EOS validated only for 0<=sigma/kT<=6')
    eta=np.asarray(eta,dtype=float)
    if not np.all(np.isfinite(eta)):raise ValueError('nonfinite chemical potential')
    x,w=normal_rule(order)
    arg=eta[...,None]-s*x
    f=expit(arg);fb=expit(-arg)
    c=f@w;dc=(f*fb)@w
    g3=np.divide(c,dc,out=np.full_like(c,np.inf),where=dc>0)
    return c,dc,g3

@dataclass(frozen=True)
class EGDM:
    sigma_eV:float
    a_cm:float
    prefactor_cm2_Vs:float # mu0*c1 in Pasveer; independent of T
    def __post_init__(self):
        if not np.all(np.isfinite([self.sigma_eV,self.a_cm,self.prefactor_cm2_Vs])) or min(self.sigma_eV,self.a_cm,self.prefactor_cm2_Vs)<=0:raise ValueError('positive material parameters required')
    @property
    def N_cm3(self):return self.a_cm**-3
    def evaluate(self,T,n_cm3,F_Vcm,strict=False):
        if not np.isfinite(T) or T<=0:raise ValueError('T>0 required')
        n,F=np.broadcast_arrays(np.asarray(n_cm3,dtype=float),np.asarray(F_Vcm,dtype=float))
        if not np.all(np.isfinite(n)) or not np.all(np.isfinite(F)):raise ValueError('nonfinite density or field')
        if np.any(n<0) or np.any(n>self.N_cm3*(1+1e-14)):raise ValueError('density outside [0,N]')
        s=self.sigma_eV/(KB*T);c=n*self.a_cm**3;u=abs(F)*self.a_cm/self.sigma_eV
        if s<=1:raise ValueError('Pasveer density formula undefined for sigma/kT<=1')
        delta=2*(np.log(s*s-s)-np.log(np.log(4)))/(s*s)
        lng1=.5*(s*s-s)*(2*c)**delta
        lng2=.44*(s**1.5-2.2)*(np.sqrt(1+.8*u*u)-1)
        lnmu=np.log(self.prefactor_cm2_Vs)-.42*s*s+lng1+lng2
        mask=(s>=2)&(s<=6)&(c>=1e-6)&(c<=1e-2)&(u<=3)
        # u<=3 is ORIGINAL FIGURE COVERAGE, not a universal proven validity bound.
        if strict and not np.all(mask):raise ValueError('outside declared literature comparison window')
        return dict(mu=np.exp(lnmu),lnmu=lnmu,s=s,occupation=c,reduced_field=u,
                    density_benchmark=(c>=1e-6)&(c<=1e-2),temperature_benchmark=2<=s<=6,
                    field_plot_coverage=u<=3,comparison_window=mask,
                    low_density_asymptote=c<1e-6)

def bernoulli(x):
    x=np.asarray(x,dtype=float)
    out=np.empty_like(x);small=abs(x)<1e-4;pos=x>50;neg=x<-50;mid=~(small|pos|neg)
    out[small]=1-x[small]/2+x[small]**2/12-x[small]**4/720
    out[pos]=x[pos]*np.exp(-x[pos])/(1-np.exp(-x[pos]))
    out[neg]=-x[neg]/(1-np.exp(x[neg]));out[mid]=x[mid]/np.expm1(x[mid])
    return out

def inverse_activity_flux(eta,psi,T,N_cm3,sigma_eV,dx_cm,mu_face,carrier,affinity=None,order=256):
    """Conventional current at each face; psi=electrostatic potential/kT[V].
    Positive q, separate hole occupation convention. Supplying exact dimensionless
    quasi-Fermi difference affinity avoids cancellation near equilibrium.
    n: affinity=deta-dpsi; p: affinity=deta+dpsi.
    This replaces single-level vacancy-blocked SG; it is not the same flux at
    finite occupation even when sigma tends to zero (mobility conventions differ).
    """
    eta=np.asarray(eta,dtype=float);psi=np.asarray(psi,dtype=float)
    if eta.ndim!=1 or len(eta)<2 or psi.shape!=eta.shape:raise ValueError('matching one-dimensional node arrays required')
    if not np.all(np.isfinite(eta)) or not np.all(np.isfinite(psi)):raise ValueError('nonfinite state')
    if not np.all(np.isfinite([T,N_cm3,dx_cm,sigma_eV])) or min(T,N_cm3,dx_cm)<=0 or sigma_eV<0:raise ValueError('invalid flux parameters')
    if not np.all(np.isfinite(mu_face)) or np.any(np.asarray(mu_face)<=0):raise ValueError('positive finite mobility required')
    dp=np.diff(psi);de=np.diff(eta);mid=(eta[:-1]+eta[1:])/2
    c=gaussian_eos(mid,sigma_eV/(KB*T),order)[0]
    if carrier=='n':sign=1;arg=-dp;aff=de-dp if affinity is None else affinity
    elif carrier=='p':sign=-1;arg=dp;aff=de+dp if affinity is None else affinity
    else:raise ValueError('carrier n or p')
    aff=np.asarray(aff)
    if aff.shape not in [(),de.shape] or not np.all(np.isfinite(aff)):raise ValueError('invalid face affinity')
    if np.any(abs(de)>500) or np.any(abs(aff)>500):raise ValueError('Newton trial gradient outside safe flux evaluation range; damp/refine')
    # Evaluate the entire signed product in logs; individual factors can overflow
    # even while the physical current is representable. Exact zero affinity stays zero.
    logb=np.empty_like(arg);pos=arg>50;neg=arg<-50;midmask=~(pos|neg)
    logb[pos]=np.log(arg[pos])-arg[pos]-np.log1p(-np.exp(-arg[pos]))
    logb[neg]=np.log(-arg[neg])-np.log1p(-np.exp(arg[neg]))
    logb[midmask]=np.log(bernoulli(arg[midmask]))
    aa=np.broadcast_to(aff,de.shape);ab=abs(aa);nonzero=ab>0
    le=np.full_like(aa,-np.inf)
    le[nonzero]=np.maximum(aa[nonzero],0)+np.log(-np.expm1(-ab[nonzero]))
    with np.errstate(divide='ignore'):
        lc=np.log(c)
    lc=np.where(c>0,lc,mid+.5*(sigma_eV/(KB*T))**2)
    logabs=np.log(Q)+np.log(mu_face)+np.log(KB*T)+np.log(N_cm3)+lc-de/2+logb+le-np.log(dx_cm)
    if np.any(logabs>np.log(np.finfo(float).max)):raise FloatingPointError('flux exceeds floating range; damp or refine trial')
    return sign*np.sign(aa)*np.exp(logabs)
