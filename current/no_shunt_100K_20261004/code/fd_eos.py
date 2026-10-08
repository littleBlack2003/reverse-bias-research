"""Gaussian Fermi–Dirac EOS for a numerical, local-equilibrium Y6 baseline.

The EOS is exact as an integral for every finite disorder; the s <= 8.7 guard
is ONLY the independently checked numerical range, not an EGDM mobility bound.
801-point Gauss–Legendre on standard-normal z in [-16,16] generates the
negative-eta half of a dense interpolation table. Positive eta uses exact
particle–hole symmetry, avoiding rounded occupations as interpolation knots.
Interpolation is C1 cubic Hermite with quadrature derivatives. Its derivative
is used consistently in the generalized Einstein factor. No material fit or
mobility model enters this module. Units: sigma in eV, T in K.
"""
from functools import lru_cache
import numpy as np
from scipy.special import roots_legendre, expit
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator

KB = 8.617333262145e-5
VALIDATED_S_MAX = 8.7
TABLE_ETA_MIN = -240.0
TABLE_STEP = 0.02
QUADRATURE_ORDER = 801
QUADRATURE_BOUND = 16.0

@lru_cache(maxsize=16)
def _rule(order=QUADRATURE_ORDER, bound=QUADRATURE_BOUND):
    z,w = roots_legendre(order)
    z = z*bound
    w = w*bound*np.exp(-0.5*z*z)/np.sqrt(2*np.pi)
    return z,w

def direct_eos(eta,s,order=QUADRATURE_ORDER,bound=QUADRATURE_BOUND):
    """Direct quadrature c, dc/deta; primarily a table-generation/reference tool."""
    eta=np.asarray(eta,dtype=float)
    if not np.isfinite(s) or not 0 <= s <= VALIDATED_S_MAX:
        raise ValueError('EOS numerical validation requires 0 <= sigma/kT <= 8.7')
    if not np.all(np.isfinite(eta)): raise ValueError('finite eta required')
    if s == 0:
        c=expit(eta); return c,c*expit(-eta)
    z,w=_rule(order,bound)
    # Work in negative half, where neither c nor dc suffers 1-c cancellation.
    flat=eta.ravel(); c=np.empty_like(flat); dc=np.empty_like(flat)
    for start in range(0,len(flat),256):
        sl=slice(start,start+256); en=-np.abs(flat[sl])
        arg=en[:,None]-s*z
        f=expit(arg)
        cn=np.sum(f*w,axis=1)
        dn=np.sum(f*expit(-arg)*w,axis=1)
        c[sl]=np.where(flat[sl]>0,1-cn,cn); dc[sl]=dn
    return c.reshape(eta.shape),dc.reshape(eta.shape)

@lru_cache(maxsize=32)
def _table(s):
    eta=np.linspace(TABLE_ETA_MIN,0,round(-TABLE_ETA_MIN/TABLE_STEP)+1)
    c,dc=direct_eos(eta,s)
    logc=np.log(c)
    dlogc=dc/c
    # Symmetry makes this value known exactly; enforce it despite rounding.
    logc[-1]=-np.log(2.)
    f=CubicHermiteSpline(eta,logc,dlogc,extrapolate=False)
    # Monotone inverse seed; Newton then inverts the very same smooth EOS.
    inv=PchipInterpolator(logc,eta,extrapolate=False)
    return f,f.derivative(),inv,float(logc[0])

class GaussianDOS:
    def __init__(self,sigma_eV,T):
        if not np.all(np.isfinite([sigma_eV,T])) or sigma_eV<0 or T<=0:
            raise ValueError('sigma >= 0 and T > 0 required')
        self.sigma_eV=float(sigma_eV); self.T=float(T); self.s=self.sigma_eV/(KB*self.T)
        if self.s>VALIDATED_S_MAX:
            raise ValueError('EOS numerical validation requires sigma/kT <= 8.7')
        if self.s:
            self._f,self._df,self._inv,self._min_logc=_table(self.s)

    def _negative_logc_derivative(self,eta):
        # eta must be <= 0; tails below table are far inside Boltzmann regime:
        # second/first moment bound exp(eta+1.5*s*s)<2e-55 at the s=8.7 limit.
        eta=np.asarray(eta,dtype=float)
        if self.s==0:
            return -np.logaddexp(0.,-eta),expit(-eta)
        inside=np.maximum(eta,TABLE_ETA_MIN)
        lc=self._f(inside); d=np.minimum(1.,self._df(inside))
        # Enforce the exact bound d(log c)/deta <= 1 against ~1e-12
        # interpolation roundoff in the asymptotic linear tail.
        low=eta<TABLE_ETA_MIN
        return np.where(low,eta+0.5*self.s*self.s,lc),np.where(low,1.,d)

    def logc(self,eta):
        eta=np.asarray(eta,dtype=float)
        if not np.all(np.isfinite(eta)): raise ValueError('finite eta required')
        lc,_=self._negative_logc_derivative(-np.abs(eta))
        return np.where(eta>0,np.log1p(-np.exp(lc)),lc)

    def dlogc(self,eta):
        eta=np.asarray(eta,dtype=float)
        if not np.all(np.isfinite(eta)): raise ValueError('finite eta required')
        lc,d=self._negative_logc_derivative(-np.abs(eta))
        # In the positive half d(log c)=c_negative/c_positive*d(log c_negative).
        ratio=np.exp(lc)/(-np.expm1(lc))
        return np.where(eta>0,ratio*d,d)

    def evaluate(self,eta):
        lc=self.logc(eta); dlc=self.dlogc(eta)
        c=np.exp(lc); dc=c*dlc
        g3=np.divide(1.,dlc,out=np.full_like(dlc,np.inf),where=dlc>0)
        return c,dc,g3

    def eta_from_logc(self,logc):
        logc=np.asarray(logc,dtype=float)
        if not np.all(np.isfinite(logc)) or np.any(logc>=0):
            raise ValueError('finite log occupation strictly below zero required')
        upper=logc>-np.log(2.)
        # log(1-exp(logc)) remains stable even for logc close to zero.
        small_logc=np.where(upper,np.log(-np.expm1(logc)),logc)
        if self.s==0:
            eta=small_logc-np.log(-np.expm1(small_logc))
        else:
            lc_clip=np.clip(small_logc,self._min_logc,-np.log(2.))
            eta=np.minimum(0.,self._inv(lc_clip))
            for _ in range(2):
                lc,dlc=self._negative_logc_derivative(eta)
                eta=np.minimum(0.,eta-(lc-lc_clip)/dlc)
            eta=np.where(small_logc<self._min_logc,
                         small_logc-0.5*self.s*self.s,eta)
        return np.where(upper,-eta,eta)

    def eta_from_c(self,c):
        c=np.asarray(c,dtype=float)
        if not np.all(np.isfinite(c)) or np.any(c<=0) or np.any(c>=1):
            raise ValueError('occupation must lie strictly between zero and one')
        return self.eta_from_logc(np.log(c))

    def secant_g(self,eta_left,eta_right):
        """Positive secant Einstein factor Δη/Δln(c), with continuous equal-state limit."""
        left,right=np.broadcast_arrays(np.asarray(eta_left,float),np.asarray(eta_right,float))
        de=right-left; dlc=self.logc(right)-self.logc(left)
        small=np.abs(de)<1e-5
        mid=0.5*(right+left)
        g=np.divide(de,dlc,out=np.ones_like(de),where=~small)
        return np.where(small,1./self.dlogc(mid),g)

# Compatibility wrapper for existing standalone closure tests.
def gaussian_eos(eta,s,order=None):
    if order is not None:
        c,dc=direct_eos(eta,s,order=order)
        return c,dc,np.divide(c,dc,out=np.full_like(c,np.inf),where=dc>0)
    return GaussianDOS(s*KB,1.).evaluate(eta)
