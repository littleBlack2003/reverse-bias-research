"""Versioned finite-filling mobility candidate, not a 100 K material validation.
Pasveer PRL94 206601(2005) Eqs3-6; sigmahat fitted range2..6.
At c>.1 use explicit conventional .1 saturation, flagged and sensitivity-tested.
No extra generalized Einstein factor: the base solver already includes FD DOS.
"""
import numpy as np
from device import KB, mobility

def density_factor(c,sigma,T,c_cut=.1,extrapolate=False):
    s=sigma/(KB*T)
    if not extrapolate and not 2<=s<=6:raise ValueError(f'EGDM disorder ratio {s:g} outside published2..6 range')
    if s<=1:raise ValueError('EGDM density expression requires sigma/kT>1')
    c=np.asarray(c)
    if np.any((c<0)|(c>1)):raise ValueError('Filling must be in[0,1]')
    delta=2*(np.log(s*s-s)-np.log(np.log(4)))/(s*s)
    return np.exp(.5*(s*s-s)*(2*np.minimum(c,c_cut))**delta)

def field_factor(F,sigma,T,a_cm):
    s=sigma/(KB*T);f=np.asarray(F)*a_cm/sigma
    return np.exp(.44*(s**1.5-2.2)*(np.sqrt(1+.8*f*f)-1))

def dilute_shape(sigma,T,mode):
    if mode=='egdm':return np.exp(-.42*(sigma/KB)**2*(1/T**2-1/300**2))
    if mode=='density_only':return mobility(1.,sigma,T)
    if mode=='baseline':return mobility(1.,sigma,T)
    raise ValueError(mode)

def mobility_factor(c,F,sigma,T,a_cm,mode,gamma=0.,c_cut=.1,field_cap=235000.,extrapolate=False):
    F=np.minimum(abs(np.asarray(F)),field_cap)
    if mode in ('density_only','baseline'):
        field=np.exp(gamma*((F*F+1)**.25-1))
    elif mode=='egdm':field=field_factor(F,sigma,T,a_cm)
    else:raise ValueError(mode)
    return field*(1 if mode=='baseline' else density_factor(c,sigma,T,c_cut,extrapolate))
