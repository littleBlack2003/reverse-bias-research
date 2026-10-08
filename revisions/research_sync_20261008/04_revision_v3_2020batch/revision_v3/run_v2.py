"""run_v2.py : short-circuit / full J-V of closure variants vs original. Uses ORIGINAL solver (device.solve/advance)."""
import sys,json,time
sys.dont_write_bytecode=True
sys.path.insert(0,'code')
import numpy as np
from scipy.optimize import brentq,minimize_scalar
from device import solve,advance,Q
from extend_device import parameters,make_device
from closure_v2 import ClosureDevice

AMP={'cal':(0.47992,0.39802),'one':(1.,1.)}   # (electron, hole) amplitudes; 'cal' = candidate SCLC-trace calibration (density_only)

def saved(T):
    f=np.load(f'reference/critical_states_T{T}_N321.npz',allow_pickle=False);names=[str(n) for n in f['names']]
    return {n:(float(l[1]),float(l[0]),z.copy()) for n,l,z in zip(names,f['labels'],f['z'])}

def build(T,amp,lam_d,lam_b,zeta=None,N=321):
    p,_,_=parameters(float(T));a=AMP[amp] if amp else (1.,1.)
    d=ClosureDevice(p,N,amp_n=a[0],amp_p=a[1],lam_d=lam_d,lam_b=lam_b,zeta=zeta)
    return d

def continue_to(d,z,V,L,final,mn_base,mp_base,amp,t0=0.,dt=.25,depth=0,log=None):
    """homotopy t: amp^t, lam=t*final; adaptive halving"""
    t=t0;lam_d,lam_b=final
    while t<1-1e-12:
        step=min(dt,1-t);t1=t+step
        d.mn=mn_base*amp[0]**t1;d.mp=mp_base*amp[1]**t1;d.lam_d=lam_d*t1;d.lam_b=lam_b*t1
        zz,info=solve(d,V,z,L,maxiter=300)
        if info['success']:z=zz;t=t1;dt=min(dt*1.5,.5)
        else:
            dt=step/2
            if dt<1e-4:raise RuntimeError(f'homotopy stalled at t={t:.4f}: {info.get("reason")}')
    return z

def solve_state(T,amp,lam_d,lam_b,zeta=None,V=0.,start=None):
    d=build(T,amp,0.,0.,zeta);mn_base,mp_base=d.mn/ (AMP[amp][0] if amp else 1),d.mp/(AMP[amp][1] if amp else 1)
    # base here is un-amplified
    d2=build(T,None,0.,0.,zeta);mn_base,mp_base=d2.mn,d2.mp
    d=build(T,None,0.,0.,zeta);d.zeta=zeta
    z0=saved(T)['sc'][2] if start is None else start
    z=continue_to(d,z0,V,1.,(lam_d,lam_b),mn_base,mp_base,AMP[amp] if amp else (1.,1.),dt=.25)
    return d,z

def metrics_at(d,z,V=0.):
    o=d.evaluate(z,V);J=np.mean(o['Jn']+o['Jp']);led=d.ledger(z,V,1.)
    bulk=(d.x>.2)&(d.x<.8)
    return dict(Jsc_mAcm2=-J*1e3,eta=-J/(Q*d.p.G_cm3s*d.p.d_cm),R_over_G=led['recombination_Acm2']/led['generation_Acm2'],
        n_bulk_median=float(np.median(o['n'][bulk])),p_bulk_median=float(np.median(o['p'][bulk])),
        mu_n_bulk_med=float(np.median((d.mn*o['mobility_factor_n'])[bulk[:-1]])),mu_p_bulk_med=float(np.median((d.mp*o['mobility_factor_p'])[bulk[:-1]])),
        gate=bool(led['gate_passed']),cont_rel=led['continuity_relative'])

def calibrate_zeta(amp):
    """zeta chosen so that Langevin-tied beta equals the 300 K calibrated beta in the bulk (median, x in 0.2..0.8) of the 300 K SC state."""
    d,z=solve_state(300,amp,1.,0.)
    mun,mup=d.mu_node(z,0.);bulk=(d.x>.2)&(d.x<.8)
    bl=np.median(Q*(mun+mup)[bulk]/d.eps)
    return d.p.beta_cm3s/bl
