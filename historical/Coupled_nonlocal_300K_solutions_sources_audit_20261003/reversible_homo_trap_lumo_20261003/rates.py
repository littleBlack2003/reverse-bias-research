"""Conditional two-charge-state HOMO/trap/LUMO spectral transfer. eV, cm, s.
No material calibration or device solver. Finite, real-valued input only.
"""
from dataclasses import dataclass
import math
import numpy as np
from scipy.special import logsumexp, expit
KB=8.617333262145e-5 # eV/K
HBAR=6.582119569e-16 # eV s
Q=1.602176634e-19 # C

def _positive(**v):
    if any(not math.isfinite(x) or x<=0 for x in v.values()):
        raise ValueError('finite positive inputs required: '+str(v))

def log_kernel(delta, t0, distance, xi, lam, temperature=300.):
    """Nonadiabatic classical Marcus; t0 eV; distance and xi same length units."""
    _positive(t0=t0,xi=xi,lam=lam,temperature=temperature)
    if not math.isfinite(distance) or distance<0: raise ValueError('distance')
    delta=np.asarray(delta,dtype=float)
    if not np.all(np.isfinite(delta)): raise ValueError('delta')
    kt=KB*temperature
    return (math.log(2*math.pi/HBAR)+2*math.log(t0)-2*distance/xi
            -.5*math.log(4*math.pi*lam*kt)-(delta+lam)**2/(4*lam*kt))

def log_exchange(energies, weights, mu, Et, *, t0, distance, xi, lam, temperature=300.):
    """weights = spectral DOS*quadrature dE, dimensionless states per trap.
    A physical bulk DOS needs a separately justified coupling-volume/spatial integral.
    Returns log(a),log(b), each rate in s^-1. mu and all energies in common eV gauge.
    """
    e=np.asarray(energies,dtype=float); w=np.asarray(weights,dtype=float)
    if e.ndim!=1 or e.shape!=w.shape or not np.all(np.isfinite(e)) or not np.all(np.isfinite(w)) or np.any(w<0) or not np.any(w>0): raise ValueError('quadrature')
    if not all(math.isfinite(v) for v in (mu,Et)): raise ValueError('energies')
    _positive(temperature=temperature)
    sel=w>0;e=e[sel];w=w[sel];x=(e-mu)/(KB*temperature)
    lf=-np.logaddexp(0,x);lv=-np.logaddexp(0,-x)
    args=dict(t0=t0,distance=distance,xi=xi,lam=lam,temperature=temperature)
    return (float(logsumexp(np.log(w)+lf+log_kernel(Et-e,**args))),
            float(logsumexp(np.log(w)+lv+log_kernel(e-Et,**args))))

@dataclass(frozen=True)
class FourRates:
    laH: float
    lbH: float
    laL: float
    lbL: float
    def __post_init__(self):
        if not all(math.isfinite(x) for x in (self.laH,self.lbH,self.laL,self.lbL)):
            raise ValueError('finite log rates required')
    @property
    def logD(self): return float(logsumexp([self.laH,self.lbH,self.laL,self.lbL]))
    @property
    def occupancy(self):
        return float(expit(np.logaddexp(self.laH,self.laL)-np.logaddexp(self.lbH,self.lbL)))
    @property
    def affinity(self): return self.laH+self.lbL-self.laL-self.lbH
    def generation(self,Nt=1.):
        """Positive generation, stable near-equilibrium steady state."""
        _positive(Nt=Nt)
        u=self.laH+self.lbL;v=self.laL+self.lbH
        if u==v:return 0.
        logabs=math.log(Nt)+max(u,v)-self.logD+math.log(-math.expm1(-abs(u-v)))
        return math.copysign(math.exp(logabs),u-v)
    def events(self,f):
        if not math.isfinite(f) or not 0<=f<=1:raise ValueError('f')
        return dict(hole_emission=math.exp(self.laH)*(1-f),hole_capture=math.exp(self.lbH)*f,
                    electron_capture=math.exp(self.laL)*(1-f),electron_emission=math.exp(self.lbL)*f)
    def sources(self,f,Nt=1.,z_empty=0.):
        _positive(Nt=Nt);e=self.events(f)
        jH=e['hole_emission']-e['hole_capture'];jL=e['electron_emission']-e['electron_capture']
        return dict(Sn=Nt*jL,Sp=Nt*jH,df=jH-jL,rho=Q*Nt*(z_empty-f),drho=-Q*Nt*(jH-jL))
    def advance(self,f,dt):
        if not math.isfinite(f) or not 0<=f<=1 or not math.isfinite(dt) or dt<0:raise ValueError('step')
        if dt==0: return f
        logstep=self.logD+math.log(dt)
        if logstep>700: return self.occupancy
        step=math.exp(logstep)
        return f*math.exp(-step)+self.occupancy*(-math.expm1(-step))

def from_positive(aH,bH,aL,bL):
    _positive(aH=aH,bH=bH,aL=aL,bL=bL)
    return FourRates(*map(math.log,(aH,bH,aL,bL)))

def electron_energy(E0,phi):return np.asarray(E0)-phi

def path_current_faces(x, xH, xt, xL, jH, jL, areal_traps=1.):
    """Conventional current A/cm² for 1D oriented point-to-point events.
    x is face coordinates, matching all positions' length unit. Event j in s^-1.
    An electron H->t contributes -q*jH; t->L contributes -q*jL, times direction.
    areal_traps is traps/cm² represented by the reaction element.
    Faces must not coincide with point-event positions (place events in cells).
    """
    x=np.asarray(x,dtype=float)
    if x.ndim!=1 or not np.all(np.isfinite(x)) or not all(math.isfinite(v) for v in (xH,xt,xL,jH,jL,areal_traps)) or areal_traps<0:
        raise ValueError('finite path inputs and nonnegative areal_traps required')
    if any(np.any(x==p) for p in (xH,xt,xL)):
        raise ValueError('point events must be inside cells, not on faces')
    def path(a,b,j):
        if a==b:return np.zeros_like(x)
        return -Q*areal_traps*j*np.sign(b-a)*((x>min(a,b))&(x<max(a,b)))
    return path(xH,xt,jH)+path(xt,xL,jL)
