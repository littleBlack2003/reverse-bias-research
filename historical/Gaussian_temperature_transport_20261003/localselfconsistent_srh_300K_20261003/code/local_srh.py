"""Conditional local, nondegenerate SRH interface. cm, s, eV; no device fit."""
from dataclasses import dataclass
import math
KB_EV=8.617333262145e-5
Q_E=1.602176634e-19

@dataclass(frozen=True)
class LocalSRH:
    cn: float
    cp: float
    n1: float
    p1: float
    Nt: float=1.0
    z_empty: float=0.0
    temperature: float=300.0
    def __post_init__(self):
        for name in ('cn','cp','n1','p1','Nt','temperature'):
            v=getattr(self,name)
            if not math.isfinite(v) or v<=0: raise ValueError(name+' must be finite and positive')
    @property
    def en(self): return self.cn*self.n1
    @property
    def ep(self): return self.cp*self.p1
    def inputs(self,n,p):
        if any(not math.isfinite(x) or x<0 for x in (n,p)): raise ValueError('n,p must be finite and nonnegative')
    def occupancy(self,n,p):
        self.inputs(n,p)
        a=self.cn*n+self.ep; b=self.en+self.cp*p
        return a/(a+b)
    def events(self,n,p,f):
        self.inputs(n,p)
        if not 0<=f<=1: raise ValueError('f outside [0,1]')
        return dict(electron_capture=self.cn*n*(1-f),electron_emission=self.en*f,
                    hole_capture=self.cp*p*f,hole_emission=self.ep*(1-f))
    def sources(self,n,p,f):
        r=self.events(n,p,f)
        Sn=self.Nt*(r['electron_emission']-r['electron_capture'])
        Sp=self.Nt*(r['hole_emission']-r['hole_capture'])
        df=(Sp-Sn)/self.Nt
        return dict(Sn=Sn,Sp=Sp,df=df,rho=Q_E*self.Nt*(self.z_empty-f),
                    drho=-Q_E*self.Nt*df)
    def U(self,n,p):
        """Positive recombination. Stable enough for declared finite test envelope."""
        self.inputs(n,p)
        return self.Nt*self.cn*self.cp*(n*p-self.n1*self.p1)/(self.cn*(n+self.n1)+self.cp*(p+self.p1))
    def steady_sources(self,n,p):
        """Avoid subtracting nearly equal gross fluxes at stationary occupancy."""
        U=self.U(n,p); f=self.occupancy(n,p)
        return dict(Sn=-U,Sp=-U,df=0.,rho=Q_E*self.Nt*(self.z_empty-f),drho=0.)
    def advance(self,n,p,f,dt):
        """Exact fixed-reservoir evolution, not an explicit Euler approximation."""
        self.inputs(n,p)
        if not 0<=f<=1 or dt<0: raise ValueError('invalid initial f or dt')
        fs=self.occupancy(n,p)
        return fs+(f-fs)*math.exp(-(self.cn*n+self.cp*p+self.en+self.ep)*dt)
    def enhance(self,gn=1.0,gp=1.0):
        if not all(math.isfinite(g) and g>0 for g in (gn,gp)): raise ValueError('invalid enhancement')
        return LocalSRH(self.cn*gn,self.cp*gp,self.n1,self.p1,self.Nt,self.z_empty,self.temperature)

# DA state labels (D electron occupancy, A electron occupancy); charge/q=1-D-A.
DA_STATES=((1,0),(0,1),(0,0),(1,1))
# All rows oriented in pair-generation direction. dn,dp count free carriers produced.
DA_EDGES=((0,1,0,0,'internal'),(1,2,1,0,'electron'),(2,0,0,1,'hole'),
          (1,3,0,1,'hole'),(3,0,1,0,'electron'))
def da_equilibrium_edges(energies,mu=0.,temperature=300.,attempt=1.):
    """Toy reversible edge rates only; energies and attempt are NOT material parameters."""
    beta=1/(KB_EV*temperature)
    w=[E-mu*sum(s) for E,s in zip(energies,DA_STATES)]
    m=min(w); weights=[math.exp(-(v-m)*beta) for v in w]; z=sum(weights)
    pi=[a/z for a in weights]
    rates=[]
    for i,j,dn,dp,kind in DA_EDGES:
        x=(w[j]-w[i])*beta/2
        rates.append((attempt*math.exp(-x),attempt*math.exp(x)))
    return pi,rates

def da_sources(probabilities,rates):
    if len(probabilities)!=4 or len(rates)!=5: raise ValueError('DA dimensions')
    dP=[0.]*4; Sn=Sp=0.
    for (i,j,dn,dp,_),(kf,kr) in zip(DA_EDGES,rates):
        flux=probabilities[i]*kf-probabilities[j]*kr
        dP[i]-=flux; dP[j]+=flux; Sn+=dn*flux; Sp+=dp*flux
    dz=sum((1-sum(s))*dp for s,dp in zip(DA_STATES,dP))
    return dict(dP=dP,Sn=Sn,Sp=Sp,dz=dz)
