"""Conservative 1D HOMO/trap/LUMO reaction block. cm, eV, V, s, C.
Conditional localized sequential Marcus candidate; no material identification.
No drift/diffusion, contacts or Poisson solve is hidden in this module.
"""
from dataclasses import dataclass
import numpy as np
from scipy.sparse import csr_matrix, lil_matrix
from rates import FourRates, log_exchange, Q

@dataclass(frozen=True)
class Mesh:
    faces: np.ndarray
    nodes: np.ndarray = None
    def __post_init__(self):
        e=np.asarray(self.faces,dtype=float).copy()
        if e.ndim!=1 or len(e)<3 or not np.all(np.isfinite(e)) or np.any(np.diff(e)<=0): raise ValueError('strictly increasing finite faces')
        x=(e[:-1]+e[1:])/2 if self.nodes is None else np.asarray(self.nodes,dtype=float).copy()
        if x.shape!=(len(e)-1,) or not np.all(np.isfinite(x)) or np.any(x<=e[:-1]) or np.any(x>=e[1:]): raise ValueError('one node strictly inside each control volume')
        e.setflags(write=False);x.setflags(write=False)
        object.__setattr__(self,'faces',e);object.__setattr__(self,'nodes',x)
    @property
    def volumes(self): return np.diff(self.faces)
    @property
    def size(self): return len(self.nodes)
    def shape(self,x):
        # No endpoint clipping/extrapolation: reservoirs must be resolved by nodes.
        if not np.isfinite(x) or x<self.nodes[0] or x>self.nodes[-1]: raise ValueError('endpoint outside resolved node interval')
        w=np.zeros(self.size); k=np.searchsorted(self.nodes,x)
        if k<self.size and self.nodes[k]==x: w[k]=1
        else:
            a=(x-self.nodes[k-1])/(self.nodes[k]-self.nodes[k-1]);w[k-1]=1-a;w[k]=a
        return w

@dataclass(frozen=True)
class Reaction:
    xH: float
    xt: float
    xL: float
    areal_traps: float # Nt*dX, NOT a volume density
    kernel: object
    z_empty: float=0.5 # explicit compatibility reference, not identified trap valence
    def __post_init__(self):
        if not np.all(np.isfinite([self.xH,self.xt,self.xL,self.areal_traps,self.z_empty])) or self.areal_traps<=0: raise ValueError('reaction geometry/count')

@dataclass(frozen=True)
class FixedRates:
    rates: FourRates
    def __call__(self,ports,reaction): return self.rates

@dataclass(frozen=True)
class SpectralRates:
    EH0: np.ndarray
    wH: np.ndarray
    EL0: np.ndarray
    wL: np.ndarray
    Et0: float
    H: dict # t0, xi, lam (no distance/temperature)
    L: dict
    temperature: float=300.
    def __call__(self,ports,reaction):
        ph,pt,pl,muh,mul=ports
        et=self.Et0-pt
        ah,bh=log_exchange(np.asarray(self.EH0)-ph,self.wH,muh,et,distance=abs(reaction.xt-reaction.xH),temperature=self.temperature,**self.H)
        al,bl=log_exchange(np.asarray(self.EL0)-pl,self.wL,mul,et,distance=abs(reaction.xt-reaction.xL),temperature=self.temperature,**self.L)
        return FourRates(ah,bh,al,bl)

class Adapter:
    """State layout [phi(N), muH(N), muL(N), f(M)].
    Residual block layout [Sn(N), Sp(N), rho_t(N), df(M)].
    muH is the HOMO electron chemical potential (= E_Fp), not a hole potential.
    Endpoint interpolation and source deposition use identical shape weights.
    """
    def __init__(self,mesh,reactions):
        self.mesh=mesh;self.reactions=tuple(reactions);self.N=mesh.size;self.M=len(self.reactions)
        if not self.M: raise ValueError('at least one reaction')
        self.W=[tuple(mesh.shape(x) for x in (r.xH,r.xt,r.xL)) for r in self.reactions]
        self.D=[tuple(w/mesh.volumes for w in ws) for ws in self.W]
    def unpack(self,state):
        s=np.asarray(state,dtype=float)
        if s.shape!=(3*self.N+self.M,) or not np.all(np.isfinite(s)): raise ValueError('state')
        phi,muh,mul=np.split(s[:3*self.N],3);f=s[3*self.N:]
        if np.any((f<0)|(f>1)): raise ValueError('occupancy')
        return phi,muh,mul,f
    def ports(self,state,k):
        phi,muh,mul,f=self.unpack(state);h,t,l=self.W[k]
        return np.array([h@phi,t@phi,l@phi,h@muh,l@mul])
    def _assemble(self,k,rate,f,steady=False):
        r=self.reactions[k];dh,dt,dl=self.D[k];wh,wt,wl=self.W[k];a=r.areal_traps
        if steady:
            f=rate.occupancy;jh=jl=rate.generation();df=0.
        else:
            src=rate.sources(f);jh=src['Sp'];jl=src['Sn'];df=src['df']
        block=np.r_[a*jl*dl,a*jh*dh,Q*a*(r.z_empty-f)*dt, np.eye(1,self.M,k).ravel()*df]
        # Each oriented branch is a discrete cumulative shape difference.
        # div J = -Q*a*(jh*(dh-dt)+jl*(dt-dl)), independently of sign/order.
        current=-Q*a*np.r_[0.,np.cumsum(jh*(wh-wt)+jl*(wt-wl))]
        drho=-Q*a*df*dt
        return block,current,drho,f
    def evaluate(self,state,steady=False):
        _,_,_,fs=self.unpack(state);b=np.zeros(3*self.N+self.M);j=np.zeros(self.N+1);dr=np.zeros(self.N);occ=[]
        for k,r in enumerate(self.reactions):
            rate=r.kernel(self.ports(state,k),r)
            bb,jj,dd,f=self._assemble(k,rate,fs[k],steady);b+=bb;j+=jj;dr+=dd;occ.append(f)
        return dict(Sn=b[:self.N],Sp=b[self.N:2*self.N],rho=b[2*self.N:3*self.N],df=b[3*self.N:],Jtransfer=j,drho=dr,f=np.array(occ),vector=b)
    def pattern(self,steady=False):
        p=lil_matrix((3*self.N+self.M,3*self.N+self.M),dtype=int)
        for k,(h,t,l) in enumerate(self.W):
            cols=np.r_[np.flatnonzero(h+t+l),self.N+np.flatnonzero(h),2*self.N+np.flatnonzero(l)]
            if not steady: cols=np.r_[cols,3*self.N+k]
            rows=np.r_[np.flatnonzero(l),self.N+np.flatnonzero(h),2*self.N+np.flatnonzero(t),3*self.N+k]
            for row in rows:p[row,cols]=1
        return p.tocsr()
    def jacobian(self,state,step=2e-5,steady=False):
        """5-point real differences on only five physical ports; analytic transient f.
        Sparse chain rule adds all endpoint couplings. No coloring assumptions.
        step is absolute energy step (eV/V); compare step and step/2 before use.
        This differentiates reaction sources/charge, not the entire device residual.
        """
        if not np.isfinite(step) or step<=0: raise ValueError('step')
        _,_,_,fs=self.unpack(state);out=lil_matrix((len(state),len(state)))
        for k,r in enumerate(self.reactions):
            ports=self.ports(state,k);h,t,l=self.W[k]
            maps=[(0,h),(0,t),(0,l),(self.N,h),(2*self.N,l)]
            for v,(offset,w) in enumerate(maps):
                vals=[]
                for m in [-2,-1,1,2]:
                    pp=ports.copy();pp[v]+=m*step
                    vals.append(self._assemble(k,r.kernel(pp,r),fs[k],steady)[0])
                d=(8*(vals[2]-vals[1])-(vals[3]-vals[0]))/(12*step)
                if not steady: d[2*self.N:3*self.N]=0. # rho has no port dependence
                for col in np.flatnonzero(w):
                    for row in np.flatnonzero(d):out[row,offset+col]+=d[row]*w[col]
            if not steady:
                rate=r.kernel(ports,r);aH,bH,aL,bL=np.exp([rate.laH,rate.lbH,rate.laL,rate.lbL]);dh,dt,dl=self.D[k];a=r.areal_traps
                d=np.r_[a*(aL+bL)*dl,-a*(aH+bH)*dh,-Q*a*dt,np.eye(1,self.M,k).ravel()*(-aH-bH-aL-bL)]
                for row in np.flatnonzero(d):out[row,3*self.N+k]+=d[row]
        return out.tocsr()
    def condensed_steady(self,fields,step=2e-5):
        """Eliminate f unknowns before stationary device assembly (no zero f rows)."""
        state=np.r_[fields,np.full(self.M,.5)]
        a=self.evaluate(state,steady=True)
        a["vector"]=a["vector"][:3*self.N].copy()
        return a,self.jacobian(state,step,True)[:3*self.N,:3*self.N]
    def continuity_defect(self,result):
        return Q*(result['Sp']-result['Sn'])+result['drho']+np.diff(result['Jtransfer'])/self.mesh.volumes
    def device_terms(self,state,steady=False):
        """Physical source API for a NEW FV device residual: Rn=-Sn,Rp=-Sp.
        Add rho to Poisson exactly once; integrate df for transient reactions.
        Existing local Nt closure must be replaced for this same trap family.
        """
        a=self.evaluate(state,steady)
        return dict(Rn=-a['Sn'],Rp=-a['Sp'],rho_trap=a['rho'],occupancy_rhs=a['df'],J_internal=a['Jtransfer'])

def boltzmann_ports(phi,n,p,Ec0,Ev0,Nc,Nv,temperature=300.):
    """Only for a declared compatible nondegenerate reservoir closure.
    Spectral wide-DOS models need DOS inversion, not this convenience mapping.
    """
    from rates import KB
    if not all(np.all(np.isfinite(v)) for v in (phi,n,p,Ec0,Ev0,Nc,Nv,temperature)): raise ValueError('finite carrier statistics')
    if np.any(np.asarray(n)<=0) or np.any(np.asarray(p)<=0) or Nc<=0 or Nv<=0 or temperature<=0: raise ValueError('positive carrier statistics')
    return np.asarray(Ev0)-phi-KB*temperature*np.log(np.asarray(p)/Nv),np.asarray(Ec0)-phi+KB*temperature*np.log(np.asarray(n)/Nc)
