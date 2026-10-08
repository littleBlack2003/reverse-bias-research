"""1D steady-state Poisson/Scharfetter-Gummel drift-diffusion prototype.
All physical inputs use cm, s, V, A. Not fitted to a specific OSC.
"""
from dataclasses import dataclass, asdict
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

@dataclass
class Params:
    d: float = 1e-5
    temperature: float = 300.
    eps_r: float = 3.5
    Nc: float = 1e19
    Eg: float = 1.3
    contact_barrier: float = .1
    mu_n: float = 1e-3
    mu_p: float = 1e-3
    beta: float = 1e-11
    capture: float = 1e-9
    Nt: float = 1e15
    Jgen: float = .025
    # Optional LOCAL, detailed-balance-preserving SRH multiplier.
    # Not a validated organic TAT model, not nonlocal tunneling.
    enhancement: float = 0.
    Fref: float = 5e5
    Bfield: float = 2e6

class Device:
    q=1.602176634e-19
    def __init__(self,p=None,nodes=121):
        self.p=p or Params(); p=self.p
        self.nodes=nodes; self.x=np.linspace(0,1,nodes); self.h=1/(nodes-1)
        self.vt=8.617333262e-5*p.temperature
        self.eps=8.8541878128e-14*p.eps_r
        self.n0=self.eps*self.vt/(self.q*p.d**2)
        self.j0=self.q*p.mu_n*self.vt*self.n0/p.d
        self.r0=self.j0/(self.q*p.d)
        self.ni=p.Nc*np.exp(-p.Eg/(2*self.vt)); self.ni_s=self.ni/self.n0
        maj=p.Nc*np.exp(-p.contact_barrier/self.vt)
        minor=self.ni**2/maj
        self.bc=np.log(np.array([[minor,maj],[maj,minor]])/self.n0)
        self.vbi=p.Eg-2*p.contact_barrier
        # Residual at node i depends only on the three adjacent nodes.
        rr=[];cc=[]
        for i in range(nodes):
            for a in range(3):
                for k in range(max(0,i-1),min(nodes,i+2)):
                    for b in range(3): rr.append(3*i+a);cc.append(3*k+b)
        self.rr=np.array(rr);self.cc=np.array(cc)
        self.colors=self.cc%9

    @staticmethod
    def bern(z):
        z=np.asarray(z); out=np.empty_like(z)
        sm=np.abs(z)<1e-4
        out[sm]=1-z[sm]/2+z[sm]**2/12-z[sm]**4/720
        out[~sm]=z[~sm]/np.expm1(z[~sm])
        return out

    def fields(self,y):
        phi,ln,lp=y.T
        n=np.exp(ln); pp=np.exp(lp)
        dp=np.diff(phi)
        bp=self.bern(dp);bm=self.bern(-dp)
        jn=(n[1:]*bp-n[:-1]*bm)/self.h
        jp=(pp[:-1]*bp-pp[1:]*bm)/self.h*(self.p.mu_p/self.p.mu_n)
        occ=(n+self.ni_s)/(n+pp+2*self.ni_s)
        e=-np.gradient(phi,self.h)*self.vt/self.p.d
        # Smooth analytic norm so complex-step Jacobians remain valid.
        fa=np.sqrt(e*e+1.)
        gamma=self.p.enhancement*(fa/self.p.Fref)**2*np.exp(-self.p.Bfield*(1/fa-1/self.p.Fref))
        srh=self.p.capture*self.p.Nt*self.n0*(n*pp-self.ni_s**2)/(n+pp+2*self.ni_s)*(1+gamma)
        bim=self.p.beta*self.n0**2*(n*pp-self.ni_s**2)
        return n,pp,occ,e,jn,jp,srh,bim,gamma

    def residual(self,y,V,light):
        n,p,f,e,jn,jp,srh,bim,gamma=self.fields(y)
        r=np.zeros_like(y)
        net=(srh+bim-light*self.p.Jgen/(self.q*self.p.d))/self.r0
        rho=p-n+self.p.Nt/self.n0*(.5-f)
        r[1:-1,0]=(y[:-2,0]-2*y[1:-1,0]+y[2:,0])/self.h**2+rho[1:-1]
        r[1:-1,1]=(jn[1:]-jn[:-1])/self.h-net[1:-1]
        r[1:-1,2]=(jp[1:]-jp[:-1])/self.h+net[1:-1]
        r[0]=[y[0,0],y[0,1]-self.bc[0,0],y[0,2]-self.bc[0,1]]
        r[-1]=[y[-1,0]-(self.vbi-V)/self.vt,y[-1,1]-self.bc[1,0],y[-1,2]-self.bc[1,1]]
        return r.ravel()

    def jacobian(self,y,V,light):
        data=np.zeros(len(self.rr))
        for color in range(9):
            z=y.astype(complex).ravel();z[color::9]+=1e-26j
            deriv=self.residual(z.reshape(y.shape),V,light).imag/1e-26
            mask=self.colors==color
            data[mask]=deriv[self.rr[mask]]
        return coo_matrix((data,(self.rr,self.cc)),shape=(y.size,y.size)).tocsc()

    def solve(self,V,light=1.,guess=None,maxiter=100):
        if guess is None:
            phi=self.x*self.vbi/self.vt
            guess=np.c_[phi,self.bc[0,0]+phi,self.bc[0,1]-phi]
        y=guess.copy()
        # Enforce new contact potentials without shifting the interior guess.
        y[0]=[0,*self.bc[0]];y[-1]=[(self.vbi-V)/self.vt,*self.bc[1]]
        for it in range(maxiter):
            r=self.residual(y,V,light)
            jac=self.jacobian(y,V,light)
            # Scale rows by Jacobian magnitude for conditioning and line search.
            scale=np.asarray(abs(jac).max(axis=1).toarray()).ravel()
            scale=np.maximum(scale,1e-20)
            norm=np.max(np.abs(r)/scale)
            if norm<2e-12:
                return y,dict(iterations=it,residual_scaled=float(norm))
            step=spsolve(jac,-r).reshape(y.shape)
            if not np.all(np.isfinite(step)): raise RuntimeError('nonfinite Newton step')
            alpha=min(1.,3./max(np.max(np.abs(step[:,1:])),1e-30),20./max(np.max(np.abs(step[:,0])),1e-30))
            for ls in range(25):
                trial=y+alpha*step
                rn=self.residual(trial,V,light)
                if np.all(np.isfinite(rn)) and np.max(np.abs(rn)/scale)<norm:
                    y=trial;break
                alpha*=.5
            else:
                if norm<2e-8: return y,dict(iterations=it,residual_scaled=float(norm))
                raise RuntimeError(f'line search failed V={V} light={light}, norm={norm}')
        raise RuntimeError(f'Newton did not converge at V={V} light={light}; norm={norm}')

    def ramp(self,V,light,y,Vprev=0.,Lprev=0.,dv=.05,dl=.05):
        steps=max(1,int(np.ceil(abs(V-Vprev)/dv)),int(np.ceil(abs(light-Lprev)/dl)))
        for k in range(1,steps+1):
            vk=Vprev+(V-Vprev)*k/steps;lk=Lprev+(light-Lprev)*k/steps
            y,diag=self.solve(vk,lk,y)
        return y,diag

    def observe(self,y):
        n,p,f,e,jn,jp,srh,bim,gamma=self.fields(y)
        jt=(jn+jp)*self.j0
        # Geometry: holes left, electrons right; Jx is negative under generation.
        return dict(J_Acm2=float(np.mean(jt)),spread_Acm2=float(np.ptp(jt)),
                    Emax_Vcm=float(np.max(np.abs(e))),n=n*self.n0,p=p*self.n0,
                    trap_occupancy=f,E=e,phi=y[:,0]*self.vt,
                    Jn_x=jn*self.j0,Jp_x=jp*self.j0,SRH=srh,bimolecular=bim,gamma=gamma)

if __name__=='__main__':
    dev=Device(nodes=81)
    y,diag=dev.solve(0,0)
    print('equilibrium',diag,{k:v for k,v in dev.observe(y).items() if np.isscalar(v)},flush=True)
    y,diag=dev.ramp(0,1,y)
    print('short circuit',diag,{k:v for k,v in dev.observe(y).items() if np.isscalar(v)},flush=True)
    prev=0
    for V in [.2,.4,.6,.8,1.]:
        y,diag=dev.ramp(V,1,y,Vprev=prev,Lprev=1);prev=V
        print(V,diag,{k:v for k,v in dev.observe(y).items() if np.isscalar(v)},flush=True)
