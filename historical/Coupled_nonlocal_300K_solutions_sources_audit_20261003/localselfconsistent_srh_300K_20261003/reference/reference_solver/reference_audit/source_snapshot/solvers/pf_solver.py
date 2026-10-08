"""PF-inspired finite-distance trap escape coupled to 1D drift diffusion.

This is a hypothesis model, not a calibrated PM6:L8-BO microscopic model.
Each capture/emission pair obeys endpoint-energy detailed balance. Electron
and hole source locations are distinct; nonlocal current is explicitly counted.
Uses the shared base_dd_solver in this package.
"""
from pathlib import Path
import sys
from .base_dd_solver import Device,Params
from dataclasses import dataclass
import numpy as np

@dataclass
class PFParams(Params):
    pf_strength: float = 1.
    escape_nm: float = 5.
    nonlocal_fraction: float = .5
    barrier_smoothing_eV: float = .002

class PFDevice(Device):
    def __init__(self,p=None,nodes=161):
        super().__init__(p or PFParams(),nodes)
        self.jump=max(1,int(round(self.p.escape_nm/(self.p.d*1e7)* (nodes-1))))
        self.actual_escape_nm=self.jump*self.h*self.p.d*1e7
        self.pairs=[]
        for offset in [-self.jump,self.jump]:
            t=np.arange(1,nodes-1)
            end=t+offset
            ok=(end>=1)&(end<nodes-1)
            self.pairs.append((t[ok],end[ok]))
        # Sources at an endpoint depend on both endpoints and on trap potential.
        # Row stencil is <= +/- (2*jump+1); coloring must span its full width.
        radius=2*self.jump+1
        rr=[];cc=[]
        for i in range(nodes):
            for a in range(3):
                for k in range(max(0,i-radius),min(nodes,i+radius+1)):
                    for b in range(3):rr.append(3*i+a);cc.append(3*k+b)
        self.rr=np.array(rr);self.cc=np.array(cc)
        self.ncolors=3*(2*radius+1);self.colors=self.cc%self.ncolors

    @staticmethod
    def smoothmax(a,b,w):
        return .5*(a+b+np.sqrt((a-b)**2+w*w))

    def trap_rates(self,y):
        p=self.p;phi=y[:,0]*self.vt
        n=np.exp(y[:,1])*self.n0;holes=np.exp(y[:,2])*self.n0
        e=-np.gradient(phi,self.h)/p.d
        F=np.sqrt(e*e+1.)
        delta=p.pf_strength*np.sqrt(self.q*F/(np.pi*self.eps)) # eV
        depth=p.Eg/2
        c0=p.capture*(1-p.nonlocal_fraction)
        an0=c0*n;ap0=c0*holes;en0=c0*self.ni;ep0=en0
        # Sum occupancy transition rates, in s^-1 for each trap.
        fill=an0+ep0;empty=ap0+en0
        links=[]
        pref=p.capture*p.Nc*p.nonlocal_fraction/2
        w=p.barrier_smoothing_eV
        base_barrier=self.smoothmax(0.,depth-delta,w)
        for t,end in self.pairs:
            du=phi[end]-phi[t]
            for carrier,energy,concentration in [('n',depth-du,n[end]),('p',depth+du,holes[end])]:
                # A common transition state cannot lie below either endpoint.
                barrier=self.smoothmax(base_barrier[t],energy,w)
                emission=pref*np.exp(-barrier/self.vt)
                capture=(pref/p.Nc)*np.exp(-(barrier-energy)/self.vt)
                incoming=capture*concentration
                if carrier=='n':
                    fill[t]+=incoming;empty[t]+=emission
                else:
                    fill[t]+=emission;empty[t]+=incoming
                links.append((carrier,t,end,incoming,emission,energy,barrier))
        f=fill/(fill+empty)
        fbar=empty/(fill+empty)
        rn=p.Nt*(an0*fbar-en0*f)
        rp=p.Nt*(ap0*f-ep0*fbar)
        # Boundary nodes are carrier reservoirs, not trap reaction volumes.
        rn[[0,-1]]=0;rp[[0,-1]]=0
        for carrier,t,end,incoming,emission,energy,barrier in links:
            net=p.Nt*(incoming*fbar[t]-emission*f[t]) if carrier=='n' else p.Nt*(incoming*f[t]-emission*fbar[t])
            np.add.at(rn if carrier=='n' else rp,end,net)
        return rn,rp,f,e,delta,links,fbar

    def transport(self,y):
        n=np.exp(y[:,1]);p=np.exp(y[:,2]);dp=np.diff(y[:,0])
        bp=self.bern(dp);bm=self.bern(-dp)
        jn=(n[1:]*bp-n[:-1]*bm)/self.h
        jp=(p[:-1]*bp-p[1:]*bm)/self.h*(self.p.mu_p/self.p.mu_n)
        bim=self.p.beta*self.n0**2*(n*p-self.ni_s**2)
        return n,p,jn,jp,bim

    def residual(self,y,V,light):
        n,p,jn,jp,bim=self.transport(y)
        rn,rp,f,e,delta,links,fbar=self.trap_rates(y)
        g=light*self.p.Jgen/(self.q*self.p.d)
        netn=(rn+bim-g)/self.r0;netp=(rp+bim-g)/self.r0
        rho=p-n+self.p.Nt/self.n0*(.5-f)
        r=np.zeros_like(y)
        r[1:-1,0]=(y[:-2,0]-2*y[1:-1,0]+y[2:,0])/self.h**2+rho[1:-1]
        r[1:-1,1]=(jn[1:]-jn[:-1])/self.h-netn[1:-1]
        r[1:-1,2]=(jp[1:]-jp[:-1])/self.h+netp[1:-1]
        r[0]=[y[0,0],y[0,1]-self.bc[0,0],y[0,2]-self.bc[0,1]]
        r[-1]=[y[-1,0]-(self.vbi-V)/self.vt,y[-1,1]-self.bc[1,0],y[-1,2]-self.bc[1,1]]
        return r.ravel()

    def jacobian(self,y,V,light):
        from scipy.sparse import coo_matrix
        data=np.zeros(len(self.rr))
        for color in range(self.ncolors):
            z=y.astype(complex).ravel();z[color::self.ncolors]+=1e-26j
            derivative=self.residual(z.reshape(y.shape),V,light).imag/1e-26
            mask=self.colors==color;data[mask]=derivative[self.rr[mask]]
        return coo_matrix((data,(self.rr,self.cc)),shape=(y.size,y.size)).tocsc()

    def observe(self,y):
        n,p,jn,jp,bim=self.transport(y)
        rn,rp,f,e,delta,links,fbar=self.trap_rates(y)
        jdd=(jn+jp)*self.j0
        # Explicit current carried by the finite-distance trap transitions.
        jnl=np.r_[0.,-self.q*self.p.d*self.h*np.cumsum((rn-rp)[1:-1])]
        jt=jdd+jnl
        # No nonlocal path crosses either external contact. Its right value
        # must vanish from trap charge balance; terminal DD currents must agree.
        return dict(J_Acm2=float((jdd[0]+jdd[-1])/2),spread_Acm2=float(np.ptp(jt)),
            terminal_difference_Acm2=float(abs(jdd[0]-jdd[-1])),
            nonlocal_boundary_error_Acm2=float(abs(jnl[-1])),
            n=n*self.n0,p=p*self.n0,trap_occupancy=f,E=e,phi=y[:,0]*self.vt,
            Rn=rn,Rp=rp,Jdd=jdd,Jnonlocal=jnl,Jtotal=jt,
            barrier_lowering_eV=delta,Emax_Vcm=float(np.max(abs(e))))

if __name__=='__main__':
    for strength in [0.,1.]:
        dev=PFDevice(PFParams(pf_strength=strength),nodes=81)
        y,diag=dev.solve(0,0)
        print('equilibrium',strength,dev.observe(y)['J_Acm2'],flush=True)
        prev=0.
        for V in [-5.,-10.,-15.,-20.,-25.]:
            y,diag=dev.ramp(V,0,y,Vprev=prev,dv=.25);prev=V
            o=dev.observe(y)
            print(strength,V,o['J_Acm2'],o['spread_Acm2'],o['terminal_difference_Acm2'],o['nonlocal_boundary_error_Acm2'],flush=True)
