"""PM6:Y6 independent 1D isothermal Poisson/drift-diffusion reference.
Units cm, eV, V, s, A. Old conditional models remain untouched.
State psi=phi/Vt, u=EFn/Vt-(b_h+V)/Vt, v=-EFh/Vt+b_h/Vt.
Contact-referenced quasi-Fermi deviations resolve tiny dark currents near
high-conductivity contacts. Gaussian centers Ev=-phi, Ec=Eg-phi. No added trap channel.
Photogeneration is an external, irreversible phenomenological drive.
Generalized Scharfetter-Gummel uses edge secant inverse compressibility.
Contact scenarios are diagnostic assumptions, NEVER inferred work functions.
"""
from dataclasses import dataclass,asdict
import numpy as np
from scipy.special import exprel
from scipy.sparse import coo_matrix,diags
from scipy.sparse.linalg import spsolve
from fd_eos import GaussianDOS
KB=8.617333262145e-5; Q=1.602176634e-19; EPS0=8.8541878128e-14

@dataclass(frozen=True)
class Parameters:
    T:float=300.
    d_cm:float=110e-7
    N0_cm3:float=2.4e20
    sigma_n_eV:float=.060
    sigma_p_eV:float=.074
    mu_n_300:float=8.4e-4
    mu_p_300:float=1.3e-4
    Eg_eV:float=1.42
    epsilon_r:float=3.5 # scenario input; not identified by this paper
    beta_cm3s:float=8e-12
    G_cm3s:float=0.
    barrier_n_eV:float=0. # DOS-center ideal majority reservoir, not fitted
    barrier_p_eV:float=0.
    minority_contact:str='blocking' # or equilibrium_reservoir sensitivity
    allow_mobility_extrapolation:bool=False

def mobility(mu300,sigma,T):
    """Effective zero-field SCLC GDM; corrected sign; NOT zero-density EGDM."""
    return mu300*np.exp(-4/9*(sigma/KB)**2*(1/T**2-1/300**2))

def log_bernoulli(x):
    x=np.asarray(x);o=np.empty_like(x);p=x>50;n=x<-50;m=~(p|n)
    o[p]=np.log(x[p])-x[p]-np.log1p(-np.exp(-x[p]))
    o[n]=np.log(-x[n])-np.log1p(-np.exp(x[n]))
    o[m]=-np.log(exprel(x[m]))
    return o

def sg_flux(eta,psi,qf,logc,dlogc,mu,N0,vt,dx,carrier):
    """Thermodynamic secant SG; exact equilibrium, classical dilute limit.
    g=Delta eta/Delta ln(c)>0; midpoint derivative when secant is tiny.
    Stable affinity form avoids subtraction of almost equal gross currents.
    """
    if carrier not in ('n','p'): raise ValueError('carrier must be n or p')
    lc=logc(eta);de=np.diff(eta);dl=np.diff(lc)
    g=np.divide(de,dl,out=1/dlogc((eta[:-1]+eta[1:])/2),where=abs(de)>1e-7)
    if np.any(g<=0) or not np.all(np.isfinite(g)):raise FloatingPointError('nonpositive inverse compressibility')
    dp=np.diff(psi);a=np.diff(qf)/g
    arg=-dp/g if carrier=='n' else dp/g
    az=abs(a);le=np.full_like(a,-np.inf);nz=az>0
    le[nz]=np.maximum(a[nz],0)+np.log(-np.expm1(-az[nz]))
    la=np.log(Q*mu*N0*vt)+np.log(g)-np.log(dx)+lc[:-1]+log_bernoulli(arg)+le
    if np.any(la>650):raise FloatingPointError('unphysical trial flux overflow')
    return (1 if carrier=='n' else -1)*np.sign(a)*np.exp(la),g

class Device:
    def __init__(self,p=Parameters(),nodes=161,mesh_strength=5.):
        positive=[p.T,p.d_cm,p.N0_cm3,p.mu_n_300,p.mu_p_300,p.Eg_eV,p.epsilon_r,p.beta_cm3s]
        if not np.all(np.isfinite(positive)) or min(positive)<=0: raise ValueError('positive finite device scales required')
        if not np.isfinite(p.G_cm3s) or p.G_cm3s<0: raise ValueError('nonnegative finite generation required')
        if not np.all(np.isfinite([p.barrier_n_eV,p.barrier_p_eV])): raise ValueError('finite contact barriers required')
        if int(nodes)!=nodes or nodes<3: raise ValueError('at least three integer nodes required')
        if not np.isfinite(mesh_strength) or mesh_strength<0: raise ValueError('finite nonnegative mesh strength required')
        if not 223<=p.T<=328 and not p.allow_mobility_extrapolation:
            raise ValueError("SCLC mobility was measured approximately223–328K; explicit allow_mobility_extrapolation=True required")
        self.p=p;self.vt=KB*p.T;self.nodes=nodes
        t=np.linspace(-1,1,nodes)
        # Smooth symmetric edge clustering resolves ideal-contact accumulation.
        self.x=(t+1)/2 if mesh_strength==0 else (1+np.tanh(mesh_strength*t)/np.tanh(mesh_strength))/2
        self.dx=np.diff(self.x)*p.d_cm
        self.vol=(self.dx[:-1]+self.dx[1:])/2
        self.en=GaussianDOS(p.sigma_n_eV,p.T);self.ep=GaussianDOS(p.sigma_p_eV,p.T)
        self.mn=mobility(p.mu_n_300,p.sigma_n_eV,p.T);self.mp=mobility(p.mu_p_300,p.sigma_p_eV,p.T)
        self.vbi=p.Eg_eV-p.barrier_n_eV-p.barrier_p_eV
        self.eps=EPS0*p.epsilon_r
        rr=[];cc=[]
        for i in range(nodes):
            for a in range(3):
                for j in range(max(0,i-1),min(nodes,i+2)):
                    for b in range(3):rr.append(3*i+a);cc.append(3*j+b)
        self.rr=np.array(rr);self.cc=np.array(cc);self.colors=self.cc%9
    def evaluate(self,z,V):
        p=self.p;psi,u,v=z.T
        en=u+psi+(p.barrier_p_eV+V-p.Eg_eV)/self.vt;ep=v-psi-p.barrier_p_eV/self.vt
        ln=self.en.logc(en);lp=self.ep.logc(ep)
        n=p.N0_cm3*np.exp(ln);h=p.N0_cm3*np.exp(lp)
        jn,gn=sg_flux(en,psi,u,self.en.logc,self.en.dlogc,self.mn,p.N0_cm3,self.vt,self.dx,'n')
        jp,gp=sg_flux(ep,psi,v,self.ep.logc,self.ep.dlogc,self.mp,p.N0_cm3,self.vt,self.dx,'p')
        A=V/self.vt+u+v
        # Sign-preserving log form of beta*n*p*(1-exp(-A)).
        la=np.full_like(A,-np.inf);nz=abs(A)>0
        la[nz]=np.maximum(-A[nz],0)+np.log(-np.expm1(-abs(A[nz])))
        lr=np.log(p.beta_cm3s)+2*np.log(p.N0_cm3)+ln+lp+la
        if np.any(lr>650):raise FloatingPointError('unphysical trial reaction overflow')
        R=np.sign(A)*np.exp(lr)
        return dict(n=n,p=h,Jn=jn,Jp=jp,R=R,A=A,eta_n=en,eta_p=ep,g_n=gn,g_p=gp)
    def residual(self,z,V,light=0.):
        p=self.p;o=self.evaluate(z,V);r=np.zeros_like(z)
        grad=np.diff(z[:,0])/self.dx
        r[1:-1,0]=np.diff(grad)+Q/(self.eps*self.vt)*self.vol*(o['p']-o['n'])[1:-1]
        source=Q*self.vol*(o['R'][1:-1]-light*p.G_cm3s)
        r[1:-1,1]=np.diff(o['Jn'])-source
        r[1:-1,2]=np.diff(o['Jp'])+source
        r[0,0]=z[0,0];r[-1,0]=z[-1,0]-(self.vbi-V)/self.vt
        r[0,2]=z[0,2]
        r[-1,1]=z[-1,1]
        if p.minority_contact=='blocking':r[0,1]=o['Jn'][0];r[-1,2]=o['Jp'][-1]
        elif p.minority_contact=='equilibrium_reservoir':
            r[0,1]=z[0,1]+V/self.vt
            r[-1,2]=z[-1,2]+V/self.vt
        else:raise ValueError('unknown contact scenario')
        return r.ravel()
    def jacobian(self,z,V,light):
        data=np.zeros(len(self.rr));h=1e-5
        for color in range(9):
            zp=z.copy().ravel();zm=z.copy().ravel();zp[color::9]+=h;zm[color::9]-=h
            der=(self.residual(zp.reshape(z.shape),V,light)-self.residual(zm.reshape(z.shape),V,light))/(2*h)
            mask=self.colors==color;data[mask]=der[self.rr[mask]]
        return coo_matrix((data,(self.rr,self.cc)),shape=(z.size,z.size)).tocsc()
    def ledger(self,z,V,light):
        o=self.evaluate(z,V);j=o['Jn']+o['Jp'];G=Q*light*self.p.G_cm3s*sum(self.vol);R=Q*np.dot(self.vol,o['R'][1:-1])
        baln=o['Jn'][-1]-o['Jn'][0]-R+G;balp=o['Jp'][-1]-o['Jp'][0]+R-G
        cell_source=Q*self.vol*(o['R'][1:-1]-light*self.p.G_cm3s)
        local_n=np.diff(o['Jn'])-cell_source; local_p=np.diff(o['Jp'])+cell_source
        local_l1=float(max(np.sum(abs(local_n)),np.sum(abs(local_p))))
        local_linf=float(max(np.max(abs(local_n)),np.max(abs(local_p))))
        scale=max(abs(G),abs(R),max(abs(j)),1e-30);err=max(np.ptp(j),abs(baln),abs(balp))
        grad=np.diff(z[:,0])*self.vt/self.dx
        charge=Q*np.dot(self.vol,(o['p']-o['n'])[1:-1]);fd=self.eps*(grad[0]-grad[-1]);gross=Q*np.dot(self.vol,(o['p']+o['n'])[1:-1])
        td=float(np.dot(o['Jn'],np.diff(z[:,1])*self.vt)-np.dot(o['Jp'],np.diff(z[:,2])*self.vt))
        rd=float(Q*np.dot(self.vol,(o['R']*o['A']*self.vt)[1:-1]))
        ld=float(Q*light*self.p.G_cm3s*np.dot(self.vol,o['A'][1:-1]*self.vt))
        terminal=float(V*np.mean(j));energy=td+rd-ld
        # Closure-based precision indicator only: excludes mesh, constitutive,
        # contact and calibration uncertainty. Include cancellation roundoff.
        current_mean=float(np.mean(j))
        uncertainty=float(np.ptp(j)+max(err,local_l1)+64*np.finfo(float).eps*max(np.max(abs(o['Jn'])+abs(o['Jp'])),1e-300))
        current_tol=float(1e-3*abs(current_mean)+256*np.finfo(float).eps*max(abs(G),abs(R),np.max(abs(o['Jn'])+abs(o['Jp'])),1e-300))
        current_gate=bool(uncertainty<=current_tol)
        return dict(T=self.p.T,V=V,light=light,nodes=self.nodes,local_continuity_L1_Acm2=local_l1,local_continuity_Linf_Acm2=local_linf,current_accuracy_gate_passed=current_gate,current_gate_tolerance_Acm2=current_tol,J_Acm2=float(np.mean(j)),current_spread_Acm2=float(np.ptp(j)),current_spread_relative_to_terminal=float(np.ptp(j)/max(abs(current_mean),1e-300)),current_closure_uncertainty_Acm2=uncertainty,numerical_sign_resolved=bool(abs(current_mean)>uncertainty),continuity_abs_Acm2=float(err),continuity_relative=float(err/scale),charge_relative=float(abs(charge-fd)/max(gross,1e-300)),generation_Acm2=float(G),recombination_Acm2=float(R),n_mean_cm3=float(np.dot(self.vol,o['n'][1:-1])/sum(self.vol)),p_mean_cm3=float(np.dot(self.vol,o['p'][1:-1])/sum(self.vol)),max_density_fraction=float(max(o['n'].max(),o['p'].max())/self.p.N0_cm3),transport_dissipation_Wcm2=td,recombination_dissipation_Wcm2=rd,light_chemical_work_Wcm2=ld,terminal_Wcm2=terminal,energy_relative=float(abs(energy-terminal)/max(abs(td),abs(rd),abs(ld),abs(terminal),1e-30)),mu_n_cm2Vs=self.mn,mu_p_cm2Vs=self.mp,mobility_extrapolated=bool(self.p.T<223 or self.p.T>328),contact_status='unvalidated idealized scenario',G_status='conditional temperature-independent generation',gate_passed=bool(current_gate and err<1e-12+1e-7*scale and abs(charge-fd)<1e-7*max(gross,1e-30)))
    def initial(self,V=0.):
        return np.c_[self.x*(self.vbi-V)/self.vt,np.zeros(self.nodes),np.zeros(self.nodes)]

def solve(device,V,z,light=0.,maxiter=160):
    z=z.copy();history=[]
    for it in range(maxiter):
        r=device.residual(z,V,light);jac=device.jacobian(z,V,light)
        scale=np.maximum(np.asarray(abs(jac).max(axis=1).toarray()).ravel(),1e-250)
        norm=float(max(abs(r)/scale));audit=device.ledger(z,V,light)
        history.append(dict(it=it,norm=norm,continuity=audit['continuity_relative']))
        if norm<3e-10 and audit['gate_passed']:
            return z,dict(success=True,iterations=it,norm=norm,history=history,**audit)
        step=spsolve(diags(1/scale)@jac,-r/scale).reshape(z.shape)
        if not np.all(np.isfinite(step)):return z,dict(success=False,reason='nonfinite Newton step',history=history,**audit)
        alpha=min(1.,4/max(np.max(abs(step[:,1:])),1e-30),8/max(np.max(abs(step[:,0])),1e-30))
        for ls in range(32):
            trial=z+alpha*step
            try:nn=float(max(abs(device.residual(trial,V,light))/scale))
            except (ValueError,FloatingPointError):nn=np.inf
            if nn<norm or (norm<3e-10 and nn<3e-10):z=trial;break
            alpha*=.5
        else:return z,dict(success=False,reason='line search stalled',norm=norm,history=history,**audit)
    return z,dict(success=False,reason='iteration limit',norm=norm,history=history,**audit)

def advance(d,z,V0,V1,L0,L1,trace=None,depth=0):
    guess=z.copy();dv=(V1-V0)/d.vt
    guess[:,0]-=dv*d.x;guess[:,1]-=dv*(1-d.x);guess[:,2]-=dv*d.x
    zz,info=solve(d,V1,guess,L1)
    if trace is not None:trace.append({k:v for k,v in info.items() if k!='history'})
    if info['success']:return zz,info
    if depth>=12:raise RuntimeError(f'Continuation unresolved T={d.p.T}, V={V1}, L={L1}: {info["reason"]}')
    vm=(V0+V1)/2;lm=(L0+L1)/2
    mid,_=advance(d,z,V0,vm,L0,lm,trace,depth+1)
    return advance(d,mid,vm,V1,lm,L1,trace,depth+1)

def poisson_seed(d,V):
    """Fixed reservoir quasi-Fermi Poisson seed, not a solved illuminated device."""
    z=d.initial(V)
    for it in range(150):
        o=d.evaluate(z,V);r=d.residual(z,V,0).reshape(-1,3)[:,0]
        dn=o['n']*d.en.dlogc(o['eta_n']);dp=o['p']*d.ep.dlogc(o['eta_p'])
        mid=np.ones(d.nodes);mid[1:-1]=-(1/d.dx[:-1]+1/d.dx[1:])-Q/(d.eps*d.vt)*d.vol*(dn+dp)[1:-1]
        lower=np.r_[1/d.dx[:-1],0];upper=np.r_[0,1/d.dx[1:]]
        jac=diags([lower,mid,upper],[-1,0,1],shape=(d.nodes,d.nodes)).tocsc()
        step=spsolve(jac,-r)
        if max(abs(step))<1e-10:return z
        z[:,0]+=min(1,4/max(abs(step)))*step
    raise RuntimeError('Poisson light seed did not converge')
