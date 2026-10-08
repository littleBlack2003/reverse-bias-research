import unittest
import numpy as np
from scipy.special import expit
from adapter import *
from rates import from_positive,KB

def mesh(n=31):
    # Graded mesh, not a uniform-grid proxy
    return Mesh(np.linspace(0,1,n+1)**1.3*1e-5)
def spectral():
    kw=dict(t0=1e-5,xi=1e-5,lam=.3)
    return SpectralRates(np.array([-.15,-.05]),np.array([.4,.6]),np.array([.08,.18]),np.array([.7,.3]),.015,kw,kw)
def setup(n=31,kernel=None,reverse=False):
    m=mesh(n);a,b=(2.31e-6,7.74e-6) if not reverse else (7.74e-6,2.31e-6)
    rs=[Reaction(a,4.62e-6,b,3e7,kernel or spectral())]
    d=Adapter(m,rs);phi=.1*(m.nodes/1e-5)**2
    return d,np.r_[phi,np.full(n,.02),np.full(n,-.03),.37]

class Tests(unittest.TestCase):
    def test_shape_moments(self):
        m=mesh()
        for x in np.linspace(m.nodes[0],m.nodes[-1],100):
            w=m.shape(x);self.assertAlmostEqual(sum(w),1);self.assertAlmostEqual(w@m.nodes,x,places=20);self.assertTrue(np.all(w>=0))
    def test_domain_reject(self):
        with self.assertRaises(ValueError):mesh().shape(0)
        with self.assertRaises(ValueError):Mesh([0,1,1])
    def test_signed_current_all_orders(self):
        import itertools
        for h,t,l in itertools.permutations([2e-6,4.62e-6,8e-6]):
            for f in [0,.37,1]:
                m=mesh();d=Adapter(m,[Reaction(h,t,l,3e7,FixedRates(from_positive(2,7,3,11)))])
                s=np.r_[np.zeros(3*m.size),f];a=d.evaluate(s)
                scale=max(abs(Q*a['Sp']).max(),abs(Q*a['Sn']).max(),abs(a['drho']).max())
                self.assertLess(max(abs(d.continuity_defect(a)))/scale,2e-14)
                self.assertLess(abs(a['Jtransfer'][-1]),1e-25);self.assertEqual(a['Jtransfer'][0],0)
    def test_integrated_number_charge(self):
        d,s=setup();a=d.evaluate(s);r=d.reactions[0];rate=r.kernel(d.ports(s,0),r);q=rate.sources(s[-1]);v=d.mesh.volumes
        np.testing.assert_allclose([a['Sn']@v,a['Sp']@v,a['rho']@v,a['drho']@v],[r.areal_traps*q['Sn'],r.areal_traps*q['Sp'],Q*r.areal_traps*(r.z_empty-s[-1]),r.areal_traps*q['drho']],rtol=2e-14)
    def test_zero_distance_local_srh(self):
        from local_srh_control import LocalSRH
        m=mesh();i=14;Nt=4e12;srh=LocalSRH(2e-9,3e-9,2e8,3e8,Nt);n=7e8;p=9e8;f=.4
        rates=from_positive(srh.ep,srh.cp*p,srh.cn*n,srh.en)
        d=Adapter(m,[Reaction(*([m.nodes[i]]*3),Nt*m.volumes[i],FixedRates(rates),srh.z_empty)])
        a=d.evaluate(np.r_[np.zeros(3*m.size),f]);b=srh.sources(n,p,f)
        for key in ['Sn','Sp','rho','drho']:np.testing.assert_allclose(a[key][i],b[key],rtol=1e-14)
        np.testing.assert_array_equal(a['Jtransfer'],0)
        ss=d.evaluate(np.r_[np.zeros(3*m.size),f],steady=True)
        np.testing.assert_allclose(ss['Sn'][i],-srh.U(n,p),rtol=1e-14)
    def test_equilibrium_common_mu(self):
        d,s=setup();mu=.021;s[d.N:3*d.N]=mu;et=.015-d.ports(s,0)[1];s[-1]=expit((mu-et)/(KB*300))
        a=d.evaluate(s);r=d.reactions[0];rr=r.kernel(d.ports(s,0),r);gross=np.exp(rr.logD)*r.areal_traps/min(d.mesh.volumes)
        self.assertLess(max(abs(a['Sn']).max(),abs(a['Sp']).max())/gross,1e-14)
        self.assertLess(abs(a['df'][0])/np.exp(rr.logD),1e-14)
    def test_gauge(self):
        d,s=setup();t=s.copy();t[:d.N]+=.73;t[d.N:3*d.N]-=.73
        for k in ['vector','Jtransfer']:np.testing.assert_allclose(d.evaluate(s)[k],d.evaluate(t)[k],rtol=1e-12,atol=1e-20)
    def test_steady_equal_integrals(self):
        d,s=setup(reverse=True);a=d.evaluate(s,steady=True)
        np.testing.assert_allclose(a['Sn']@d.mesh.volumes,a['Sp']@d.mesh.volumes,rtol=2e-15)
        self.assertEqual(a['df'][0],0)
        scale=max(abs(Q*a['Sn']).max(),abs(Q*a['Sp']).max());self.assertLess(max(abs(d.continuity_defect(a)))/scale,1e-14)
    def test_transient_integrated_exact(self):
        d,s=setup(kernel=FixedRates(from_positive(2,7,3,11)));r=d.reactions[0];rate=r.kernel(None,r);f0=s[-1];dt=.03;f1=rate.advance(f0,dt);D=np.exp(rate.logD);fi=rate.occupancy*dt+(f0-rate.occupancy)*(-np.expm1(-D*dt))/D
        # Time integrals are analytic; resulting cells obey charge + face transport.
        aH,bH,aL,bL=np.exp([rate.laH,rate.lbH,rate.laL,rate.lbL]);h=aH*dt-(aH+bH)*fi;l=(aL+bL)*fi-aL*dt
        wh,wt,wl=d.W[0];a=r.areal_traps
        dq=Q*a*(h*wh-l*wl-(f1-f0)*wt);integralJ=-Q*a*np.r_[0,np.cumsum(h*(wh-wt)+l*(wt-wl))]
        self.assertLess(max(abs(dq+np.diff(integralJ)))/(Q*a),1e-15)
    def test_jacobian_transient(self):self._jac(False)
    def test_jacobian_steady(self):self._jac(True)
    def _jac(self,steady):
        d,s=setup(n=13);j=d.jacobian(s,steady=steady).toarray();j2=d.jacobian(s,step=1e-5,steady=steady).toarray();sc=np.maximum(np.max(abs(j),axis=1),1e-30)
        self.assertLess(np.max(abs(j-j2)/sc[:,None]),3e-9)
        rng=np.random.default_rng(21);v=rng.normal(size=len(s));eps=1e-6
        fd=(d.evaluate(s+eps*v,steady)['vector']-d.evaluate(s-eps*v,steady)['vector'])/(2*eps)
        self.assertLess(max(abs(fd-j@v)/sc),2e-7)
        self.assertFalse(np.any((abs(j)>0)&(d.pattern(steady).toarray()==0)))
        self.assertGreater(d.pattern().nnz,0)
    def test_multiple_reactions(self):
        d,s=setup();r=d.reactions[0];r2=Reaction(8e-6,3e-6,2e-6,7e7,spectral(),0)
        both=Adapter(d.mesh,[r,r2]);a=both.evaluate(np.r_[s[:-1],.37,.8]);sep=Adapter(d.mesh,[r2]).evaluate(np.r_[s[:-1],.8]);one=d.evaluate(s)
        for k in ['Sn','Sp','rho','Jtransfer']:np.testing.assert_allclose(a[k],one[k]+sep[k])
        self.assertLess(max(abs(both.continuity_defect(a)))/max(abs(Q*a['Sn'])),1e-14)
    def test_grid_fixed_physical_distance(self):
        errors=[]
        for n in [21,41,81,161,321,641,1281]:
            d,s=setup(n);r=d.reactions[0];exact=np.array([.1*(r.xH/1e-5)**2,.1*(r.xt/1e-5)**2,.1*(r.xL/1e-5)**2,.02,-.03])
            reference=r.kernel(exact,r).generation();computed=d.evaluate(s,True)['Sn']@d.mesh.volumes/r.areal_traps
            self.assertLessEqual(max(abs(d.ports(s,0)-exact)),.1/1e-10*max(np.diff(d.mesh.nodes))**2/4+1e-15)
            errors.append(abs(computed/reference-1));self.assertAlmostEqual(r.xL-r.xH,5.43e-6,places=20)
            # Source first moments preserve physical endpoints on every grid.
            a=d.evaluate(s,True);self.assertAlmostEqual((a['Sn']*d.mesh.volumes)@d.mesh.nodes/(a['Sn']@d.mesh.volumes),r.xL,places=19)
        self.assertLess(errors[-1],1e-7)
        self.assertLess(errors[-1],errors[0]/100)
        print('grid_relative_errors',errors)
    def test_fixed_rates_zero_port_jacobian(self):
        d,s=setup(kernel=FixedRates(from_positive(1e8,2e8,3e8,4e8)))
        for steady in [False,True]:np.testing.assert_array_equal(d.jacobian(s,steady=steady).toarray()[:,:3*d.N],0)
    def test_statistics_reject_nonfinite(self):
        for n in [np.nan,np.inf]:
            with self.assertRaises(ValueError):boltzmann_ports(0,n,1,0,0,1,1)
    def test_condensed_steady(self):
        d,s=setup();a,j=d.condensed_steady(s[:3*d.N]);self.assertEqual(j.shape,(3*d.N,3*d.N))
        self.assertEqual(a["vector"].shape,(3*d.N,))
        np.testing.assert_allclose(a["vector"],d.evaluate(s,True)["vector"][:3*d.N])
    def test_device_signs(self):
        d,s=setup();a=d.evaluate(s);b=d.device_terms(s)
        np.testing.assert_array_equal(b['Rn'],-a['Sn']);np.testing.assert_array_equal(b['Rp'],-a['Sp'])
    def test_statistics_gauge(self):
        phi=np.array([.1,.3]);args=(np.array([2e8,3e8]),np.array([4e8,5e8]),.2,-.2,1e20,1e20)
        h,l=boltzmann_ports(phi,*args);hh,ll=boltzmann_ports(phi+.7,*args)
        np.testing.assert_allclose(hh,h-.7);np.testing.assert_allclose(ll,l-.7)
if __name__=='__main__':unittest.main(verbosity=2)
