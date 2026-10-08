import math, unittest
import numpy as np
from scipy.special import expit
from rates import *
ARGS=dict(t0=.001,distance=2.,xi=1.,lam=.3,temperature=300.) # diagnostic only; nm

def rates(muH=-.1,muL=.1,phiH=0.,phit=0.,phiL=0.):
    Et=-phit
    H=log_exchange(electron_energy(np.array([-.5,-.4,-.3]),phiH),[.2,.6,.2],muH,Et,**ARGS)
    L=log_exchange(electron_energy(np.array([.3,.4,.5]),phiL),[.3,.4,.3],muL,Et,**ARGS)
    return FourRates(*H,*L)

class Tests(unittest.TestCase):
    def test_kernel_detailed_balance(self):
        for d in [-1.,-.3,0,.2,1.]:self.assertAlmostEqual(float(log_kernel(d,**ARGS)-log_kernel(-d,**ARGS)),-d/(KB*300),places=11)
    def test_spectral_db(self):
        for mu in [-.6,0,.6]:
            a,b=log_exchange([-.3,.1,.5],[.1,3.,.5],mu,.07,**ARGS)
            self.assertAlmostEqual(a-b,(mu-.07)/(KB*300),places=11)
    def test_equilibrium_nonzero_field(self):
        r=rates(.03,.03,.15,-.1,-.3); f=expit((.03-.1)/(KB*300))
        self.assertAlmostEqual(r.occupancy,f,places=14)
        e=r.events(f)
        for a,b in [('hole_emission','hole_capture'),('electron_emission','electron_capture')]:
            self.assertLess(abs(e[a]-e[b])/max(e[a],e[b]),1e-12)
    def test_affinity_and_sign(self):
        for h,l in [(-.1,.1),(.1,-.1),(0,0)]:
            r=rates(h,l); self.assertAlmostEqual(r.affinity,(h-l)/(KB*300),places=11)
            self.assertLessEqual(-r.generation()*r.affinity,1e-20)
    def test_gauge(self):
        r=rates();shift=.73
        s=rates(-.1-shift,.1-shift,shift,shift,shift)
        np.testing.assert_allclose([r.laH,r.lbH,r.laL,r.lbL],[s.laH,s.lbH,s.laL,s.lbL],atol=1e-12)
    def test_distance(self):
        a=float(log_kernel(.1,**ARGS));kw=dict(ARGS,distance=3.)
        self.assertAlmostEqual(float(log_kernel(.1,**kw))-a,-2.)
    def test_field_sign_units(self):
        F=1.5e6;dx=-12e-7 # V/cm and cm; downhill electron displacement
        self.assertAlmostEqual(F*dx,-1.8)
        self.assertAlmostEqual(float(electron_energy(0.,-F*dx)),-1.8)
    def test_charge(self):
        for f in [0,.2,1]:
            s=rates().sources(f,1e15); scale=max(abs(Q*s['Sn']),abs(Q*s['Sp']),1e-300)
            self.assertLess(abs(Q*(s['Sp']-s['Sn'])+s['drho'])/scale,1e-14)
    def test_stationary_sources(self):
        r=rates();s=r.sources(r.occupancy)
        self.assertAlmostEqual(s['Sn']/r.generation(),1.,places=12)
        self.assertAlmostEqual(s['Sp']/r.generation(),1.,places=12)
    def test_probability(self):
        r=rates()
        for f in [0,.7,1]:
            for dt in [0,1e-12,1]:self.assertTrue(0<=r.advance(f,dt)<=1)
    def test_tiny_initial_occupation_zero_step(self):
        r=rates()
        self.assertEqual(r.advance(1e-20,0),1e-20)
        self.assertGreater(r.advance(1e-20,1e-50),0.)
    def test_large_finite_lograte_relaxation(self):
        r=FourRates(1000.,1000.,1000.,1000.)
        self.assertEqual(r.advance(0.,1.),.5)
    def test_large_lograte_small_affinity(self):
        r=FourRates(712.,712.,712.,712.+1e-10)
        self.assertTrue(math.isfinite(r.generation()))
        self.assertGreater(r.generation(),0.)
    def test_path_invalid_inputs(self):
        with self.assertRaises(ValueError):path_current_faces([0.,1.],0,1,2,1,1)
        with self.assertRaises(ValueError):path_current_faces([math.nan],0,1,2,1,1)
    def test_restore_bottleneck(self):
        r=from_positive(2,1e-20,1e-20,1000)
        self.assertAlmostEqual(r.generation(),2*1000/1002,places=12)
        self.assertLess(r.generation(),2.)
        s=from_positive(1e-15,1e-20,1e-20,1000)
        self.assertLess(abs(s.generation()),1.1e-15)
    def test_pair_enhancement_preserves_affinity(self):
        r=rates();s=FourRates(r.laH+3,r.lbH+3,r.laL,r.lbL)
        self.assertAlmostEqual(r.affinity,s.affinity,places=12)
    def test_single_emission_breaks_equilibrium(self):
        r=rates(0,0);s=FourRates(r.laH,r.lbH,r.laL,r.lbL+1)
        self.assertAlmostEqual(s.affinity,1.,places=12)
        self.assertGreater(s.generation(),0)
    def test_srh_reduction(self):
        cn,cp,n1,p1,n,p=1e-8,2e-8,1e10,2e10,3e9,4e9
        r=from_positive(cp*p1,cp*p,cn*n,cn*n1)
        G=cn*cp*(n1*p1-n*p)/(cn*(n+n1)+cp*(p+p1))
        self.assertAlmostEqual(r.generation()/G,1,places=12)
    def test_path_local_continuity(self):
        # nodes H=0,t=1,L=2, sources per area q*[jH,-(jH-jL),-jL]
        faces=np.array([-.5,.5,1.5,2.5]);jH=7.;jL=3.
        J=path_current_faces(faces,0,1,2,jH,jL)
        source=Q*np.array([jH,-jH+jL,-jL])
        np.testing.assert_allclose(np.diff(J),-source,rtol=1e-15,atol=0)
    def test_path_reverse_orientation(self):
        J=path_current_faces([-.5,.5,1.5,2.5],2,1,0,7.,3.)
        source=Q*np.array([-3.,-7.+3.,7.])
        np.testing.assert_allclose(np.diff(J),-source,rtol=1e-15,atol=0)
    def test_equilibrium_bias_not_field(self):
        r=rates(0,0,1.,.2,-.7)
        self.assertLess(abs(r.affinity),1e-10)
    def test_input_validation(self):
        with self.assertRaises(ValueError):log_kernel(0,**dict(ARGS,lam=-1))
        with self.assertRaises(ValueError):log_exchange([1],[0],0,0,**ARGS)
        with self.assertRaises(ValueError):rates().advance(-.1,1)
    def test_small_affinity_stability(self):
        r=FourRates(0,0,0,1e-12)
        self.assertAlmostEqual(r.generation()/2.5e-13,1,places=11)

if __name__=='__main__':unittest.main(verbosity=2)
