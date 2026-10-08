import unittest, math, random
from local_srh import LocalSRH,Q_E,DA_STATES,DA_EDGES,da_equilibrium_edges,da_sources
class InterfaceTests(unittest.TestCase):
    def setUp(self): self.m=LocalSRH(1e-8,1e-8,1e10,1e10,Nt=1e12)
    def test_equilibrium_each_channel(self):
        for n in (1e6,1e10,1e14):
            p=1e20/n; f=self.m.occupancy(n,p); e=self.m.events(n,p,f)
            self.assertAlmostEqual(e['electron_capture'],e['electron_emission'],delta=1e-9)
            self.assertAlmostEqual(e['hole_capture'],e['hole_emission'],delta=1e-9)
            self.assertEqual(self.m.U(n,p),0.)
    def test_positive_invariant_occupancy(self):
        rng=random.Random(310)
        for _ in range(1000):
            n,p=[10**rng.uniform(0,18) for _ in range(2)]; f=rng.random()
            for dt in (0,1e-12,1e-3,1e5): self.assertTrue(0<=self.m.advance(n,p,f,dt)<=1)
            self.assertTrue(all(x>=0 for x in self.m.events(n,p,f).values()))
    def test_transient_charge_conservation(self):
        for f in (0,.1,.5,1):
            s=self.m.sources(2e12,3e8,f)
            scale=Q_E*self.m.Nt*(self.m.cn*2e12+self.m.cp*3e8+self.m.en+self.m.ep)
            self.assertLess(abs(Q_E*(s['Sp']-s['Sn'])+s['drho']),scale*1e-14)
    def test_steady_sources_equal_negative_U(self):
        for n,p in ((1e4,1e4),(2e12,3e12),(1e16,1e-2)):
            s=self.m.steady_sources(n,p); U=self.m.U(n,p)
            for k in ('Sn','Sp'): self.assertAlmostEqual(s[k],-U,delta=max(abs(U)*1e-7,1))
    def test_trap_charge_reference(self):
        a=self.m.sources(1,1,.3)['rho']; d=LocalSRH(1e-8,1e-8,1e10,1e10,1e12,1).sources(1,1,.3)['rho']
        self.assertAlmostEqual((d-a)/(Q_E*1e12),1.)
    def test_depletion_capture_counterexample(self):
        n,p=1e16,1e-2
        self.assertEqual(n*p/1e20,1e-6)
        eps=(self.m.cn*n+self.m.cp*p)/(self.m.en+self.m.ep)
        self.assertGreater(eps,4.9e5)
        full=-self.m.U(n,p); limit=self.m.Nt*self.m.en*self.m.ep/(self.m.en+self.m.ep)
        self.assertLess(full/limit,2.01e-6)
    def test_strong_depletion_serial_limit(self):
        for ratio in (1e-6,1,1e6):
            m=LocalSRH(1e-8*ratio,1e-8,1e10,1e10)
            G=-m.U(0,0); serial=m.en*m.ep/(m.en+m.ep)
            self.assertAlmostEqual(G/serial,1)
            self.assertLessEqual(G,min(m.en,m.ep)*(1+1e-15))
    def test_single_branch_saturates_both_branches_scale(self):
        g=1e6; base=-self.m.U(0,0)
        self.assertAlmostEqual(-self.m.enhance(g,1).U(0,0)/base,2*g/(g+1))
        self.assertAlmostEqual(-self.m.enhance(g,g).U(0,0)/base,g)
    def test_consistent_enhancement_equilibrium_and_recombination(self):
        e=self.m.enhance(50,30)
        self.assertEqual(e.U(1e12,1e8),0)
        self.assertGreater(e.U(1e12,1e12),self.m.U(1e12,1e12))
    def test_emission_only_breaks_cycle_balance(self):
        n=p=1e10; cn=self.m.cn; cp=self.m.cp; en=10*self.m.en; ep=self.m.ep
        f=(cn*n+ep)/(cn*n+ep+en+cp*p)
        Sn=en*f-cn*n*(1-f)
        self.assertGreater(Sn,0)
        self.assertAlmostEqual((cn*n*cp*p)/(en*ep),.1)
    def test_da_each_event_charge_and_particle_balance(self):
        for i,j,dn,dp,_ in DA_EDGES:
            zi=1-sum(DA_STATES[i]); zj=1-sum(DA_STATES[j])
            self.assertEqual(-dn+dp+zj-zi,0)
            # Total electron change free conduction + valence + localized =0.
            self.assertEqual(dn-dp+sum(DA_STATES[j])-sum(DA_STATES[i]),0)
    def test_da_probability_charge_conservation(self):
        s=da_sources([.1,.2,.3,.4],[(2,3),(5,7),(11,13),(17,19),(23,29)])
        self.assertAlmostEqual(sum(s['dP']),0)
        self.assertAlmostEqual(-s['Sn']+s['Sp']+s['dz'],0)
    def test_da_equilibrium_each_edge(self):
        pi,r=da_equilibrium_edges([0,.2,.1,.15],mu=.03)
        for (i,j,*_), (kf,kr) in zip(DA_EDGES,r): self.assertAlmostEqual(pi[i]*kf,pi[j]*kr)
        s=da_sources(pi,r)
        self.assertLess(max(abs(v) for v in s['dP']),1e-14)
    def test_da_complete_cycle_and_reset_bottleneck(self):
        # Generation cycle 10->01->00->10, irreversible diagnostic limit.
        for reset in (1e-6,1,1e6):
            k,en,ep=1e3,1e3,reset; R=1/(1/k+1/en+1/ep)
            P=[R/k,R/en,R/ep,0]; s=da_sources(P,[(k,0),(en,0),(ep,0),(0,0),(0,0)])
            self.assertLessEqual(R,min(k,en,ep))
            self.assertAlmostEqual(s['Sn'],R); self.assertAlmostEqual(s['Sp'],R)
            self.assertLess(max(abs(v) for v in s['dP']),1e-10)
    def test_invalid_inputs_rejected(self):
        for n,p in ((-1,1),(1,float('nan'))):
            with self.assertRaises(ValueError): self.m.occupancy(n,p)
if __name__=='__main__': unittest.main(verbosity=2)
