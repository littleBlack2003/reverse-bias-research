import unittest,numpy as np
from coupled import Coupled,ROOT
class CoupledTests(unittest.TestCase):
 def test_common_mu_each_edge(self):
  for n in [81,161,321]:
   for dist in [0,1]:
    d=Coupled(n,dist);z=np.load(ROOT/'data'/f'd{dist}_N{n}_V0.npz')['z'];o=d.evaluate(z,0);f=o['actual_f'];fb=o['actual_fb']
    # Gross channel balance at the common electrochemical potential.
    np.testing.assert_allclose(o['aH']*fb,o['bH']*f,rtol=2e-12,atol=1e-25)
    np.testing.assert_allclose(o['aL']*fb,o['bL']*f,rtol=2e-12,atol=1e-25)
    self.assertTrue(np.all(o['j']==0));self.assertEqual(d.physical(z,0)['J_Acm2'],0)
 def test_invalid_parameters(self):
  for kwargs in [dict(distance_nm=-1),dict(distance_nm=50),dict(t0=np.nan),dict(lam=np.inf),dict(xi_nm=0)]:
   with self.assertRaises(ValueError):Coupled(81,**kwargs)
 def test_fixed_trap_count_and_no_contact_source(self):
  for n in [81,161,321]:
   for dist in [0,1]:
    d=Coupled(n,dist);self.assertAlmostEqual(sum(d.a)/9.6e9,1,places=12)
    for W in [d.H,d.L,d.T]:self.assertEqual(W[:,[0,-1]].nnz,0)
 def test_all_saved_ledgers(self):
  for n in [81,161,321]:
   for dist in [0,1]:
    d=Coupled(n,dist)
    for v in [0,-1,-5,-15]:
     z=np.load(ROOT/'data'/f'd{dist}_N{n}_V{v}.npz')['z'];r,_=d.ledger(z,v)
     self.assertTrue(r['gate_passed']);self.assertLess(r['energy_relative_error'],1e-8);self.assertLess(r['charge_relative_error'],1e-8)
if __name__=='__main__':unittest.main(verbosity=2)
