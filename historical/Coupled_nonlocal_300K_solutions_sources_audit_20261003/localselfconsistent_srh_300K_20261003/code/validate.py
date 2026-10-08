from integrate import *
from local_srh import LocalSRH
import unittest
class IntegrationTests(unittest.TestCase):
 def test_reference_equivalence(self):
  d=make(81);ref=DensityDevice(d.p,81,trap_e_depth=.65)
  for V in [0,-1,-5,-15]:
   z=np.load(ROOT/'data'/f'baseline_N81_V{V:g}.npz')['z'];a=d.evaluate(z,V);b=ref.evaluate(z,V)
   for k in ['f','fb','Rn','Rp','pair','Jn','Jp','bim']:np.testing.assert_allclose(a[k],b[k],rtol=2e-13,atol=1e-28)
   np.testing.assert_allclose(d.jacobian(z,V,0).toarray(),ref.jacobian(z,V,0).toarray(),rtol=2e-12,atol=2e-11)
 def test_scalar_four_rate_adapter(self):
  d=make(81);V=-1;z=np.load(ROOT/'data'/'baseline_N81_V-1.npz')['z'];a=d.local(z,V)
  for i in [10,30,40,60,70]:
   m=LocalSRH(d.p.capture,d.p.capture,d.ni,d.ni,d.p.Nt)
   self.assertAlmostEqual(m.occupancy(a['n'][i],a['p'][i]),a['f'][i],places=14)
   np.testing.assert_allclose(m.U(a['n'][i],a['p'][i]),a['U'][i],rtol=2e-13)
 def test_equilibrium_each_channel(self):
  for pf in [False,True]:
   d=make(161,pf);z=np.load(ROOT/'data'/f'{"singlePF" if pf else "baseline"}_N161_V0.npz')['z'];a=d.local(z,0)
   for fwd,bwd in [(a['cn']*a['n']*a['fb'],a['en']*a['f']),(a['cp']*a['p']*a['f'],a['ep']*a['fb'])]:np.testing.assert_allclose(fwd,bwd,rtol=3e-13)
   self.assertEqual(d.physical(z,0)['J_Acm2'],0.)
 def test_jacobian_direction_and_trap_charge(self):
  rng=np.random.default_rng(13)
  for pf in [False,True]:
   d=make(81,pf);V=-1;z=np.load(ROOT/'data'/f'{"singlePF" if pf else "baseline"}_N81_V-1.npz')['z'];v=rng.normal(size=z.shape);v[[0,-1]]=0;h=1e-5
   J=d.jacobian(z,V,0);a=J@v.ravel();b=(d.residual(z+h*v,V)-d.residual(z-h*v,V))/(2*h)
   self.assertLess(np.linalg.norm(a-b)/np.linalg.norm(a),2e-8)
   # Includes occupation feedback through n,p and field enhancement, not frozen rho_t.
   cs=d.local(z+1e-26j*v,V)['f'].imag/1e-26
   fd=(d.local(z+h*v,V)['f']-d.local(z-h*v,V)['f'])/(2*h)
   np.testing.assert_allclose(cs,fd,rtol=1e-5,atol=3e-10)
 def test_archive_and_grid(self):
  old=json.loads((ROOT/'reference_archived_points.json').read_text())
  lookup={r['V']:r for r in old['rows']}
  new=json.loads((ROOT/'data/baseline_N81.json').read_text())
  for r in new['rows']:np.testing.assert_allclose(r['J_Acm2'],lookup[r['V']]['J_Acm2'],rtol=2e-8,atol=1e-27)
  for mode in ['baseline','singlePF']:
   a=json.loads((ROOT/'data'/f'{mode}_N161.json').read_text())['rows'];b=json.loads((ROOT/'data'/f'{mode}_N321.json').read_text())['rows']
   for r,s in zip(a[1:],b[1:]):self.assertLess(abs(r['J_Acm2']/s['J_Acm2']-1),.002)
 def test_ledgers_all_points(self):
  for p in (ROOT/'data').glob('*N*.json'):
   for r in json.loads(p.read_text())['rows']:
    self.assertTrue(r['gate_passed']);self.assertLess(r['charge_relative_error'],1e-8)
    self.assertEqual(r['nonlocal_max_Acm2'],0)
if __name__=='__main__':unittest.main(verbosity=2)
