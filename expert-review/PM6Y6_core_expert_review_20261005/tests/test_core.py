"""NEW lightweight closure tests; no full JV scan and no network dependency."""
from pathlib import Path
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import unittest
import numpy as np
import review
from device import mobility, sg_flux, KB
from fd_eos import GaussianDOS, direct_eos


class CoreTests(unittest.TestCase):
    def test_01_original_source_identity(self):
        self.assertEqual(review.check_core_identity(), 10)

    def test_02_mobility_reference_and_temperature(self):
        self.assertEqual(mobility(8.4e-4, .060, 300), 8.4e-4)
        np.testing.assert_allclose(mobility(8.4e-4, .060, 100), 4.041352977282537e-12, rtol=1e-13)
        np.testing.assert_allclose(mobility(1.3e-4, .074, 100), 2.8954602150509694e-17, rtol=1e-13)

    def test_03_fd_eos(self):
        for T in (100, 300):
            dos = GaussianDOS(.074, T)
            eta = np.array([-60., -15., -3., 0., 3., 15.])
            c, dc, g = dos.evaluate(eta)
            cq, dq = direct_eos(eta, dos.s)
            np.testing.assert_allclose(c, cq, rtol=2e-7, atol=1e-14)
            np.testing.assert_allclose(dc, dq, rtol=2e-6, atol=1e-14)
            np.testing.assert_allclose(c+np.exp(dos.logc(-eta)), 1., atol=2e-15)
            self.assertTrue(np.all(dc > 0) and np.all(g >= 1.))
            np.testing.assert_allclose(dos.eta_from_c(c), eta, atol=1e-9)
            self.assertAlmostEqual(float(dos.evaluate(-250.)[2]), 1.)
        with self.assertRaises(ValueError): GaussianDOS(.074, 90)

    def test_04_sg_zero_affinity(self):
        dos = GaussianDOS(.060, 300)
        psi = np.linspace(0, 1, 5); eta = psi-10; qf = np.zeros(5)
        for carrier in ('n', 'p'):
            j, g = sg_flux(eta, psi, qf, dos.logc, dos.dlogc, 8.4e-4, 2.4e20,
                           KB*300, np.full(4, 1e-6), carrier)
            np.testing.assert_array_equal(j, np.zeros(4))
            self.assertTrue(np.all(g > 0))

    def test_05_field_factor_no_density_multiplier(self):
        d, states, _ = review.load_case(300)
        z = states['sc'][2].copy(); o = d.evaluate(z, 0.)
        E = np.minimum(o['field_abs_Vcm'], d.field_cap)
        sq = (E*E+d.field_regularization**2)**.25-np.sqrt(d.field_regularization)
        np.testing.assert_allclose(o['mobility_factor_n'], np.exp(d.gamma_n*sq), rtol=1e-14)
        np.testing.assert_allclose(o['mobility_factor_p'], np.exp(d.gamma_p*sq), rtol=1e-14)
        z[:, 1] += .01; different_density = d.evaluate(z, 0.)
        np.testing.assert_array_equal(o['mobility_factor_n'], different_density['mobility_factor_n'])
        np.testing.assert_array_equal(o['mobility_factor_p'], different_density['mobility_factor_p'])

    def test_06_recombination_detailed_balance(self):
        d, states, _ = review.load_case(300)
        z = d.initial(0.)
        np.testing.assert_array_equal(d.evaluate(z, 0.)['R'], np.zeros(d.nodes))
        for sign in (-1., 1.):
            z[:, 1] = sign*.1
            o = d.evaluate(z, 0.)
            expected = d.p.beta_cm3s*o['n']*o['p']*(-np.expm1(-o['A']))
            np.testing.assert_allclose(o['R'], expected, rtol=3e-14)
            self.assertTrue(np.all(o['R']*o['A'] >= 0))

    def test_08_boundaries(self):
        d, states, _ = review.load_case(300)
        z = states['sc'][2]
        o = d.evaluate(z, 0.)
        np.testing.assert_allclose(z[0, 0], 0., atol=1e-12)
        np.testing.assert_allclose(z[-1, 0], d.vbi/d.vt, atol=1e-12)
        np.testing.assert_allclose(o['p'][0]/d.p.N0_cm3, .5, atol=1e-13)
        np.testing.assert_allclose(o['n'][-1]/d.p.N0_cm3, .5, atol=1e-13)
        self.assertLess(abs(o['Jn'][0]), 1e-12)
        self.assertLess(abs(o['Jp'][-1]), 1e-12)

    def test_09_zero_field_and_cap(self):
        d, states, _ = review.load_case(300)
        z = np.zeros_like(states['sc'][2])
        o = d.evaluate(z, 0.)
        np.testing.assert_array_equal(o['mobility_factor_n'], np.ones(d.nodes-1))
        np.testing.assert_array_equal(o['mobility_factor_p'], np.ones(d.nodes-1))
        z[:, 0] = d.x*d.p.d_cm/d.vt*(2*d.field_cap)
        o = d.evaluate(z, 0.)
        cap_factor = np.exp(d.gamma_p*((d.field_cap**2+1)**.25-1))
        np.testing.assert_allclose(o['mobility_factor_p'], cap_factor, rtol=1e-14)

    def test_07_saved_states_residual_and_no_shunt(self):
        report = review.smoke()
        self.assertTrue(report['passed']); self.assertEqual(len(report['states']), 14)


if __name__ == '__main__':
    unittest.main()
