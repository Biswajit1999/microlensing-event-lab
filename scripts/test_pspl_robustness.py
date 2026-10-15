import unittest
import numpy as np
from fit_ogle_event import flux_model
from pspl_robustness import (
    analytic_flux_jacobian, fisher_covariance,
    observing_condition_diagnostics, parametric_bootstrap,
)


class TestRobustness(unittest.TestCase):
    def setUp(self):
        self.t=np.linspace(2459900.,2460100.,101)
        self.p=(2460001.5,24.,0.31,0.9,0.1)
        self.s=np.full(len(self.t),0.02)

    def test_analytic_jacobian_matches_central_difference(self):
        j=analytic_flux_jacobian(self.t,self.p)
        for k,h in enumerate((1e-4,1e-4,1e-6,1e-6,1e-6)):
            a=list(self.p); b=list(self.p)
            a[k]+=h;b[k]-=h
            finite=(flux_model(self.t,*a)-flux_model(self.t,*b))/(2*h)
            np.testing.assert_allclose(j[:,k],finite,rtol=1e-5,atol=2e-8)

    def test_fisher_is_positive_and_full_rank(self):
        f=fisher_covariance(self.t,self.s,self.p)
        self.assertTrue(f['identifiable'])
        self.assertEqual(f['rank'],5)
        self.assertTrue(all(v>0 for v in f['formal_sigma'].values()))
        eig=np.linalg.eigvalsh(f['covariance'])
        self.assertTrue(np.all(eig>0))

    def test_fisher_uncertainty_scales_with_error(self):
        a=fisher_covariance(self.t,self.s,self.p)
        b=fisher_covariance(self.t,2*self.s,self.p)
        for name in a['formal_sigma']:
            self.assertAlmostEqual(b['formal_sigma'][name]/a['formal_sigma'][name],2.,places=9)

    def test_fisher_rejects_nonpositive_errors(self):
        with self.assertRaises(ValueError):
            fisher_covariance(self.t,-self.s,self.p)

    def test_quality_associations(self):
        f=flux_model(self.t,*self.p)
        seeing=np.linspace(2,8,len(self.t))
        sky=np.linspace(200,600,len(self.t))
        flux=f+np.linspace(-0.03,0.03,len(self.t))
        r=observing_condition_diagnostics(flux,f,self.s,seeing,sky)
        self.assertEqual(r['n'],len(self.t))
        self.assertGreater(r['high_seeing_count'],0)
        self.assertGreaterEqual(r['high_seeing_mean_squared_z'],0)

    def test_quality_rejects_bad_data(self):
        f=flux_model(self.t,*self.p)
        with self.assertRaises(ValueError):
            observing_condition_diagnostics(f,f,self.s,np.ones(len(f)),np.ones(len(f)))

    def test_bootstrap_determinism_and_recovery(self):
        a=parametric_bootstrap(self.t,self.s,self.p,n_draws=12,seed=42)
        b=parametric_bootstrap(self.t,self.s,self.p,n_draws=12,seed=42)
        self.assertEqual(a,b)
        self.assertEqual(a['draws_successful'],12)
        for name,truth in zip(('t0_hjd','tE_days','u0'),self.p[:3]):
            low,high=a['central_68_percentile_intervals'][name]
            self.assertLess(low,high)
            self.assertLess(low,truth)
            self.assertGreater(high,truth)

    def test_bootstrap_rejects_too_few(self):
        with self.assertRaises(ValueError):
            parametric_bootstrap(self.t,self.s,self.p,n_draws=3)


if __name__=='__main__':
    unittest.main()
