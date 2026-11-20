import unittest
from finite_exposure import (magnification, pspl_flux, exposure_flux,
                             quadrature_convergence, blend_peak_degeneracy)


class ExposureTests(unittest.TestCase):
    def test_zero_width_exact(self):
        self.assertEqual(exposure_flux(1, 0, 0, 20, .2, 2, .3),
                         pspl_flux(1, 0, 20, .2, 2, .3))

    def test_peak_smeared_down(self):
        instantaneous = pspl_flux(0, 0, 1, .1, 1)
        self.assertLess(exposure_flux(0, .4, 0, 1, .1, 1), instantaneous)

    def test_convergence(self):
        _, errors = quadrature_convergence(0, .4, 0, 1, .1, 1)
        self.assertLess(errors[-1], errors[0])

    def test_blend_family_recovers_flux(self):
        source, blend = blend_peak_degeneracy(3, 1, .2)
        self.assertAlmostEqual(source + blend, 1)
        self.assertAlmostEqual(source * magnification(.2) + blend, 3)

    def test_invalid_parameters(self):
        with self.assertRaises(ValueError):
            exposure_flux(0, -1, 0, 1, .1, 1)


if __name__ == "__main__":
    unittest.main()
