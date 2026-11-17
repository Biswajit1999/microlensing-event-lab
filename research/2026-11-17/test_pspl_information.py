import numpy as np
import pytest
from pspl_information import magnification, blended_flux, numerical_jacobian, fisher_diagnostics, cadence_deletion

P = (0.0, 0.3, 12.0, 2.0, 0.4)

def test_magnification_limit():
    assert magnification(1000.0) == pytest.approx(1, abs=1e-10)

def test_peak_flux():
    assert blended_flux([0.0], *P)[0] == pytest.approx(2*magnification(.3)+.4)

def test_flux_derivative():
    j = numerical_jacobian(np.array([-2., 0., 2.]), P)
    assert np.allclose(j[:, 4], 1)
    assert np.allclose(j[:, 3], magnification(np.hypot(.3, np.array([-2., 0., 2.])/12)), rtol=1e-7)

def test_information_scaling():
    t = np.linspace(-30, 30, 61)
    a = fisher_diagnostics(t, np.ones(t.size), P)
    b = fisher_diagnostics(t, np.full(t.size, 2.), P)
    assert np.allclose(b["fisher"], a["fisher"]/4)
    assert a["rank"] == 5

def test_cadence_groups():
    t = np.linspace(-30, 30, 61)
    out = cadence_deletion(t, np.ones(t.size), P, np.arange(t.size)%3)
    assert len(out) == 3 and min(out.values()) >= 0

def test_invalid_parameters():
    with pytest.raises(ValueError):
        magnification(0)
    with pytest.raises(ValueError):
        fisher_diagnostics([1,2], [1,0], P)
