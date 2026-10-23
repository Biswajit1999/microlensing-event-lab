"""Synthetic-only tests for OGLE residual and robust-fit diagnostics."""
import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from fit_ogle_event import flux_model
from residual_robustness import (adjacent_residual_association,
                                 covariate_residual_association, fit_sensitivity)


def lightcurve(seed=10):
    rng = np.random.default_rng(seed)
    t = np.linspace(2459955, 2460050, 130)
    f = flux_model(t, 2460001.5, 24., 0.31, .9, .1)
    e = np.full_like(t, 0.02)
    return t, f + rng.normal(0, .02, len(t)), e


def test_permutation_detects_synthetic_correlated_residuals():
    rng = np.random.default_rng(42)
    t = np.arange(150.)
    noise = rng.normal(size=150)
    residual = np.zeros(150)
    for i in range(1, 150):
        residual[i] = 0.75 * residual[i-1] + noise[i]
    report = adjacent_residual_association(t, residual, max_gap_days=1.1,
                                           permutations=400, seed=77)
    assert report["adjacent_pearson_r"] > .5
    assert report["iid_permutation_p_two_sided"] < .02
    assert report["close_pair_count"] == 149


def test_permutation_deterministic():
    t = np.arange(50.)
    r = np.sin(t * .3)
    a = adjacent_residual_association(t, r, permutations=100, seed=12)
    b = adjacent_residual_association(t, r, permutations=100, seed=12)
    assert a == b


def test_long_gaps_not_bridged():
    t = np.array([0., 1., 2., 3., 4., 100., 101., 102., 103., 104.])
    r = np.sin(t)
    report = adjacent_residual_association(t, r, max_gap_days=2)
    assert report["close_pair_count"] == 8


def test_duplicate_time_rejected():
    with pytest.raises(ValueError):
        adjacent_residual_association([0,1,1,2,3,4], [0,1,2,3,4,5])


def test_few_pairs_rejected():
    with pytest.raises(ValueError):
        adjacent_residual_association([0,5,10,15,20], [1,2,3,4,5], max_gap_days=1)


def test_covariate_regression_recovers_known_slopes():
    rng = np.random.default_rng(19)
    seeing = rng.normal(5, 1, 400)
    sky = rng.normal(400, 40, 400)
    residual = 0.5 + .8*(seeing-seeing.mean())/seeing.std() - .3*(sky-sky.mean())/sky.std()
    residual += rng.normal(0, .01, len(seeing))
    report = covariate_residual_association(residual, seeing, sky)
    assert report["seeing_slope_per_sd"] == pytest.approx(.8, abs=.01)
    assert report["sky_slope_per_sd"] == pytest.approx(-.3, abs=.01)
    assert report["descriptive_r_squared"] > .99


def test_covariate_collinearity_rejected():
    s = np.arange(10.)
    with pytest.raises(ValueError):
        covariate_residual_association(np.sin(s), s, 2*s)


def test_covariate_zero_variance_rejected():
    with pytest.raises(ValueError):
        covariate_residual_association(np.arange(10.), np.ones(10), np.arange(10.))


def test_baseline_fit_recovers_injected_parameters():
    t, f, e = lightcurve()
    report = fit_sensitivity(t, f, e)
    assert report["baseline"]["tE_days"] == pytest.approx(24, abs=3)
    assert report["baseline"]["u0"] == pytest.approx(.31, abs=.06)
    assert report["soft_l1"]["success"]
    assert report["clipped"]["success"]
    assert report["floor"]["success"]


def test_robust_fit_resists_single_extreme_outlier():
    t, f, e = lightcurve(20)
    f[len(f)//2] += 1.2
    report = fit_sensitivity(t, f, e)
    assert report["clipped_points"] >= 1
    assert abs(report["soft_l1"]["tE_days"]-24) < abs(report["baseline"]["tE_days"]-24)


def test_zero_floor_reproduces_baseline():
    t, f, e = lightcurve()
    report = fit_sensitivity(t, f, e, floor_fraction=0)
    assert report["floor"]["tE_days"] == pytest.approx(report["baseline"]["tE_days"], rel=1e-6)


def test_invalid_error_and_floor_rejected():
    t, f, e = lightcurve()
    e[0] = 0
    with pytest.raises(ValueError):
        fit_sensitivity(t, f, e)
    e[0] = .02
    with pytest.raises(ValueError):
        fit_sensitivity(t, f, e, floor_fraction=-.1)


def test_original_inputs_not_modified():
    t, f, e = lightcurve()
    original = [a.copy() for a in (t, f, e)]
    fit_sensitivity(t, f, e)
    for x, y in zip((t, f, e), original):
        np.testing.assert_array_equal(x, y)
