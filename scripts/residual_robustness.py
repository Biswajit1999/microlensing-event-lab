"""PSPL residual structure and robust-fit sensitivity for OGLE-like light curves.

Three checks: irregular-cadence adjacent residual association, seeing/sky
covariates, and PSPL sensitivity to robust loss, clipping and error floors.
Not a lens-mass or planet-detection pipeline.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

from fit_ogle_event import flux_model, fit_model


def _vector(values: np.ndarray, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1 or len(arr) < 3 or not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must be a finite vector with >=3 elements")
    return arr


def adjacent_residual_association(time: np.ndarray, standardized_residual: np.ndarray,
                                  max_gap_days: float = 3.0, *,
                                  permutations: int = 1000,
                                  seed: int = 20261023) -> dict:
    """Lag-one Pearson r over only temporally adjacent close pairs.

    Permutation p is a descriptive iid-exchangeability benchmark, NOT a
    false-alarm probability for microlensing anomalies.
    """
    t = _vector(time, "time")
    r = _vector(standardized_residual, "residual")
    if t.shape != r.shape or not (np.isfinite(max_gap_days) and max_gap_days > 0):
        raise ValueError("Aligned vectors and positive maximum gap required")
    if not isinstance(permutations, int) or permutations < 1:
        raise ValueError("Need at least one permutation")
    order = np.argsort(t)
    t, r = t[order], r[order]
    gap = np.diff(t)
    if np.any(gap <= 0):
        raise ValueError("Times must be unique")
    mask = gap <= max_gap_days
    if np.count_nonzero(mask) < 4:
        raise ValueError("Too few adjacent close pairs")
    left, right = np.flatnonzero(mask), np.flatnonzero(mask) + 1

    def correlation(values: np.ndarray) -> float:
        a, b = values[left], values[right]
        aa, bb = a - np.mean(a), b - np.mean(b)
        den = np.linalg.norm(aa) * np.linalg.norm(bb)
        return float((aa @ bb) / den) if den > 0 else float("nan")

    observed = correlation(r)
    if not np.isfinite(observed):
        raise ValueError("Residuals have no measurable variation")
    rng = np.random.default_rng(seed)
    null = np.array([correlation(rng.permutation(r)) for _ in range(permutations)])
    null = null[np.isfinite(null)]
    if len(null) == 0:
        raise ValueError("No valid permutation correlations")
    empirical = (1 + int(np.sum(np.abs(null) >= abs(observed)))) / (len(null) + 1)
    return {"close_pair_count": int(np.count_nonzero(mask)), "max_gap_days": max_gap_days,
            "adjacent_pearson_r": observed, "iid_permutation_p_two_sided": empirical,
            "permutations_valid": len(null), "seed": seed}


def covariate_residual_association(residual: np.ndarray, seeing: np.ndarray,
                                   sky: np.ndarray) -> dict:
    """Descriptive OLS slopes per standard deviation of seeing and sky.

    Not a noise-systematics correction: seeing and sky may correlate with
    cadence, magnitude and observing season.
    """
    r, s, b = (_vector(x, name) for x, name in
               ((residual, "residual"), (seeing, "seeing"), (sky, "sky")))
    if not (r.shape == s.shape == b.shape):
        raise ValueError("Vectors must be aligned")
    std_s, std_b = np.std(s), np.std(b)
    if std_s <= 0 or std_b <= 0:
        raise ValueError("Covariates need nonzero variance")
    x = np.column_stack((np.ones(len(r)), (s - np.mean(s)) / std_s,
                         (b - np.mean(b)) / std_b))
    coef, _, rank, _ = np.linalg.lstsq(x, r, rcond=None)
    if rank < 3:
        raise ValueError("Covariates are collinear")
    prediction = x @ coef
    variance = float(np.sum((r - np.mean(r)) ** 2))
    r2 = 1 - float(np.sum((r - prediction) ** 2)) / variance if variance > 0 else 0.0
    return {"n": len(r), "intercept": float(coef[0]),
            "seeing_slope_per_sd": float(coef[1]),
            "sky_slope_per_sd": float(coef[2]),
            "descriptive_r_squared": float(r2)}


def _fit_with_error(time: np.ndarray, flux: np.ndarray, sigma: np.ndarray,
                    *, loss: str = "linear") -> dict:
    """Bounded blended PSPL fit with existing parameterization and priors."""
    initial = np.array([2460001.5, np.log(24.0), np.log(0.31), 0.9, 0.1])
    lower = np.array([2459980.0, np.log(0.5), np.log(0.01), 0.0, 0.0])
    upper = np.array([2460020.0, np.log(300.0), np.log(3.0), 3.0, 3.0])

    def unpack(p):
        return p[0], np.exp(p[1]), np.exp(p[2]), p[3], p[4]

    result = least_squares(lambda p: (flux_model(time, *unpack(p)) - flux) / sigma,
                           initial, bounds=(lower, upper), loss=loss,
                           f_scale=1.0, max_nfev=100_000)
    t0, te, u0, fs, fb = unpack(result.x)
    return {"success": bool(result.success), "t0_hjd": float(t0),
            "tE_days": float(te), "u0": float(u0),
            "source_flux": float(fs), "blend_flux": float(fb),
            "n": len(time), "cost": float(result.cost)}


def fit_sensitivity(time: np.ndarray, flux: np.ndarray, flux_error: np.ndarray,
                    *, clip_sigma: float = 3.0,
                    floor_fraction: float = 0.01) -> dict:
    """Compare baseline, soft-L1, baseline-residual clipping and flux floor.

    All fits use a blended PSPL and identical physical parameter bounds.
    The floor is in units of relative baseline flux and is hypothetical.
    """
    t, f, e = (_vector(x, name) for x, name in
               ((time, "time"), (flux, "flux"), (flux_error, "flux_error")))
    if not (t.shape == f.shape == e.shape) or np.any(e <= 0):
        raise ValueError("Aligned flux vectors and positive errors required")
    if not (clip_sigma > 0 and floor_fraction >= 0 and
            np.isfinite(clip_sigma) and np.isfinite(floor_fraction)):
        raise ValueError("Invalid sensitivity thresholds")
    baseline = fit_model(t, f, e, blended=True)
    mask = np.abs(baseline["residual"]) <= clip_sigma
    if np.sum(mask) < 10:
        raise ValueError("Too few points retained after clipping")
    robust = _fit_with_error(t, f, e, loss="soft_l1")
    clipped = _fit_with_error(t[mask], f[mask], e[mask])
    floor = _fit_with_error(t, f, np.sqrt(e*e + floor_fraction**2))
    standard = {k: float(baseline["parameters"][k]) for k in
                ("t0_hjd", "tE_days", "u0", "source_flux", "blend_flux")}
    return {"baseline": standard, "soft_l1": robust, "clipped": clipped,
            "floor": floor, "clip_sigma": clip_sigma,
            "clipped_points": int(np.sum(~mask)),
            "hypothetical_flux_error_floor": floor_fraction,
            "baseline_reduced_chi2": float(baseline["reduced_chi2"]),
            "note": "Clipping uses baseline residuals; robust cost is not Gaussian chi2."}
