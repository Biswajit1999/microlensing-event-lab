"""Conditional point-source point-lens (PSPL) sensitivity calculations.

All times and errors are supplied by the caller. No OGLE data are loaded.
"""
import numpy as np


def magnification(u):
    u = np.asarray(u, dtype=float)
    if np.any(~np.isfinite(u)) or np.any(u <= 0):
        raise ValueError("separation must be finite and positive")
    return (u*u + 2) / (u * np.sqrt(u*u + 4))


def blended_flux(times, t0, u0, t_e, source_flux, blend_flux):
    times = np.asarray(times, dtype=float)
    if not np.all(np.isfinite(times)) or t_e <= 0 or u0 <= 0:
        raise ValueError("invalid times or lens parameters")
    u = np.hypot(u0, (times - t0) / t_e)
    return source_flux * magnification(u) + blend_flux


def numerical_jacobian(times, params, relative_step=1e-5):
    """Central differences for (t0,u0,tE,Fs,Fb), with physical constraints."""
    p = np.asarray(params, dtype=float)
    if p.shape != (5,) or not np.all(np.isfinite(p)) or p[1] <= 0 or p[2] <= 0:
        raise ValueError("invalid parameter vector")
    if not np.isfinite(relative_step) or not 0 < relative_step < 1e-2:
        raise ValueError("invalid step")
    cols = []
    for k in range(5):
        h = relative_step * max(abs(p[k]), 1.0)
        if k in (1, 2):
            h = min(h, 0.25 * p[k])
        plus, minus = p.copy(), p.copy()
        plus[k] += h
        minus[k] -= h
        cols.append((blended_flux(times, *plus) - blended_flux(times, *minus)) / (2*h))
    return np.stack(cols, axis=-1)


def fisher_diagnostics(times, sigma, params):
    """Weighted-Jacobian singular values and Fisher condition for independent Gaussian errors."""
    t = np.asarray(times, dtype=float)
    s = np.asarray(sigma, dtype=float)
    if t.ndim != 1 or t.size == 0 or s.shape != t.shape or not np.all(np.isfinite(s)) or np.any(s <= 0):
        raise ValueError("invalid uncertainties")
    j = numerical_jacobian(t, params)
    weighted = j / s[:, None]
    singular = np.linalg.svd(weighted, compute_uv=False)
    rank = int(np.linalg.matrix_rank(weighted))
    condition = np.inf if rank < weighted.shape[1] else (singular[0] / singular[-1]) ** 2
    return {"fisher": weighted.T @ weighted, "rank": rank,
            "singular_values": singular,
            "condition": float(condition)}


def cadence_deletion(times, sigma, params, groups):
    """Leave-one-group-out smallest weighted-Jacobian singular value."""
    t, s, g = np.asarray(times), np.asarray(sigma), np.asarray(groups)
    if g.shape != t.shape:
        raise ValueError("group shape mismatch")
    result = {}
    for label in np.unique(g):
        keep = g != label
        if np.count_nonzero(keep) < 5:
            raise ValueError("too few retained epochs")
        result[str(label)] = float(fisher_diagnostics(t[keep], s[keep], params)["singular_values"][-1])
    return result
