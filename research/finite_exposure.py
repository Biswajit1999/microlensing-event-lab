"""Finite-exposure point-source point-lens photometric diagnostics.

Times and exposure widths share units (normally days). The quadrature
models a uniform shutter response, not OGLE's calibrated exposure weighting.
"""
import math


def magnification(u):
    if not math.isfinite(u) or u <= 0:
        raise ValueError("separation must be finite and positive")
    return (u * u + 2) / (u * math.sqrt(u * u + 4))


def pspl_flux(t, t0, timescale, impact, source, blend=0):
    if not all(math.isfinite(x) for x in (t, t0, timescale, impact, source, blend)):
        raise ValueError("nonfinite model parameter")
    if timescale <= 0 or impact <= 0 or source < 0:
        raise ValueError("invalid timescale, impact parameter or source flux")
    u = math.hypot(impact, (t - t0) / timescale)
    return source * magnification(u) + blend


def exposure_flux(t, width, t0, timescale, impact, source, blend=0, n=64):
    """Composite midpoint quadrature of the flux over a centered exposure."""
    if not math.isfinite(width) or width < 0 or not isinstance(n, int) or n < 1:
        raise ValueError("invalid exposure width or quadrature order")
    if width == 0:
        return pspl_flux(t, t0, timescale, impact, source, blend)
    h = width / n
    return sum(pspl_flux(t - width / 2 + (k + 0.5) * h,
                         t0, timescale, impact, source, blend)
               for k in range(n)) / n


def quadrature_convergence(t, width, t0, timescale, impact, source, blend=0,
                           orders=(8, 16, 32, 64)):
    if len(orders) < 2 or any(not isinstance(n, int) or n < 1 for n in orders):
        raise ValueError("two or more positive integer orders required")
    values = [exposure_flux(t, width, t0, timescale, impact, source, blend, n)
              for n in orders]
    return values, [abs(b - a) for a, b in zip(values, values[1:])]


def blend_peak_degeneracy(peak_flux, baseline_flux, impact):
    """Solve Fs,Fb for a fixed impact parameter and measured peak/baseline.

    This is a conditional algebraic family, not a unique determination of u0.
    """
    if not all(math.isfinite(x) for x in (peak_flux, baseline_flux)):
        raise ValueError("nonfinite flux")
    a = magnification(impact)
    source = (peak_flux - baseline_flux) / (a - 1)
    blend = baseline_flux - source
    return source, blend
