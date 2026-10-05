"""Fit a point-source/point-lens model to public OGLE-IV photometry.

The output is deliberately browser-ready: observed fluxes, best-fit curves,
standardised residuals and information-criterion comparisons.  No parameter is
fitted in the browser, so the public page is reproducible from this script.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data" / "ogle-2023-blg-0100-phot.dat"
OUTPUT = ROOT / "data" / "ogle-2023-blg-0100-fit.json"
REFERENCE_MAG = 18.669


def magnification(u: np.ndarray) -> np.ndarray:
    return (u * u + 2.0) / (u * np.sqrt(u * u + 4.0))


def flux_model(time: np.ndarray, t0: float, t_e: float, u0: float,
               source_flux: float, blend_flux: float) -> np.ndarray:
    u = np.sqrt(u0 * u0 + ((time - t0) / t_e) ** 2)
    return source_flux * magnification(u) + blend_flux


def fit_model(time: np.ndarray, flux: np.ndarray, error: np.ndarray,
              blended: bool) -> dict:
    if blended:
        initial = np.array([2460001.5, np.log(24.0), np.log(0.31), 0.9, 0.1])
        lower = np.array([2459980.0, np.log(0.5), np.log(0.01), 0.0, 0.0])
        upper = np.array([2460020.0, np.log(300.0), np.log(3.0), 3.0, 3.0])

        def unpack(p):
            return p[0], np.exp(p[1]), np.exp(p[2]), p[3], p[4]
    else:
        initial = np.array([2460001.5, np.log(24.0), np.log(0.31), 1.0])
        lower = np.array([2459980.0, np.log(0.5), np.log(0.01), 0.0])
        upper = np.array([2460020.0, np.log(300.0), np.log(3.0), 3.0])

        def unpack(p):
            return p[0], np.exp(p[1]), np.exp(p[2]), p[3], 0.0

    result = least_squares(
        lambda p: (flux_model(time, *unpack(p)) - flux) / error,
        initial,
        bounds=(lower, upper),
        loss="linear",
        max_nfev=100_000,
    )
    parameters = unpack(result.x)
    prediction = flux_model(time, *parameters)
    residual = (flux - prediction) / error
    chi2 = float(np.sum(residual**2))
    k = len(result.x)
    n = len(time)
    return {
        "parameters": {
            "t0_hjd": float(parameters[0]),
            "tE_days": float(parameters[1]),
            "u0": float(parameters[2]),
            "source_flux": float(parameters[3]),
            "blend_flux": float(parameters[4]),
            "peak_magnification": float(magnification(np.array([parameters[2]]))[0]),
        },
        "chi2": chi2,
        "dof": n - k,
        "reduced_chi2": chi2 / (n - k),
        "aic": chi2 + 2 * k,
        "bic": chi2 + k * np.log(n),
        "prediction": prediction,
        "residual": residual,
    }


def main() -> None:
    raw = np.loadtxt(SOURCE)
    time, magnitude, magnitude_error, seeing, sky = raw.T
    flux = 10 ** (-0.4 * (magnitude - REFERENCE_MAG))
    flux_error = flux * np.log(10.0) * 0.4 * magnitude_error

    blended = fit_model(time, flux, flux_error, blended=True)
    unblended = fit_model(time, flux, flux_error, blended=False)

    weights = 1.0 / flux_error**2
    constant = float(np.sum(weights * flux) / np.sum(weights))
    constant_residual = (flux - constant) / flux_error
    constant_chi2 = float(np.sum(constant_residual**2))
    constant_bic = constant_chi2 + np.log(len(time))

    dense_time = np.linspace(
        blended["parameters"]["t0_hjd"] - 110,
        blended["parameters"]["t0_hjd"] + 110,
        700,
    )
    dense_flux = flux_model(
        dense_time,
        blended["parameters"]["t0_hjd"],
        blended["parameters"]["tE_days"],
        blended["parameters"]["u0"],
        blended["parameters"]["source_flux"],
        blended["parameters"]["blend_flux"],
    )

    t0 = blended["parameters"]["t0_hjd"]
    observations = [
        {
            "day": round(float(t - t0), 5),
            "hjd": round(float(t), 5),
            "magnitude": round(float(m), 4),
            "magnitude_error": round(float(me), 4),
            "relative_flux": round(float(f), 7),
            "flux_error": round(float(fe), 7),
            "seeing_px": round(float(s), 3),
            "sky": round(float(sk), 3),
            "model_flux": round(float(mod), 7),
            "standardised_residual": round(float(res), 5),
        }
        for t, m, me, f, fe, s, sk, mod, res in zip(
            time, magnitude, magnitude_error, flux, flux_error, seeing, sky,
            blended["prediction"], blended["residual"]
        )
    ]

    payload = {
        "schema_version": "2.0.0",
        "event": {
            "name": "OGLE-2023-BLG-0100",
            "field": "BLG510.09",
            "star_number": 99023,
            "ra_j2000": "17:59:29.62",
            "dec_j2000": "-35:27:55.6",
            "survey": "OGLE-IV Early Warning System",
            "band": "I",
            "reference_magnitude": REFERENCE_MAG,
            "n_observations": len(time),
        },
        "provenance": {
            "event_page": "https://ogle.astrouw.edu.pl/ogle4/ews/2023/blg-0100.html",
            "photometry_url": "https://www.astrouw.edu.pl/ogle/ogle4/ews/2023/blg-0100/phot.dat",
            "finding_chart_url": "https://ogle.astrouw.edu.pl/ogle4/ews/data/2023/blg-0100/fchart.jpg",
            "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
            "retrieved_utc": "2026-10-05",
            "columns": ["HJD", "I magnitude", "magnitude error", "seeing [px]", "sky level"],
        },
        "model": {
            "name": "point-source point-lens (PSPL) with constant blend flux",
            "equations": {
                "separation": "u(t)=sqrt(u0^2+((t-t0)/tE)^2)",
                "magnification": "A(u)=(u^2+2)/(u*sqrt(u^2+4))",
                "flux": "F(t)=Fs*A(u(t))+Fb",
            },
            "fit_method": "bounded scipy.optimize.least_squares, Gaussian reported flux errors",
        },
        "fits": {
            "constant": {
                "parameters": 1,
                "chi2": constant_chi2,
                "dof": len(time) - 1,
                "reduced_chi2": constant_chi2 / (len(time) - 1),
                "bic": constant_bic,
            },
            "pspl_unblended": {key: value for key, value in unblended.items() if key not in ("prediction", "residual")},
            "pspl_blended": {key: value for key, value in blended.items() if key not in ("prediction", "residual")},
        },
        "comparison": {
            "delta_bic_constant_minus_pspl": constant_bic - blended["bic"],
            "delta_bic_unblended_minus_blended": unblended["bic"] - blended["bic"],
            "outliers_abs_residual_gt_3": int(np.sum(np.abs(blended["residual"]) > 3)),
            "median_abs_standardised_residual": float(np.median(np.abs(blended["residual"]))),
        },
        "observations": observations,
        "curve": [
            {"day": round(float(t - t0), 5), "relative_flux": round(float(f), 7)}
            for t, f in zip(dense_time, dense_flux)
        ],
        "scope": {
            "inferred": ["event timescale", "impact parameter", "source and blend flux", "goodness of fit"],
            "not_inferred": ["unique lens mass", "lens distance", "planetary companion"],
            "reason": "A photometric PSPL light curve alone leaves the mass-distance-velocity degeneracy unresolved.",
        },
    }
    OUTPUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({
        "observations": len(time),
        "tE_days": blended["parameters"]["tE_days"],
        "u0": blended["parameters"]["u0"],
        "peak_magnification": blended["parameters"]["peak_magnification"],
        "reduced_chi2": blended["reduced_chi2"],
        "delta_bic_vs_constant": payload["comparison"]["delta_bic_constant_minus_pspl"],
        "delta_bic_blend": payload["comparison"]["delta_bic_unblended_minus_blended"],
        "outliers_gt_3sigma": payload["comparison"]["outliers_abs_residual_gt_3"],
    }, indent=2))


if __name__ == "__main__":
    main()
