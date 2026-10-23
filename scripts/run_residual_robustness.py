#!/usr/bin/env python3
"""OGLE residual diagnostics, pinned observational data or explicit synthetic fixture."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from fit_ogle_event import REFERENCE_MAG, fit_model, flux_model
from residual_robustness import (adjacent_residual_association,
                                 covariate_residual_association, fit_sensitivity)

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data" / "ogle-2023-blg-0100-phot.dat"
EXPECTED_SHA256 = "63f39592308728c46747400f5a608fb7ecf1367014ac54658221ccfdd2878a2f"


def synthetic_fixture():
    rng = np.random.default_rng(20261023)
    t = np.linspace(2459940, 2460080, 180)
    sigma = np.full(len(t), 0.025)
    flux = flux_model(t, 2460001.5, 24, 0.31, .9, .1)
    seeing = 5 + rng.normal(0, 0.8, len(t))
    sky = 400 + rng.normal(0, 45, len(t))
    noise = rng.normal(0, 1, len(t))
    for i in range(1, len(noise)):
        noise[i] += 0.25 * noise[i-1]
    flux += sigma * noise
    return t, flux, sigma, seeing, sky


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--synthetic", action="store_true", help="Explicitly use seeded simulated photometry")
    parser.add_argument("--output", type=Path, default=Path("results/ogle_residual_audit.json"))
    parser.add_argument("--plot", action="store_true", help="Requires matplotlib")
    args = parser.parse_args()
    if args.synthetic:
        t, flux, error, seeing, sky = synthetic_fixture()
        label, checksum = "Synthetic OGLE-like PSPL fixture", None
    else:
        digest = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
        if digest != EXPECTED_SHA256:
            raise ValueError("OGLE source SHA-256 mismatch; refusing to analyse altered input")
        t, mag, mag_err, seeing, sky = np.loadtxt(SOURCE).T
        flux = 10 ** (-0.4 * (mag - REFERENCE_MAG))
        error = flux * np.log(10.0) * 0.4 * mag_err
        label, checksum = "OGLE-2023-BLG-0100, OGLE-IV I band", digest
    baseline = fit_model(t, flux, error, blended=True)
    association = adjacent_residual_association(t, baseline["residual"],
        max_gap_days=3.0, permutations=1000, seed=20261023)
    covariates = covariate_residual_association(baseline["residual"], seeing, sky)
    sensitivity = fit_sensitivity(t, flux, error, clip_sigma=3.0,
                                  floor_fraction=0.01)
    report = {
        "dataset": label, "source_sha256": checksum,
        "n_observations": len(t),
        "baseline_reduced_chi2": float(baseline["reduced_chi2"]),
        "adjacent_residual_association": association,
        "seeing_sky_residual_association": covariates,
        "robust_pspl_sensitivity": sensitivity,
        "limitations": ["Permutation null assumes iid exchangeability and is not an anomaly false-alarm probability",
                        "Seeing/sky slopes are descriptive, not causal calibration coefficients",
                        "Clipping is selected from the fitted residuals and cannot establish formal uncertainty",
                        "Error floor is hypothetical, not inferred", "No planetary companion or lens mass inferred"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if args.plot:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.errorbar(t-baseline["parameters"]["t0_hjd"], baseline["residual"],
                    fmt=".", markersize=3, alpha=.65)
        ax.axhline(0, color="black", lw=.8)
        ax.axhline(3, color="gray", lw=.7, ls="--")
        ax.axhline(-3, color="gray", lw=.7, ls="--")
        ax.set_xlabel("Days from fitted peak")
        ax.set_ylabel("Standardized PSPL residual")
        ax.set_title(label + ": residual diagnostic")
        fig.tight_layout()
        fig.savefig(args.output.with_suffix(".png"), dpi=160)
        plt.close(fig)
    print(f"Wrote {args.output} for {len(t)} measurements")


if __name__ == "__main__":
    main()
