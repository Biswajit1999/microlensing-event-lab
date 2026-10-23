# OGLE-2023-BLG-0100: residual and robustness audit

Run from the repository root using the pinned OGLE photometry:

```bash
python -m pip install -r requirements-research.txt
python -m pip install pytest
python -m pytest -q tests/test_residual_robustness.py
python scripts/run_residual_robustness.py --output results/ogle_residual_audit.json
# Optional: pip install matplotlib; append --plot for a figure.
# Explicitly synthetic pilot: append --synthetic.
```

The script verifies the existing raw `phot.dat` SHA-256 before an observational
run and never modifies the original photometry or existing browser fit.

Three independent sensitivity checks:

1. **Irregular-cadence residual association:** Pearson correlation of
   adjacent standardized residuals separated by at most three days, plus a
   seeded shuffle benchmark under *iid exchangeability*. The permutation
   probability is not an event-anomaly false-alarm probability; overlapping
   adjacent pairs and observing conditions complicate interpretation.
2. **Seeing and sky dependence:** Multiple linear regression of PSPL
   standardized residuals on centered, unit-standard-deviation seeing and
   sky covariates. Report slopes and descriptive R²; no causal inference or
   automatic photometric correction is applied.
3. **PSPL fit robustness:** Compare the original Gaussian blended fit against
   soft-L1 loss, removal of |baseline residual|>3 points, and a *hypothetical*
   0.01 relative-flux error floor added in quadrature. All use the same
   parameter bounds as `scripts/fit_ogle_event.py`. Soft-L1 cost is not
   Gaussian chi-squared; clipping is data-dependent, so differences are
   sensitivity diagnostics rather than formal parameter errors.

The point-source point-lens model constrains event shape only. These tests
cannot identify the lens mass, detect a planet, or establish a photometric
systematic without independent calibration and follow-up. Cite the original
OGLE-IV EWS source and preserve its provenance in any publication or outreach.
