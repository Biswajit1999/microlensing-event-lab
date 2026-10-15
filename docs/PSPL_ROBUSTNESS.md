# OGLE-2023-BLG-0100: conditional PSPL robustness diagnostics

From the repository root, install `requirements-research.txt`, then run:

```bash
python -m unittest discover -s scripts -p 'test_pspl_robustness.py' -v
python scripts/pspl_robustness.py --bootstrap 64
```

The new script reads the existing OGLE-IV photometry and reuses `fit_ogle_event.py` without changing the browser's baseline JSON. It writes a **separate** `data/ogle-2023-blg-0100-robustness.json` file with three diagnostic families:

1. Exact analytic Jacobian of `F=Fs*A(sqrt(u0²+((t-t0)/tE)²))+Fb` and an SVD-based, local Fisher covariance in physical parameter units. Degenerate designs return no finite formal uncertainties.
2. Residual magnitude association with seeing and sky background, plus a top-seeing-quartile versus other squared standardized-residual comparison. These are descriptive, not causal tests; p-values are deliberately not presented.
3. Seeded Gaussian parametric bootstrap, refitting the original bounded blended PSPL model to each synthetic realization, reporting 16th–84th percentile intervals and fit failures.

**Limitations:** All uncertainty estimates are conditional on a PSPL light curve, the supplied photometric errors and independent Gaussian noise. The original fit has a reduced chi-square above unity; the formal Fisher and bootstrap intervals do not include correlated photometric systematics, selection effects, model mismatch, parallax or finite-source physics. These diagnostics do **not** establish a planetary lens, a unique lens mass, or a discovery. The OGLE-IV EWS photometry provenance remains in the original README.
