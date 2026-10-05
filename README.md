# One Star, Briefly Brighter

A reproducible, beginner-readable microlensing case study built from **360 public OGLE-IV I-band observations** of **OGLE-2023-BLG-0100**. The site does not ask visitors to tune unexplained sliders. It leads them through three falsifiable questions:

1. Is the star temporarily brighter than its long-term baseline?
2. Does a point-source/point-lens (PSPL) curve describe the event?
3. What do the residuals and the decade-long baseline reveal?

## Why the event changed

The previous version used OGLE-2023-BLG-0001. The official OGLE record flags that object as **“Variable, Mira”**, so it was not a defensible clean microlensing anchor. This version replaces it with OGLE-2023-BLG-0100, whose EWS record has no such warning and reports a conventional point-lens solution.

## Result

The offline fit uses

\[
u(t)=\sqrt{u_0^2+\left(\frac{t-t_0}{t_E}\right)^2},\quad
A(u)=\frac{u^2+2}{u\sqrt{u^2+4}},\quad
F(t)=F_sA[u(t)]+F_b.
\]

The repository fit finds:

- Einstein timescale: **24.06 days** (OGLE EWS: 23.802 ± 2.584 days)
- closest approach: **u₀ = 0.307** (OGLE EWS: 0.312 ± 0.052)
- peak magnification: **3.37×**
- PSPL+blend reduced χ²: **1.58**
- seven measurements with |standardised residual| > 3
- the no-blend PSPL model has the lower BIC; the extra blend term is not justified by this implementation

The constant-brightness model is decisively worse (ΔBIC ≈ 12,778 relative to PSPL+blend). That supports a transient lensing-shaped excursion; it does not determine a unique lens mass.

## Data provenance

- Event record: [OGLE-2023-BLG-0100](https://ogle.astrouw.edu.pl/ogle4/ews/2023/blg-0100.html)
- Raw photometry: [official OGLE `phot.dat`](https://www.astrouw.edu.pl/ogle/ogle4/ews/2023/blg-0100/phot.dat)
- Columns: HJD, I magnitude, magnitude error, seeing in pixels, sky level
- Source SHA-256: `63f39592308728c46747400f5a608fb7ecf1367014ac54658221ccfdd2878a2f`
- Retrieved: 2026-10-05

The original file is preserved at `data/ogle-2023-blg-0100-phot.dat`. `scripts/fit_ogle_event.py` converts magnitude to relative flux, propagates the reported magnitude errors, fits constant/unblended/blended hypotheses, and writes `data/ogle-2023-blg-0100-fit.json`.

## Reproduce

```bash
python -m pip install -r requirements-research.txt
python scripts/fit_ogle_event.py
npm run check
npm run validate:research
```

Serve the repository through HTTP because the page fetches JSON:

```bash
python -m http.server 8000
```

## Scientific boundary

The photometric PSPL fit constrains the event shape: t₀, tE, u₀ and flux terms. It does **not** uniquely infer lens mass, distance or transverse speed. Breaking that degeneracy requires additional information such as microlens parallax, finite-source effects, astrometric lensing, or lens-flux constraints. The page makes that limit part of the main argument rather than hiding it in a disclaimer.

## Citation

Data: OGLE-IV Early Warning System. Survey reference: Udalski, A., Szymański, M. K., & Szymański, G. (2015), *OGLE-IV: Fourth Phase of the Optical Gravitational Lensing Experiment*, Acta Astronomica, 65, 1.

## Author

Biswajit Jana, 2026. MIT License.
