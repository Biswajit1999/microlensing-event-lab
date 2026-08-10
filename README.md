# Microlensing Event Lab

Paczynski light curves, blending and finite cadence residuals.

Created and maintained by Biswajit Jana.

## Scientific Purpose

This zero-build browser laboratory puts a compact reference-data bundle in front of the simulation. The app loads `data/reference.json`, renders those published anchors first, then sends the adjustable model to `physicsWorker.js` so numerical work stays off the UI thread.

## Architecture

- `index.html`: mission-control interface.
- `styles.css`: dense dark scientific dashboard.
- `app.js`: UI state, Canvas rendering and worker orchestration.
- `physicsWorker.js`: numerical model and heatmap generation.
- `data/reference.json`: small auditable reference-data bundle.
- `scripts/validate.js`: no-dependency repository validation.

## Run

```bash
python -m http.server 8080
```

Open `http://localhost:8080`.

## Validate

```bash
npm run check
```

The validation script checks required files, JSON reference data, worker syntax, citations and absence of unfinished scaffold tokens.

## Reference Data

Reference points following the Paczynski single-lens morphology used by OGLE microlensing analyses.

## References

- Paczynski, B., 1986. Gravitational microlensing by the galactic halo. The Astrophysical Journal, 304, pp.1-5.
- Udalski, A., Szymanski, M.K. and Szymanski, G., 2015. OGLE-IV: Fourth phase of the Optical Gravitational Lensing Experiment. Acta Astronomica, 65, pp.1-38.
