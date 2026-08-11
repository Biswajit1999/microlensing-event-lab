# Microlensing Event Lab

Paczynski light curves, blending and finite cadence residuals.

Created and maintained by Biswajit Jana.

## Scientific Purpose

This zero-build browser laboratory puts a compact reference-data bundle in front of the simulation. The app loads `data/reference.json`, renders those published anchors first, then sends the adjustable model to `physicsWorker.js` so numerical work stays off the UI thread. It is a hands-on tool for building intuition about single-lens gravitational microlensing: how the Einstein-ring geometry of a foreground "lens" star turns into the smooth, symmetric brightening of a background source, and how survey-level effects like blending and cadence distort the observed curve.

## Background: Gravitational Microlensing

Gravitational microlensing happens when a massive foreground object (the **lens**) passes close to the line of sight to a more distant background star (the **source**). The lens's gravity bends the source's light, and because the light bending conserves surface brightness while changing the apparent solid angle of the source, the observer sees the source's flux *magnified* rather than displaced into resolvable multiple images (the images are typically far too close together, at milliarcsecond scales, to resolve from the ground).

Unlike lensing by galaxies or clusters, no persistent image splitting or arc is visible. The only observable is a smooth, achromatic, time-symmetric brightening and dimming of the source as the lens-source alignment sweeps past — a **microlensing event**. Because the effect depends only on mass and geometry, not on the lens emitting any light itself, microlensing is one of the few techniques that can detect completely dark objects: brown dwarfs, black holes, free-floating planets, and (via light curve anomalies) planets bound to the lens star.

### The Einstein radius and impact parameter

The characteristic angular scale of the effect is the **Einstein radius** `theta_E`, set by the lens mass `M`, and the distances to lens (`D_L`) and source (`D_S`):

```
theta_E = sqrt( (4 G M / c^2) * (D_S - D_L) / (D_L D_S) )
```

The lens-source separation is usually expressed in units of `theta_E` as the dimensionless **impact parameter** `u`. If the lens has relative proper motion `mu` with respect to the source, the source's line-of-sight offset evolves as:

```
u(t) = sqrt( u0^2 + ((t - t0) / tE)^2 )
```

where:
- `u0` — the minimum impact parameter, i.e. how close (in Einstein radii) the alignment gets at closest approach. Smaller `u0` means a sharper, higher peak.
- `t0` — the time of closest alignment (the peak of the event).
- `tE` — the **Einstein timescale**, `tE = theta_E / mu`, the characteristic duration of the event (days to months for stellar-mass lenses in the Galaxy).

### The Paczynski (1986) point-source, point-lens magnification

For a single point-mass lens and a point source, the total magnification of the two (unresolved) lensed images as a function of `u` is the Paczynski law:

```
A(u) = (u^2 + 2) / (u * sqrt(u^2 + 4))
```

`A(u) -> 1` as `u -> infinity` (no lensing, baseline flux) and `A(u) -> infinity` as `u -> 0` (formally infinite magnification for a true point source at exact alignment; in practice finite source size caps the peak). Combining the two equations above gives the full Paczynski light curve `A(t)` for a given `(t0, u0, tE)`.

### Blending

Ground-based microlensing fields are crowded, so the photometric aperture around the source star almost always contains additional, unrelated flux from unresolved neighbouring stars (and often the lens itself, if it is luminous). This constant extra light dilutes the observed magnification. If `f_S` is the true source flux and `f_B` is the blend flux, the observed relative flux is:

```
F(t) = (f_S * A(t) + f_B) / (f_S + f_B)
```

which can be rewritten with a blend fraction `b = f_B / (f_S + f_B)` as:

```
F(t) = 1 + b + (1 - b) * (A(t) - 1)
```

Blending suppresses the apparent peak magnification and is one of the main systematics that has to be modelled when fitting real survey photometry.

### Surveys: OGLE, MOA, KMTNet

Because any single microlensing event is unrepeatable and unpredictable (it depends on a chance alignment), detecting them requires continuously monitoring tens to hundreds of millions of stars, mainly toward the Galactic bulge, where source and lens star densities are highest. Dedicated survey collaborations do this:

- **OGLE** (Optical Gravitational Lensing Experiment, Warsaw University Observatory, Las Campanas, Chile) — the longest-running survey, now in its fourth phase (OGLE-IV), detecting thousands of events per year.
- **MOA** (Microlensing Observations in Astrophysics, Mt. John Observatory, New Zealand) — a Southern Hemisphere survey optimized for high-cadence bulge monitoring, providing longitude coverage complementary to OGLE.
- **KMTNet** (Korea Microlensing Telescope Network) — three identical telescopes in Chile, South Africa and Australia giving near-continuous, 24-hour coverage of the bulge fields, critical for catching short anomalies.

These surveys alert the community in near real time so that follow-up telescopes worldwide can densely sample events, in particular the short-lived light-curve anomalies that reveal a planet orbiting the lens star.

### Relevance to exoplanet detection

If the lens star hosts a planet, the planet's own much smaller Einstein radius produces a brief (hours-to-days), localized perturbation on top of the smooth single-lens curve as the source crosses near the planet's caustic structure. Because the perturbation's strength does not depend on the planet emitting light and the sensitivity peaks near the lens's projected Einstein radius (roughly the snow line, a few AU, for typical bulge-lens distances), microlensing is particularly good at finding:

- Cold, wide-orbit planets (beyond the snow line) that transit and radial-velocity surveys are far less sensitive to.
- Free-floating planets, seen as isolated short (hours-long) events with no bound host.
- Low-mass planets down to sub-Earth mass in favourable, high-magnification events.

This complements transit and radial-velocity surveys, which are biased toward close-in planets around bright, nearby stars, and is a large part of the science case for the *Nancy Grace Roman Space Telescope*'s planned Galactic Bulge Time Domain Survey.

## How It Works

This repository is a self-contained, zero-build browser app (no bundler, no framework):

1. **`app.js`** builds the control panel from the `LAB.controls` definition (`tE`, `u0`, `blend`, `cadence`), fetches `data/reference.json`, and on every slider change posts the current parameter set to a Web Worker so the UI thread never blocks on the numerics. It then draws the returned series and a parameter-response heatmap on two `<canvas>` elements and renders the scalar metrics panel.
2. **`physicsWorker.js`** runs off the main thread and implements the `microlensing(params)` model:
   - Builds a time axis spanning `+/-4 tE` around the peak.
   - Computes `u(t) = sqrt(u0^2 + (t/tE)^2)` and the Paczynski magnification `A(u) = (u^2+2) / (u*sqrt(u^2+4))`.
   - Applies the blend-fraction dilution `F(t) = 1 + blend + (1-blend) * (A(u(t)) - 1)`.
   - Reports scalar metrics: `peak_flux` (maximum of `F(t)`), `magnification` (`A(u0)`, the peak magnification at closest approach), and an approximate `fwhm_days` (full width at half of the excess magnification, from the analytic Paczynski FWHM relation).
   - Also generates a 72x72 normalized "magnification field" heatmap, sweeping impact-parameter offsets in both plotted dimensions with the same `A(u)` law, as a quick visual for how magnification falls off around the `(u0, 0)` alignment point.
   - The same worker file also hosts the physics models for this repository's sibling labs (CMB, supernovae, FRB dispersion, rotation curves, asteroseismology, weak lensing shear, spectrograph RV precision, clustering, exoplanet transmission spectra) behind a `lab` id dispatch table — only the `microlensing` branch is exercised by this app.
3. **`data/reference.json`** now ships the actual I-band photometry (160-epoch stratified subsample of a 631-point light curve) of a real OGLE-IV Early Warning System event, **OGLE-2023-BLG-0001**, pulled directly from `ogle.astrouw.edu.pl/ogle4/ews/`. Magnitudes are converted to baseline-normalised relative flux (`flux = 10^(-0.4*(mag - baseline))`, with the baseline taken as the median of the faintest 150 epochs), and time is referenced to the epoch of peak brightness. The real event shows genuine ~8x amplification at peak. Each point carries its propagated 1-sigma flux uncertainty, rendered as an error bar. This replaces the previous 5-point synthetic morphology sketch with an actual observed event — drag `tE`/`u0` to see how well a single-lens Paczynski curve can (or cannot) match real survey photometry.
4. **`research-overlay.js`** adds a non-invasive status panel that surfaces the repository's validation state, and **`data/research-reference.json`** plus **`scripts/validate_repository.mjs`** hold the auditable benchmark anchors used for that check (see `RESEARCH_QUALITY.md`).
5. **`scripts/validate.js`** is a dependency-free repository check (`npm run check`): verifies required files exist, that reference points are finite, that required citations appear in the README, and that no unfinished scaffold tokens remain.

The `cadence` control does not currently resample or add noise to the simulated curve inside `physicsWorker.js` — it is exposed as a parameter for future work on finite-cadence sampling residuals (see Limitations below).

## Architecture

- `index.html`: mission-control interface.
- `styles.css`: dense dark scientific dashboard.
- `app.js`: UI state, Canvas rendering and worker orchestration.
- `physicsWorker.js`: numerical model and heatmap generation.
- `data/reference.json`: small auditable reference-data bundle.
- `data/research-reference.json`: benchmark anchors for the research-quality validation layer.
- `research-overlay.js`: non-invasive validation/telemetry panel.
- `scripts/validate.js`: no-dependency repository validation.
- `scripts/validate_repository.mjs`: research-quality validation (see `RESEARCH_QUALITY.md`).

## Usage

Run a static file server from the repository root (no build step, no dependencies) and open the app in a browser:

```bash
python -m http.server 8080
```

Open `http://localhost:8080`.

Drag the **tE** (Einstein timescale), **u0** (impact parameter), **blend** (blend fraction) and **cadence** sliders to see the light curve, magnification-field heatmap and telemetry metrics update live via the Web Worker. The **Reset model** button restores the default parameters. The published-morphology reference points from `data/reference.json` are plotted as fixed yellow markers for comparison against the live model curve.

## Validate

```bash
npm run check
npm run validate:research
```

`npm run check` verifies required files, JSON reference data, worker syntax, citations and the absence of unfinished scaffold tokens. `npm run validate:research` runs the additional research-quality benchmark check described in `RESEARCH_QUALITY.md`.

## Math Appendix

Symbols: `t` observation time, `t0` time of peak alignment, `tE` Einstein timescale, `u0` minimum impact parameter (in units of the Einstein radius), `b` blend fraction.

Instantaneous impact parameter:

```
u(t) = sqrt( u0^2 + ((t - t0) / tE)^2 )
```

Paczynski (1986) point-lens, point-source magnification:

```
A(u) = (u^2 + 2) / (u * sqrt(u^2 + 4))
```

Blended observed relative flux (`b` = blend fraction of total baseline flux):

```
F(t) = 1 + b + (1 - b) * (A(u(t)) - 1)
```

Peak magnification at closest approach:

```
A_peak = A(u0)
```

These are exactly the relations implemented in the `microlensing()` function of `physicsWorker.js`.

## Reference Data

`data/reference.json` holds the actual 631-epoch I-band light curve of the real OGLE-IV
Early Warning System event **OGLE-2023-BLG-0001** (160-point stratified subsample, with real
per-point flux uncertainties), not an illustrative sketch -- see "Reference Data: A Real
Event" above for the full provenance.

## What the Timescale Tells Us: An Order-of-Magnitude Lens Mass

`bulge_lens_mass_Msun_estimate` converts the current `tE` slider value into a rough lens mass,
assuming typical Galactic-bulge lensing geometry (median relative parallax `pi_rel ~ 0.1 mas`
and relative proper motion `mu_rel ~ 4 mas/yr`, both representative OGLE-bulge-survey values;
Gould, 2000, ApJ, 542, 785): `M = (tE * mu_rel / 365.25)^2 / (kappa * pi_rel)`, `kappa = 8.14
mas/Msun`. At the default `tE = 24 d`, this lands around `0.08 Msun` -- a low-mass star or
brown-dwarf-class lens, consistent with the well-known statistical result that most Galactic
bulge microlensing lenses are low-mass M dwarfs (short `tE`), not massive stars.

**This is deliberately labelled an estimate, not a measurement.** A single-epoch Paczynski
light curve genuinely cannot determine lens mass uniquely: `tE` alone is degenerate with the
lens distance and the relative proper motion. Breaking that degeneracy in a real analysis
requires either microlensing parallax (a second vantage point, e.g. a satellite, or the annual
parallax signal in a long-duration event) or a resolved finite-source effect during the peak.
Reporting a plausible-but-uncertain number with its assumptions stated is more honest than
either omitting mass entirely or presenting a false-precision single value.

## Limitations

This is a teaching/exploration lab, not a research-grade fitting pipeline: it implements the
idealized point-source, point-lens Paczynski model only (no finite-source effects, no
parallax, no binary-lens/planetary caustics), and the `cadence` control is not yet wired into
resampling or synthetic noise. The lens-mass estimate above assumes typical bulge lensing
geometry rather than fitting it from the data, which is a real methodological limitation
inherent to single-epoch photometric microlensing, not an implementation shortcut.

## References

- Paczynski, B., 1986. Gravitational microlensing by the galactic halo. The Astrophysical Journal, 304, pp.1-5.
- Gould, A., 2000. A natural formalism for microlensing. The Astrophysical Journal, 542, pp.785-788.
- Udalski, A., Szymanski, M.K. and Szymanski, G., 2015. OGLE-IV: Fourth phase of the Optical Gravitational Lensing Experiment. Acta Astronomica, 65, pp.1-38.
- Gaudi, B.S., 2012. Microlensing surveys for exoplanets. Annual Review of Astronomy and Astrophysics, 50, pp.411-453.
- Mao, S. and Paczynski, B., 1991. Gravitational microlensing by double stars and planetary systems. The Astrophysical Journal, 374, pp.L37-L40.

## Research Quality Upgrade

See [RESEARCH_QUALITY.md](RESEARCH_QUALITY.md) for the validation layer, reference anchors, equations and research boundaries added to this repository.
