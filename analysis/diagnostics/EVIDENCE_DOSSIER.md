# Latitude bias in a light-based geolocation estimator: evidence dossier

Self-contained brief. No preferred hypothesis is stated, deliberately.

## The system

`invTwilightFree` estimates an animal's track from a light-level time series
recorded by an archival tag. Light intensity is compared against a modelled
clear-sky expectation as a function of solar zenith angle; position is inferred
from the timing (longitude) and shape/duration (latitude) of the daily light
cycle. No twilight thresholding: the whole light curve is used.

The engine under test is a **grid HMM**. Candidate positions are the cells of a
regular lon/lat grid (1 degree here). Latent state = cell occupied at each knot.
Knots are 12 hours apart. Forward-backward smoothing gives an exact marginal
posterior over cells at every knot. The reported point estimate is the joint
posterior mode; the posterior mean is also available.

Components:
- **Emission**: per observation, a spike-and-slab density. Spike = asymmetric
  exponential around the expected light, normalised over `[0, max_light]`. Slab =
  uniform, weight `prob_slab`, absorbing false-light events. The spike absorbs
  shading (an animal diving, cloud). Rate parameter `lambda`; upper arm rate =
  `shade_ratio * lambda`. About 24 observations per 12-hour knot.
- **Movement prior**: isotropic Gaussian in great-circle distance,
  `exp(-d^2 / 2 sigma^2)`, `sigma = diffusion * sqrt(dt_days)`, truncated at
  5 sigma. `diffusion` in km/sqrt(day), default 110 here.
- **Area correction** (on by default): `log(cos(latitude))` of the destination
  cell added to the forward variable once per knot, OUTSIDE the sum over source
  cells. Converts a density on the sphere to a probability per cell.
- **Light response**: expected light vs zenith. Here a 4-parameter logistic
  `floor + amp / (1 + exp((z - z50) / scale))`.

## The measurement apparatus (audit this if you are the adversarial reviewer)

`scratch/nes_calibration/build_light_battery.R` builds synthetic light records
with **exactly known truth**:
- true positions are real Argos fixes from northern elephant seal deployments,
  interpolated to observation times (30-minute spacing, ~240-day tracks);
- the `perfect` light level is the clear-sky logistic response evaluated at the
  true position at each observation time. **Noiseless. No attenuation.** All
  results below use this level;
- date is shifted in WHOLE DAYS by 0/90/180/270 to place identical movement in
  four different seasons while preserving local time of day;
- 6 tracks x 4 offsets = 24 records.

`fit_battery_null.R` then fits that light with **the same logistic response that
generated it**, and with both endpoints pinned to the true positions.

Truth is Argos-derived, but the light is GENERATED at those same Argos positions
and SCORED against them, so Argos error is common to signal and reference.

## The finding

Fitting noiseless, correctly-modelled light with the generating response yields a
**systematic southward (equatorward) latitude bias**.

48 of 48 fits:

| level | bias (by track) | rmse | median error | coverage |
|---|---|---|---|---|
| perfect (noiseless) | **-0.643 deg** | 1.03 | 73 km | 0.999 |
| noisy (clear sky + the engine's own spike-and-slab) | -0.455 deg | 1.92 | 208 km | 0.993 |

All 6 tracks negative (-0.46 to -0.81). Longitude is essentially unbiased.

For scale: on real data (29 deployments, against Argos) the latitude bias is
-1.54 deg. So this accounts for roughly 40% of it, with no shading, no response
misspecification and no clock error involved.

## Ruled out, with the method that ruled it out

| candidate | how it was tested | result |
|---|---|---|
| Argos error in truth | by construction: light generated AT the Argos positions, scored against them | cancels exactly; not a candidate |
| solar zenith algorithm | independent NOAA/Meeus implementation vs the package's USNO formula, 200k samples, 2021-23, 30-60 N | bias -0.0003 deg, rms 0.003, max 0.009. Also cancels by construction (same function generates and fits) |
| atmospheric refraction | the package computes GEOMETRIC zenith, no refraction | real: +0.221 deg in the twilight band (85-96), up to +0.578 at horizon. Cancels in this battery; a real-data issue only |
| point-estimate readout | compared joint mode / marginal-latitude mode / median / posterior mean | mode -0.639, marg mode -0.692, median -0.758, MEAN -0.797. The MEAN is MORE biased. Ordering explained by measured skewness -0.145. The posterior itself is displaced |
| grid quantisation | kernel drift computed at cell 1.0 and 0.5 deg | identical to 5 dp |
| cell-area correction | refit with `area_correction` off, 6 tracks | bias -0.567 -> -0.313 (+0.254, 6/6 tracks, p=0.031). Removing it moves estimates NORTH; it does not null the effect, it partly cancels it. Form is right (see below), magnitude is not |
| spherical movement-prior drift | derived, measured, corrected, retested | see next section. FALSIFIED |
| knot-grid alignment | knot phase shifted 0/3/6/9 h, identical data, common scoring window | bias -0.743/-0.738/-0.730/-0.739. Within-track spread median 0.01 deg, max 0.04. Fraction of knot boundaries landing in twilight varied 0.021-0.053 (2.5x) with no effect, so the test was sensitive. EXONERATED |

## The falsified hypothesis, in full (read this before proposing anything)

An isotropic kernel on a sphere pulls latitude equatorward before any data:
`sin(phi) = sin(phi0)cos(d) + cos(phi0)sin(d)cos(theta)`; averaging over uniform
bearing gives `E[sin phi] < sin(phi0)`. To second order, with `u = phi - phi0`,
`u ~= d cos(theta) - (1/2)tan(phi0)d^2 + (1/2)tan(phi0)d^2 cos^2(theta)`, and
`E[cos^2 theta] = 1/2` leaves **`E[u] ~= -tan(phi0) sigma^2 / (2R^2)`** per step.

Verified numerically against the engine's own kernel: correlation **0.99994**,
sigma^2 scaling **16.05** against 16 predicted, resolution-independent. Corroborated
at fit level: bias -0.567 at diffusion 110 vs -1.670 at 440 (6/6 tracks, p=0.031).

**It was implemented exactly and it changed nothing.** The correction nulls the
one-step drift to 5 dp. `log_z` moves, so it is applied. But the posterior MEAN
latitude moved **0.00009 deg** at diffusion 110 and **0.0008** at 440, against
biases of -0.586 and -2.025. Zero of 476 posterior modes changed. It is also inert
in the FORWARD FILTER alone (moves filtered bias by 0.023 deg), so smoother
cancellation is not the explanation either -- the light re-anchors the filter at
every knot and a 0.004 deg/step prior drift never accumulates.

**Two lessons that constrain what follows.**
1. The sigma^2 agreement was a coincidence. The real mechanism also widens with the
   kernel, so "worse at larger diffusion" is consistent with many stories and
   discriminates none.
2. A mechanism must be shown to survive into the SMOOTHED POSTERIOR, not merely to
   exist in a component, before a correction is designed for it.

## What the mechanism must look like (the strongest constraints available)

### It is seasonal. Track position explains nothing.

The battery places identical movement in four seasons. The bias profile MOVES with
the date offset:

| offset | 0-20% | 20-40% | 40-60% | 60-80% | 80-100% |
|---|---|---|---|---|---|
| 0 | 0.08 | -0.03 | **-1.65** | -1.17 | -0.88 |
| 90 | **-1.53** | -1.32 | -0.89 | **-1.75** | -0.87 |
| 180 | -0.88 | **-1.48** | -1.13 | 0.03 | 0.16 |
| 270 | -0.99 | 0.01 | 0.08 | -0.10 | **-1.76** |

Pooled over offsets, `e ~ frac` gives **R2 = 0.0003**.

### It is monotone in solar declination, and so is the posterior width.

| \|declination\| | bias | posterior sd |
|---|---|---|
| <5 deg (equinox) | **-1.672** | 1.68 |
| 5-10 | -1.188 | 1.45 |
| 10-15 | -0.790 | 1.19 |
| 15-20 | -0.534 | 1.00 |
| >20 deg (solstice) | **-0.388** | 0.91 |

Holds within every offset separately. Spearman(posterior sd, |declination|) =
**-0.807**. Context: latitude is read from the light's even harmonics, whose
content nulls at the equinox, so the light genuinely goes blind to latitude there.

### The bias is proportional to posterior VARIANCE, with a form that matches and a magnitude that does not.

A log-linear tilt `c` on a knot's marginal shifts its mean by `c * sigma_post^2`.
For the area factor `c = -tan(lat)`.

| model | R2 |
|---|---|
| `e ~ -sigma_post^2 tan(lat)` (the area tilt) | **0.6848** |
| `e ~ sigma_post` | 0.6061 |
| `e ~ \|declination\|` | 0.2681 |
| `e ~ frac` | 0.0003 |

Spearman(e, tilt) = +0.768. **But the fitted slope is +48.9 where 1.0 would be
exact**, and the ratio of means is 30.5x. Mean predicted tilt -0.0263 deg vs mean
observed bias -0.8024 deg.

So: a systematic equatorward FORCE, resisted by the likelihood, winning wherever
the light is weak. It has the area tilt's functional form and roughly **30-50 times
its strength**.

One structural note offered as fact, not interpretation: there are about **24
observations per 12-hour knot**, and the area factor is applied ONCE per knot while
emission terms are applied PER OBSERVATION.

## Open questions

1. What produces an equatorward force of order 30-50x `tan(lat)` per unit latitude
   (i.e. roughly 0.5-0.9 nats per degree of latitude) in this model?
2. Is the apparatus sound? Every conclusion above rests on
   `build_light_battery.R`. It has not been independently audited.
3. Anything that predicts a result the leading alternatives do NOT predict.

## Relevant files

- `src/rust/src/lib.rs` - engine. `run_grid_hmm`, `log_move_kernel`,
  `log_cell_area`, `ito_lat_shift`, `spike_density`, `spike_normaliser`,
  `expected_light`, `LightMix`.
- `R/TwilightFreeGrid.R` - R interface, knot construction (line ~220),
  `grid_posterior()`.
- `scratch/nes_calibration/build_light_battery.R` - the apparatus.
- `notes/latitude_bias_investigation.md` - narrative record.
- Result CSVs in `scratch/nes_calibration/`: `battery_null.csv`,
  `mode_vs_mean.csv`, `area_diffusion.csv`, `equinox_tilt.csv`,
  `knot_phase.csv`, `filtered_vs_smoothed.csv`.

## What a useful answer looks like

A DISCRIMINATING EXPERIMENT: a proposal that predicts an outcome the leading
alternatives do not. Every hypothesis that has failed here fit the existing
numbers well, including matching sign, scaling and interaction terms. Fitting the
numbers is not evidence. Cheap experiments are strongly preferred; a full fit is
~2.5 minutes per record and pure numerics are free.
