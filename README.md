# invTwilightFree

`invTwilightFree` reconstructs animal tracks from the continuous light records
of archival geolocator (GLS) tags. It follows the *twilight-free* approach of
Bindoff et al. (2018): rather than detecting discrete sunrise/sunset events, it
models the likelihood of the entire light curve directly, which removes the
manual twilight-annotation step and makes the method robust to shading noise and
sensor anomalies. The numerical core is implemented in Rust (via `extendr`).

## What it does

The package provides three estimators over a shared likelihood and movement
model, in increasing order of scope:

* **`TwilightFreeSMC`** — a Sequential Monte Carlo (particle filter) for a single
  track. Two smoothers are available: a fast guided Brownian bridge (`"guided"`)
  and a Forward-Filtering Backward-Smoothing pass (`"ffbs"`) for full posterior
  smoothing.
* **`TwilightFreeGrid`** — a grid-based hidden Markov model on the same
  likelihood. Because it evaluates every cell explicitly, it is the natural place
  for **sensor fusion**: any number of auxiliary constraints (bathymetry, sea
  surface temperature, spatial masks or priors) can be added as terms in an
  additive log-likelihood.
* **`TwilightFreeHier`** — a partial-pooling hierarchical model over a *panel* of
  tags. Each track is reconstructed with a surrogate-posterior block sampler
  (delayed acceptance against the exact light likelihood), and the per-tag
  movement scale is pooled through a conjugate inverse-gamma / gamma population
  model, so short tracks borrow strength from the panel.

The likelihood is a spike-and-slab model of the light curve. The spike is an
asymmetric density around the light expected from solar geometry: its gentler
lower tail tolerates *shading*, which can only reduce measured light. The slab is
a flat component that absorbs gross false-light events unrelated to geometry, such
as artificial light at night or sensor glitches. Together they let the method
handle these anomalies without preprocessing, which is what underlies its
robustness under adverse conditions.

## Performance

The Rust engine makes the single-track filter fast enough for routine use. On a
180-day simulated deployment the guided filter completes in about 15 seconds with
1000 particles, against of order two hours for comparable SGAT/FLightR runs on the
same data. The engine is single-threaded, so this does not depend on core count.
It also degrades more gracefully than twilight-detection methods under heavy
shading, where discrete twilight events become unreliable.

On real data the grid HMM is the slower and more accurate engine: about 0.6 s per
12-hour knot, or 1.2 s per tag-day, measured over 29 elephant seal deployments.
Benchmarks and the code that produces them are in the package vignettes; the
comparator timings are being re-measured back to back on a quiesced machine before
publication, so treat the ratio above as indicative.

## Sensor fusion

`TwilightFreeGrid` and `TwilightFreeHier` accept auxiliary constraints through
`location_term()`, which pairs a data *source* with a *rule* mapping an
observation to a per-cell log-likelihood:

* **Sources** — `sst_source()`, `bathy_source()`, `raster_source()`,
  `sea_mask_source()`, `prior_raster()`, `function_source()`.
* **Rules** — `gaussian_rule()`, `student_rule()`, `mask_rule()`,
  `identity_rule()`, and threshold rules (`floor_rule()`, `ceiling_rule()`).

Each term contributes additively to the log-likelihood of every candidate cell,
so several sensors can be combined. In the hierarchical fit, `terms` may be a
single shared constraint, a list applied to every tag, or a
`function(id, df)` that builds per-tag terms (for example an SST term drawn from
each tag's own temperature column).

## Seasonal (hemisphere) priors

Near the equinoxes, day length alone cannot separate north from south and the
latitude likelihood is genuinely bimodal ("hemisphere swapping").
`hemisphere_prior()` encodes a date-dependent seasonal constraint: given a
function mapping each date to `"N"`, `"S"`, or `"both"`, it softly downweights
the wrong hemisphere. The default weighting is strong enough to break an equinox
tie yet weak enough to be overridden by decisive light data, and a migration
window can be left unconstrained.

The complementary equal-area (cosine-latitude) weighting of grid cells is no
longer a prior you supply. It belongs to the movement kernel, which is a density
on the sphere evaluated on a grid that is uniform in degrees, and it is applied by
`TwilightFreeGrid(area_correction = TRUE)` and by
`TwilightFreeHier(metric = "spherical")`, both defaults. The former `area_prior()`
term is deprecated; passing it as well applies the factor twice and is an error.

## Calibration verification

The samplers are checked by **simulation-based calibration** (SBC): data are
simulated from the model's own prior, refit, and the rank of each true parameter
within its posterior is accumulated. Correct posteriors give uniform ranks, which
is tested against simultaneous ECDF confidence bands. SBC scripts covering the
grid HMM, the hierarchical population parameters, and the spherical movement
metric are in `inst/sbc/` and `notes/da_block/`. This is a stronger guarantee
than a single point-estimate recovery check: it verifies that the *uncertainty*
is calibrated, not just that the mean is close.

## Installation

`invTwilightFree` compiles a Rust component, so it needs a working
[Rust toolchain](https://rustup.rs/) (rustc >= 1.65).

```r
# 1. Install Rust from https://rustup.rs/ for your OS.

# 2. Install the package.
# if (!requireNamespace("remotes")) install.packages("remotes")
remotes::install_github("ABindoff/invTwilightFree")
```

## Quick start

### Single track

```r
library(invTwilightFree)

# track_data has a Date/Time column and a light-intensity column.
fit <- TwilightFreeSMC(
  date_time   = track_data$time,
  light       = track_data$light,
  start_time  = as.POSIXct("2024-01-05"),  # optional: drop pre-deployment data
  end_time    = as.POSIXct("2024-06-15"),  # optional: drop post-retrieval data
  start_lat   = -45.0, start_lon = 140.0,
  end_lat     = -45.0, end_lon  = 140.0,
  method      = "guided",                  # or "ffbs" for full smoothing
  n_particles = 1000
)

plot(fit, type = "track")        # estimated path
plot(fit, type = "diagnostics")  # anomalous-light diagnostics
```

### Panel of tags (hierarchical)

```r
# `data` is a list of per-tag data frames (or one long frame with an id column).
# `locations` has one row per tag: id, deploy_lon, deploy_lat, and optionally
# retrieve_lon / retrieve_lat.
hfit <- TwilightFreeHier(
  data       = tracks,
  locations  = deployments,
  step_hours = 12,
  metric     = "spherical"
)

print(hfit)          # population movement scale + per-tag summaries
plot(hfit)           # per-tag movement posteriors, pooled
```

## Citation

If you use this method, please cite the foundational paper:

> Bindoff AD, Wotherspoon SJ, Guinet C, Hindell MA. Twilight-free geolocation
> from noisy light data. *Methods in Ecology and Evolution.* 2018;9(5):1190-1198.
> https://doi.org/10.1111/2041-210X.12953

## License

GPL (>= 3).
