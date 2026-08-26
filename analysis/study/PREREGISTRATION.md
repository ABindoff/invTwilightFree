# Analysis plan — northern elephant seal reporting run

**Study version** `v1` · **Written** 2026-08-26, before the run · **Panel** 29
double-tagged deployments, 15 loggers, three seasons

This plan is written and committed before the run so the write-up cannot choose
its statistics after seeing them. The campaign has fitted these 29 deployments
many times under different scorers, kernels and calibrations; the manuscript will
quote this run and no other. Where a result contradicts a prediction below, the
prediction stands as written and the contradiction is reported.

---

## 1. What is fixed before the run

Every constant lives in `study_config.R`. Nothing downstream may redefine one.
Three decisions were settled on 2026-08-26 and are the reason this run exists.

### 1.1 Calibration

These tags were deployed to record depth. **There was no field calibration** — no
period on a post at a known position, before deployment or after. A purpose-built
GLS study would have roughly a week ashore at a known site at each end. This study
substitutes the pre-departure haul-out, and the manuscript says so in Methods and
again in Discussion.

The geometry window is `[t0, min(t0 + 7 d, departure)]`, where departure is dated
from the dive record alone — a hauled-out seal does not dive, so the last day
before continuous diving is the last day the release position is true. No Argos
enters this.

The 7-day cap is chosen *a priori* to match standard field practice and is **not
tuned against the Argos truth**. The truncation is what makes the window honest.
Measured on this panel:

| | days ashore |
|---|---|
| median | 6.8 |
| IQR | 3.2 – 11.8 |
| range | 0.0 – 16.9 |

At day 7 only 13 of 29 animals are still ashore, and the 90th percentile of
distance from the colony is already 537 km (max 664 km). A fixed window with the
position assumed is the recipe measured at 473 km against 237 km in
Supplementary S2.1.

Nine of 29 deployments have fewer than 200 usable pre-departure samples and cannot
calibrate themselves under any window. Pooling within tag family is therefore not
a refinement — it is what makes the panel analysable.

The intensity scale is each tag's own first 15 days, unchanged from the campaign,
because it needs a window whose shading regime matches the record being tracked
(S2.2) rather than one where the position is known.

### 1.2 Ground truth

**One scorer**: `argos_at()`, which returns `NA` when the bracketing Argos gap
exceeds 24 h. The permissive `approx(rule = 2)` variant that four campaign scripts
defined locally is forbidden — it extrapolates past both track ends and adds the
interpolator's error to the model's, preferentially where Argos coverage is worst.
Preflight asserts the gap rule is live.

`frac_scored` is reported next to every error. The dropped knots are not random:
a seal transmits when it surfaces, so long Argos gaps concentrate on the offshore
foraging leg. Errors from tags with very different scored fractions describe
different parts of a trip.

### 1.3 Endpoints

**Both ends pinned at Año Nuevo for every deployment.** For a central-place
forager tagged and recovered at one colony, that is prior knowledge. The
Argos-derived deploy and recover positions are recorded as diagnostics and are
never fitted. This removes the campaign's per-season endpoint rule, under which
6 of 19 deployments in 2022–23 were anchored at a position ≥ 20 km from the
colony (max 190 km) derived from the same Argos record used for scoring.

Pinning both ends suppresses the sloppy latitude direction roughly twentyfold
(fitted sensitivity to a calibration error 0.07 °/° against 1.4–4.4 °/° free).
**The panel is therefore a favourable case, and the reported errors are a floor
rather than a typical expectation.** The write-up states this.

---

### 1.4 The build, and four defects fixed before the run

The study runs against an installed build, asserted by `preflight.R` before
launch. Setting the run up surfaced four defects, three in the harness and one in
the package. All were found by a two-tag, forty-day smoke test that exercises
every arm in about seven minutes, and none would have announced itself in a
twelve-hour run — each produced plausible-looking output.

1. **Bathymetry sign (harness).** `bathy.tif` stores ETOPO altitude, negative
   below sea level, while `bathy_source()` and `floor_rule()` document
   positive-down depth. Passing altitude straight through made `floor_rule()`
   return `-Inf` for every ocean cell, so the hard arm excluded the entire ocean
   (`log_z = -Inf`, no scorable knot) and the soft arm penalised deep water about
   103 log-units more than the shore, driving the hierarchical fit 12.7° south.
2. **Endpoint versus hard floor (harness).** Año Nuevo's grid cell holds 16 m of
   water. A knot in which the animal dives 500 m excludes that cell, and since
   both endpoints are pinned there the whole track collapses to `log_z = -Inf` —
   at 66 of 80 knots in the smoke slice. The pinned endpoint knots are now exempt
   from the depth constraint, through `floor_rule()`'s own missing-data path. No
   auxiliary sensor can add information at a knot whose position is known, and
   one that contradicts it is wrong.
3. **Short checkpoint rows (harness).** A fit that scored nothing wrote a
   narrower CSV row, which made the checkpoint file unreadable and would have
   caused a resumed run to silently refit completed work. `summarise_knots()` now
   always returns the full schema.
4. **SMC auxiliary longitude indexing (package, now fixed).** Particle longitudes
   are held in −180..180 while the auxiliary raster carries the caller's
   convention. An analysis in 0–360 supplies an extent like 150..250, so
   `(p_lon - xmin) / cell_w` was negative for every particle and a negative
   `f64 as usize` saturates to 0 in Rust. **Every particle read column 0, and all
   sensor-fusion terms silently lost their entire longitude dependence.**
   Latitude still worked, which is what made the feature look alive: a
   latitude-only term moved the fit 38.5° → 48.9°, a longitude-only term left it
   bit-identical. This affects the guided filter as well as FFBS, so it reaches
   the engine the benchmark tables time. Fixed by wrapping the longitude into the
   raster's range before indexing, with `tests/testthat/test-smc_aux_longitude.R`
   pinning it in both conventions.

Defect 4 means the manuscript's §2.5 claim that sensor fusion reaches all three
engines was **false for the SMC engine in a 0–360 analysis** until this fix, and
any previously cached SMC fusion result is light-only in disguise.

## 2. Arms

| arm | engine | calibration | terms |
|---|---|---|---|
| `grid_light` | grid HMM | study | — |
| `grid_fusion` | grid HMM | study | hard bathymetry + coastline |
| `grid_cal_untrunc` | grid HMM | fixed 7 d, no truncation, pooled | — |
| `grid_cal_pertag` | grid HMM | fixed 7 d, no truncation, per tag | — |
| `hier_light` | hierarchical | study (panel median) | — |
| `hier_fusion` | hierarchical | study (panel median) | soft bathymetry + coastline |
| `ffbs_light` | FFBS, 3 tags | study | — |
| `ffbs_fusion` | FFBS, 3 tags | study | soft bathymetry + coastline |

Each contrast changes exactly one thing:

- `grid_light` − `grid_cal_untrunc` → **the departure truncation**
- `grid_cal_untrunc` − `grid_cal_pertag` → **pooling within tag family**
- `grid_fusion` − `grid_light` → **the tag's own depth channel**, in accuracy and
  in wall clock
- `hier_light` − `grid_light` → **joint estimation across the panel**
- `ffbs_light` − `grid_light` → **the cost of discretising space**

---

## 3. Endpoints and how they are analysed

**Primary.** Median great-circle error per deployment, km, `grid_light`.
Reported as the distribution across 29 deployments: median, IQR, range — not as a
single number.

**Secondary.** Latitude bias and RMSE; longitude bias and RMSE; empirical coverage
of the nominal 95% interval; wall clock per tag-day; `frac_scored`.

**Rules, fixed in advance.**

1. **Every within-panel contrast is paired by deployment** and tested with a
   Wilcoxon signed-rank test. Marginal summaries of two arms are reported
   alongside but never used as the test.
2. **When the mean and the median of per-tag values disagree in direction, both
   are reported.** On this panel they have done so before: the kernel correction
   moves the mean of per-tag medians 242 → 276 km (worse) and the median
   253 → 241 km (better).
3. **Clustered at the logger.** 29 deployments come from 15 loggers, three used
   three times. Every claim in this study is about the sensor, so the primary
   significance statement uses the logger as the unit. Deployment-level figures
   are reported beside them and are the optimistic ones.
4. **No claim rests on a single wall clock.** Timing is reported as seconds per
   tag-day with the raw seconds and the machine probe alongside.

**The season question is a null result unless it survives a test.** Season is
reported with a Kruskal–Wallis test and its effect size, and the three confounds
(Argos product, scored fraction, endpoint rule — the last now removed) are named.
On existing fits the gradient is not significant (p = 0.256, R² = 0.065), and the
default expectation is that it will not be here either.

---

## 4. Predictions

Recorded so the run can falsify them.

1. `grid_light` beats `grid_cal_untrunc` on median error, on a majority of tags,
   paired p < 0.05. *This is the calibration claim; if it fails, S2.1 does not
   generalise from ten tags to twenty-nine.*
2. `grid_cal_untrunc` beats `grid_cal_pertag`, and the gap is largest on the nine
   deployments that cannot calibrate themselves.
3. `grid_cal_pertag` fails outright on some deployments — they get no response and
   no track. That failure is the result, not a bug.
4. `grid_fusion` is **faster** than `grid_light` by roughly the fraction of cells
   the hard constraint excludes, and **not detectably more accurate**. The speed
   claim is the narrow, testable one.
5. Empirical coverage is well below the nominal 0.95, near 0.6–0.7, and is
   **bimodal across tags**. Simulation-based calibration passes on both arms, so a
   real-data shortfall is misspecification against real light residuals, not an
   implementation error.
6. Latitude bias is negative on a large majority of deployments, near −1.5°, with
   a spread of several degrees that **neither logger nor season explains**
   (measured between-logger ICC on existing fits: 0.00).
7. Season is not significant.

Prediction 4 and prediction 7 could each go the other way; 3 and 6 are stated
precisely so that a near miss counts as a miss.

---

## 5. Dropped, with reasons

**Pooled-evidence `z50` profiling.** Recommended earlier as the primary
calibration and withdrawn on 2026-08-26. It recovers a known `z50` exactly on the
synthetic battery, but `z50_profile_real.csv` — an incomplete run, 3 tags of a
planned 12 — puts the evidence maximiser 2° from the accuracy optimum on the two
tags complete across the sweep, at a cost of **+127 km**:

| dz50 (°) | pooled log_z | median km |
|---|---|---|
| 0 | −42.9 | **249** |
| +2 | **0.0** (max) | 377 |

That falsifies pre-registered prediction 4 of `diag_profile_recovery.R`, which the
synthetic had passed. n = 2, so it is not decisive; it is reported in the paper as
a negative result with its sample size, and the route is not used.

**Sea surface temperature.** Mk9 temperature as an SST proxy degrades latitude on
these tags (analysis §3.6). `USE_SST = FALSE`.

**Clock correction.** No per-tag estimate passes the acceptance rule (offset and
drift both within 20 min); applying them regardless costs +863 km and improved
0 of 10 tags. Raw timestamps throughout.

---

## 6. What is reported regardless of outcome

- The distribution across 29 deployments, not a headline number alone.
- `frac_scored` beside every error.
- The deployment-level scatter that nothing explains, as a stated finding.
- The design's favourability: both endpoints pinned at one colony.
- That these tags carried no field calibration, and the recommended protocol for
  a study designed for geolocation (about a week ashore at a known position at
  each end).
- Every arm that ran, including the controls that failed.
