# Calibration-window investigation, 2026-08-06

Scratch diagnostics behind commits 177d934, 2f28746 and 2bc32e0. Gitignored, as
everything under `scratch/` is; kept because the negative results cost more to
find than the positive one and are easy to re-derive by accident.

Run from the package root. `nes_common.R` carries the crosswalk, the Argos
reader, the departure detector and the scoring; every other script sources it.
All of them need `data/nes_untracked/` and the cached archives.

## What each one settled

| script | question | answer |
|---|---|---|
| `diag_overcount.R` | does a composite-likelihood weight fix latitude coverage? | **No.** sd grows as kappa^0.26 not sqrt(kappa), and the point estimate moves 5.7 deg. Posterior sd already matched error sd, so the miscalibration is a heavy tail, not a narrow interval. |
| `diag_daylength.R` | is the bias a day-length error? | Inconclusive and **biased by its own design** — thresholding each day at half of THAT day's amplitude makes a clouded day cross later and manufactures the effect. Superseded. |
| `diag_daylength2.R` | same, on the response's absolute half level | Fitting the response at Argos positions instead of the colony halves the day-length scatter (sd 7.2 -> 3.5 deg of latitude). First real evidence that the WINDOW is the problem. |
| `diag_calwindow.R` | can the haul-out replace the colony assumption? | Partly. Rescues 2021025 (543 -> 168 km) and wrecks 2021027 (211 -> 843). Geometry yes, scale no. |
| `diag_calsplit.R` | take the scale from the whole record instead? | **No.** Every record-scale variant drifts north (486, 331, 364 km). Deep dives drag q05 down, the light looks brighter, days look longer. |
| `diag_scale.R` | use the envelope's floor/amp as the scale? | **No**, and now tested under a CORRECT calibration too: 572-894 km, biased south. Vindicates the existing q05 choice. |
| `diag_bootstrap.R` | bootstrap geometry from a first-pass track? | **No.** Diverges monotonically, 321 -> 391 -> 755 km, z50 marching down each turn, even confined to the calibration window. |
| `diag_split2.R` | haul-out geometry + 15-day scale | **Yes.** 321 -> 104 km on three tags. The combination the earlier failures pointed at. |
| `diag_all10.R` | does it hold on all ten, full deployments? | Yes where it applies: 508 -> 259 km on the seven with a haul-out. Three animals left within a day of tagging and have none. |
| `diag_pool.R` | pool geometry across tags of one model? | **Yes, and it beats per-tag fitting** (209 vs 259 km) even for tags that have their own haul-out. 473 -> 237 km over all ten. |
| `diag_slope.R` | tangent amplitude: max_light or the envelope's amp? | max_light, 259 vs 411 km, once the windows differ. |
| `diag_default.R` | does that hold on the default path too? | Yes, 473 -> 363 km, so the conditional was removed. |
| `diag_report_config.R` | does any of it survive the report's 1-degree grid and diffusion 110, with raw timestamps? | **Yes, all of it.** 429 -> 208 km. See below. |

`smoke_chunk.R` runs the rewritten section 5.2 without any fitting.
`report_results.R` summarises `analysis/output/nes/*.csv`.
`show_clocks.R` prints the estimated clock offsets from the cache.

## Resolved: it was the clock, and only the clock

The full render came back at 557-857 km against 224 km previously, and the
report's own naive-versus-pooled comparison was a wash (704 vs 730) where this
harness gave 473 vs 237 on the same animals. That looked like the calibration
change failing at the report's resolution.

It was not. `diag_report_config.R`, in the report's own grid, domain and
diffusion, with RAW timestamps, crossing the two halves of the change
(`RESULT_report_config.txt`):

| recipe | median km | worst | q90 | bias_lon | rmse_lat |
|---|---|---|---|---|---|
| a 15d geometry, amp slope (previous release) | 429 | 1038 | 788 | +0.20 | 4.46 |
| b 15d geometry, max_light slope | 318 | 631 | 630 | +0.25 | 3.56 |
| c pooled geometry, amp slope | 345 | 445 | 562 | +0.27 | 3.41 |
| **d pooled geometry, max_light slope (current)** | **208** | **396** | **447** | +0.30 | **2.44** |

Both changes help, independently and together, at 1 degree as at 2. And
`bias_lon` is +0.2 to +0.3 in EVERY recipe, so the +2.5 to +7 degrees seen in
the render came entirely from the clock correction.

The mechanism: `calibrate_clock_from_endpoints()` calls `fit_light_response()`.
Its profile likelihood used to run to the +/-120 min search bound on all ten
tags, and the acceptance rule (reject above 60 min) threw every estimate away.
Better calibration gave the profile interior optima, which were then accepted and
applied -- one tag got 9 min at deployment ramping to 41.7 at recovery, and
`correct()` spreads that linearly over the record, about 8 degrees of longitude
by the end of the trip.

The rule was never a plausibility test, only a did-it-hit-the-wall test. The
report's own SPOT table says a clock of this vintage drifts between -1.4 and
+13.3 minutes over a whole deployment, so `clock_max_accept_min` is now 20 and
applies to the implied DRIFT as well as the offset. `show_clocks.R` prints the
estimates from `analysis/cache/nes/clocks_v2.rds` if this needs rechecking.

Lesson worth keeping: fixing one component moved a second component from
failing-loudly to failing-quietly. The estimates were always bad; they only
became harmful once they stopped hitting the bound.

## Second round: leakage, dilation, and the depth channel

| script | question | answer |
|---|---|---|
| `diag_leakage.R` | where else could Argos be helping rather than only scoring? | Recovery endpoint is **not** leakage (median 1 km from the deployment fix -- it IS the colony). Mesh pad is soft: 74/26 deg against an observed range of 65/19. |
| `diag_domain.R` | is the grid built from the Argos envelope flattering the result? | **No.** Three domains -- Argos-derived, a truth-free North Pacific box, and a deliberately over-wide one -- give **bit-identical** fits (208 km, rmse_lat 2.439, max difference 0.000). Fitted latitude spans 33.5-56.5 against bounds of 20-70, and no knot ever lands within a cell of an edge. The box never binds. Domain changed to the truth-free one anyway, because an argued domain beats a harmless one. |
| `diag_dilation.R` | is the loop dilation a declination sign-flip? | **No, falsified.** At a self-crossing (15474 pairs, same place >60 days apart) the latitude-error difference is uncorrelated with the declination difference (r = -0.001, p = 0.88). Speed is no better (p = 0.53). |
| `diag_dilation2.R` | then what is it? | **Noise, not deformation.** Latitude span of the estimate is 1.99x truth's at matched longitude, but est-on-true regression slope is 0.85-1.02 (*compression*), and only 5 of 10 tags have opposite-signed limb bias (mean +/- difference -0.12 deg). Error sd 3.11 deg against a true latitude span of 2.9 deg: the noise is as large as the structure, so the eye reads two independent perturbations on a closed loop as a coherent deformation. |
| `diag_dive_movement.R` | does dive effort predict step length? | Pooled correlation **+0.013** -- looks like nothing. |
| `diag_dive_movement2.R` | ...after removing the haul-out? | **Simpson's paradox.** At-sea only (94% of knots) the correlation is **-0.170**, step length falls 31.9 -> 24.6 km monotonically across dive-depth quintiles, p = 2.6e-43. Real, but worth only 1.33x on sigma and a 2.5% reduction in residual sd -- against an 8x insensitivity of accuracy to D. Good biology, poor geolocation lever. |
| `diag_states.R` | is latent state-switching worth having? | **The data want it; geolocation does not.** d_log_z +156 (2-state 44/220) and +431 (30/300), better on 10/10 tags. Accuracy flat (-5 to +11 km). Coverage gets WORSE, 0.77 -> 0.65, because switching tightens an already-too-narrow posterior. 41 of 60 fits; the last two configs were skipped as redundant and can be resumed. |
| `diag_residence.R` | what is a guided/Langevin attraction worth? | **Promising, not established.** Weight 0.25: 209 -> 180 km, rmse_lat 2.44 -> 1.90, lat_span 14.3 -> 12.3 against truth's 10.8 (information, not shrinkage). But **p = 0.375** paired over 10 tags, and it damages the individuals that deviate from the population. |

### State-switching: the accuracy/calibration trade

| config | median km | rmse_lat | post_sd_lat | cover_lat | mean log_z |
|---|---|---|---|---|---|
| 1-state D=44 (Argos one-step scale) | **190** | **1.82** | 0.61 | 0.44 | -49863 |
| 2-state 44/220 | 203 | 2.10 | 1.05 | 0.65 | -49785 |
| 1-state D=110 (report) | 209 | 2.44 | 1.42 | **0.77** | -49941 |
| 2-state 30/300 | 220 | 2.03 | 1.02 | 0.66 | **-49510** |

The single-state fit at the TRUE one-step scale gives the best accuracy in the
sweep and the worst calibration. Accuracy and coverage pull opposite ways
because the residual error is per-tag latitude BIAS, which no movement model can
represent: D = 110 is not "right", it is loose enough to swallow a bias the model
does not know it has. **So the coverage problem is not under-dispersion to be
fixed by a better movement model.** It is unmodelled bias papered over by an
over-loose prior, and it points back at the observation model again.

Corollary for the dive-effort covariate: state-switching structure is real
(log_z is emphatic) but it is not the geolocation bottleneck, so a covariate that
improves state identification improves the wrong component.

### Residence prior: helps the weak, hurts the distinctive

| weight | median km | q90 | bias_lat | rmse_lat | rmse_lon | cover_lat | lat_span |
|---|---|---|---|---|---|---|---|
| 0.00 | 209 | 448 | +0.66 | 2.44 | 0.75 | 0.77 | 14.3 |
| 0.25 | **180** | **353** | +0.25 | **1.90** | 0.81 | 0.78 | 12.3 |
| 0.50 | 194 | 382 | +0.12 | 1.99 | 0.95 | 0.72 | 11.3 |
| 1.00 | 250 | 663 | -0.08 | 2.25 | 1.42 | 0.66 | 10.5 |

True latitude span is 10.8, so weight 0.25 removes excess spread without
over-shrinking and weight 1.00 crosses the line (span 10.5, error up, longitude
wrecked). The `lat_span` guard did its job.

Per tag at weight 0.25: 7/10 better on median error, 8/10 on latitude RMSE,
**paired Wilcoxon p = 0.375**. The gain is concentrated in the animals the light
was failing (2021034 -181 km, 2021035 -88) and the losses in the ones it was not
(**2021033 +121 km**, 2021023 +69). 2021033 is the most northerly track in the
panel at 56 N, and a field built from the other nine does not cover where it
went, so it is dragged toward the population. **A population prior penalises the
individuals who deviate from the population** -- often the reason for tagging.

What it would need: gentle weighting; ANISOTROPY, since it is isotropic and
spends longitude accuracy that was not needed (rmse_lon 0.75 -> 1.42 at weight
1); and per-individual deviation rather than a fixed field, which is partial
pooling applied to occupancy instead of to movement scale.

Read all of it as an UPPER BOUND: the field is built from other animals' ARGOS
tracks, not from light-based ones.

### Per-tag observation noise: falsified, and the default vindicated

`diag_lambda.R`. Latitude scatter varies fourfold between tags (1.47 to 6.16 deg)
while the model hands every animal the same uncertainty, so the obvious fix is a
per-tag noise parameter. A tempering constant is not identifiable from the
likelihood it tempers, so the candidate was the spike rate `lambda`, which is a
parameter INSIDE the likelihood -- the argument being that log_z therefore has an
interior optimum and can select it. Sweep of a multiplier on the default
`1/(max_light*0.5)`:

| mult | median km | rmse_lat | post_sd_lat | cover_lat | mean log_z |
|---|---|---|---|---|---|
| 0.25 | 382 | 6.85 | 2.19 | 0.62 | -54697 |
| 0.50 | 228 | 3.38 | 1.76 | **0.78** | -52773 |
| **1.00 (default)** | **209** | **2.44** | 1.42 | 0.77 | -49941 |
| 2.00 | 250 | 2.76 | 1.16 | 0.58 | -47067 |
| 4.00 | 283 | 3.40 | 0.96 | 0.42 | **-46065** |

Three pre-registered criteria, all failed. **No interior optimum: 0 of 10 tags**
-- every animal picks the largest lambda by 2886-4758 log units. **No per-tag
variation**: all ten select the same value. **Selecting by log_z is worse**:
209 -> 283 km, coverage 0.77 -> 0.42.

WHY THE ARGUMENT WAS WRONG, since it looked sound. The claim was that as lambda
grows the spike becomes a delta, every real observation falls to the slab, and
information is lost -- so log_z must turn over. But the TRACK is a latent state
with freedom to move: a sharper spike lets the model find positions where the
expected curve passes close to the observations, and the peak density grows
faster than the misfit penalises. **Marginal likelihood is not a safe criterion
for the observation-noise scale when the latent state can absorb the residuals.**

Two things salvaged from it. `lambda` DOES control interval width as intended
(posterior sd 2.19 -> 0.96, monotone), so the mechanism is sound and only the
selection rule fails. And the hard-coded default, never previously tested, turns
out to be near-optimal on both axes at once: best for accuracy and within 0.01
of best for coverage.

### The coverage problem: three failures with one cause

`overcount`, latent state-switching and per-tag `lambda` have now all failed to
fix per-tag interval calibration, and for the same reason: **nothing in the light
tells you how wrong the light is.** Position-free noise indices reach only
r ~ 0.54, p = 0.11 (`diag_noiseindex.R`). This looks like a real limit rather
than three dead ends, and it is worth stating as such: a model of this kind
cannot self-diagnose its own reliability, so honest PER-TAG uncertainty needs
external information. That is the strongest argument for double-tagging a subset
-- stronger than the occupancy prior, because this gap cannot be closed from the
light at all.

2021032 is untouched by every one of them (coverage 0.13 at the default, 0.36,
0.31, 0.27 across the lambda sweep, 0.29 under switching). Whatever is wrong
with that animal is not calibration, not movement and not observation noise.

### 2021032: seven mechanisms excluded, cause unknown

`diag_alan.R`, `diag_truegeom.R`, `diag_2021032.R`, `diag_seasonal.R`.

The persistent outlier: coverage 0.13, latitude scatter 6.16 deg against 1.16-3.11
for the rest, and unmoved by recalibration, the movement model, state-switching,
the lambda sweep and the search domain.

**Not ALAN.** Tested because the North Pacific Transition Zone carries a large
squid-jigging fleet. 2021032 is the DARKEST tag at night of the ten: night light
above its own dark level `frac_above_10` = 0.058 and `frac_above_30` = 0.000,
both lowest, max excess 27 against 36-50 elsewhere. Correlations run against the
hypothesis (more night light, LESS scatter, r = -0.42).

**Not its response calibration.** Fitted against its true Argos positions over
230 days at sea, its z50 is 92.26 against the pooled 92.09 -- an error of +0.17
deg, the SMALLEST of the ten. The pooled curve fits this animal better than any
other.

**The failure is seasonal, not episodic.** Monthly latitude bias: -3.8, -4.8,
-3.2 through Jul-Sep, then +8.5, +10.8, +8.2 through Oct-Dec. A 12 degree
reversal at the equinox, more than twice any other tag's. Only 28% of its squared
error lies in its worst tenth of knots (second-LOWEST concentration in the panel),
so this is two sustained opposite biases, not a bad patch.

**Its twilight light got noisier in the second half**: twilight scatter 0.161-0.164
(Jul-Sep) rising to 0.215-0.224 (Oct-Dec), while a clean control stays 0.165-0.193
throughout. Sensor fouling was checked and is NOT indicated -- daytime maximum
declines no faster than the control's.

Cause unidentified. Further sweeps on one animal with no candidate mechanism
left would be fishing.

### The declination mechanism: proposed, falsified, resurrected, falsified again

Worth recording in full because the idea is seductive and will occur to the next
person. A fixed error in the assumed threshold zenith produces opposite-signed
latitude bias either side of the equinox, because
`f'(phi) = cos(phi) sin(delta) - sin(phi) cos(delta) cos(H)` changes sign with
declination. It predicts that a loop track comes back DILATED, which is exactly
what the figures show.

1. Tested pooled across tags and at self-crossings: r = -0.001, p = 0.88. Dead.
2. Resurrected on 2021032, which shows the pattern unmistakably (-3.7 summer,
   +8.0 winter) and whose pooled transition width is 4.7 deg too wide, putting
   its assumed zero-crossing about 2 deg too deep -- the exact error required.
3. Tested per tag with the signed, quantitative prediction (swing should be
   positive whenever the pooled zero-crossing is deeper than the tag's true one,
   and should scale with it): **7 of 9 tags swing the WRONG WAY**, r = +0.459,
   p = 0.21, and 2021033 has the largest zero-crossing error with the smallest
   swing. Dead again.

Two of ten animals behave as the mechanism says; the panel norm is the opposite
sign. It is not the general driver of the dilation, whatever it does for 2021032.

### Why pooling beats per-tag fitting: the measurement is noisier than the signal

`diag_truegeom.R`. Comparing each tag's haul-out calibration with the response it
really had at sea:

- haul-out z50 vs at-sea z50: **r = -0.43**; width: **r = -0.57** (both negative)
- mean |haul-out z50 minus at-sea z50| = 0.87 deg
- total spread of true at-sea z50 across all ten animals = 1.85 deg (90.66-92.51)

The per-tag measurement error is comparable to the entire between-animal
variation, so per-tag estimates are noise around a common value. That is why
pooling beat per-tag fitting (209 km against 259), and it means **partial pooling
would not help either**: with signal-to-noise below 1, full pooling is the right
answer rather than a compromise.

### Ranking of everything tried, by what it was worth

| change | median km | verdict |
|---|---|---|
| light-response calibration (window + pooling + slope) | 473 -> 208 | the whole game |
| leave-one-out residence prior, weight 0.25 | 208 -> 180 | promising, p = 0.375 |
| truth-free search domain | no change | bit-identical; do it for defensibility |
| sensor fusion (bathymetry) | no change | buys SPEED (-31%), not accuracy |
| latent state-switching | no change | data prefer it; coverage worsens |
| dive effort as a movement covariate | untested on position | 1.33x on sigma; too small to matter |
| per-tag `lambda` by marginal likelihood | worse | falsified 0/10; default vindicated |
| composite-likelihood weight (`overcount`) | worse | falsified; recommend removal |

### Ideas that were tested and should NOT be retried without new information

Recorded because each looked sound enough to cost a run, and the next person
will have the same idea:

- **Bootstrapping the response geometry from a first-pass track.** Diverges,
  321 -> 391 -> 755 km, even confined to the calibration window.
- **Envelope-based light scale** (floor/amp instead of q05/rng). 572-894 km.
- **Whole-record intensity scale.** North bias, 486 km.
- **A composite-likelihood weight** to widen the posterior. Moves the point
  estimate; the miscalibration is a heavy tail, not under-dispersion.
- **Latent state-switching** for better coverage. Data prefer it emphatically
  (+431 log units, 10/10) but it TIGHTENS an already-narrow posterior.
- **Per-tag lambda by marginal likelihood.** No interior optimum, 0/10.
- **A hand-built position-free noise index.** r ~ 0.54, p = 0.11 at n = 10.

The pattern across all of them: every attempt to fix the OBSERVATION model's
calibration using only the observations has failed, while the one fix that
worked (the calibration window and pooling) used external structure -- a known
release site, the depth channel, and other tags of the same model.

Two traps this round, both worth remembering:

1. **`Select-Object -Last 5` in the launcher truncated a five-hour run's stdout.**
   The results survived only because the script also called `saveRDS`. Always
   persist results to disk, never to the terminal.
2. **Scripts that write only at the end lose everything to an interrupted run.**
   `diag_domain.R` ran through a laptop suspend and nearly lost 160 CPU-minutes.
   Append per-fit instead.
