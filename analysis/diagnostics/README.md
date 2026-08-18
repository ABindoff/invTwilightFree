# Diagnostics: the northern elephant seal calibration and latitude-bias investigation

The code and the reasoning behind a long investigation into a systematic latitude
bias in the grid HMM. Kept in the repository deliberately: the negative results cost
far more to find than the positive ones, and every one of them is easy to re-derive
by accident. If you are about to try something on the latitude bias, **read the
"Do not re-run" list below first.**

Excluded from the package build by `^analysis$` in `.Rbuildignore`, so none of this
reaches a CRAN tarball.

## The data is NOT here

These scripts need `data/nes_untracked/` and the cached archives under
`analysis/cache/nes/`, neither of which is in the repository. The northern elephant
seal records are the property of the TOPP programme and are not publicly
shareable; derived per-knot position tables are the same data in another form and
are excluded for the same reason.

So **the scripts will not run on a fresh clone.** They are here as a record of
method and reasoning, and as a starting point if you have equivalent data of your
own. Everything they concluded is written down — in this file, in
`LOG_calibration_window.md`, in `SESSION_STATE.md`, and above all in
`../../notes/latitude_bias_investigation.md`.

Run from the **package root**, not from this directory. `nes_common.R` carries the
crosswalk, the Argos reader, the departure detector and the scoring; nearly every
other script sources it.

## Start here

| document | what it is |
|---|---|
| `../../notes/latitude_bias_investigation.md` **§0** | **the current state.** Read this first; the rest of that note is layered with superseded material and marked as such |
| `LOG_calibration_window.md` | the 2026-08-06 calibration-window investigation, script by script, with what each one settled |
| `SESSION_STATE.md` | the long-form narrative record (~1900 lines), in rough chronological order |
| `EVIDENCE_DOSSIER.md` | a self-contained brief written for independent reviewers; the tightest summary of the evidence, though its framing of the mechanism has since been superseded by §0 |

## THE TRAP THAT COST A DAY

`build_light_battery.R` produces two light levels. **They are not equally valid.**

- **`perfect`** — noiseless: the clear-sky expectation itself. This is **NOT a valid
  null test.** Noiseless light is the *expectation* of the emission, not a draw from
  it, and evaluating a likelihood at `E[y]` is a different object from evaluating it
  at a sample. It manufactures a latitude tilt that is not there: measured per-knot
  tilt is **-0.113 nats/deg on the perfect arm and -0.0024 on the noisy arm**.
- **`noisy`** — clear-sky plus a draw from the engine's own spike-and-slab mixture,
  KS-verified against `LightMix`. The emission is correct by construction while the
  track stays a real Argos track. **This is the arm to use.**

An entire day was spent characterising a mechanism on the `perfect` arm, deriving a
correction for it, implementing that correction in Rust, and confirming it against
three pre-registered predictions — all of it measuring an artefact of the test. Every
"fix" that helped the perfect arm made the noisy arm worse.

**If you use the battery, use the noisy arm, and run any result past the other arm
before believing it.**

## Do not re-run: measured, and neutral or worse

| approach | result |
|---|---|
| Censoring the emission instead of truncating | mean tilt improves but the seasonal swing gets **20x worse** on the noisy arm |
| Symmetric spike arms (`shade_ratio = 1`) | bias worsens (-1.017 -> -1.291); the sweep is monotone and ratio 4 only trades one tilt for another |
| Ito drift correction for spherical movement | exact (nulls the one-step drift to 5 dp) and **inert** (posterior mean moves 9e-5 deg). Shipped as `drift_correction`, default FALSE |
| Tightening the movement prior to measured movement | D=48 gives bias **-2.95**, coverage 0.27, 409 km. The prior's excess width is doing real work |
| Per-axis anisotropic prior | no better than isotropic on bias; `aniso77_148` is slightly better *calibrated* (coverage 0.952 vs 0.996) at equal rmse |
| Posterior mean / median / marginal mode instead of the joint mode | the mean is **more** biased; the joint mode is the best available readout |
| Unclamping the expected curve | seasonal swing unchanged |
| Knot-grid phase (0/3/6/9 h) | **invariant** — bias moves by 0.01 deg. Binning alignment is exonerated |
| Composite-likelihood weight (`overcount`) | reverted; see `diag_overcount.R` |

Also excluded with evidence: Argos error (cancels by construction in the battery),
the solar zenith algorithm (agrees with an independent NOAA/Meeus implementation to
0.003 deg rms), grid quantisation (identical at 1.0 and 0.5 deg cells), and the
cell-area correction — which is not a bias source but is **preventing a 5-degree
poleward error**; removing it is catastrophic on the noisy arm.

## The live conclusions, and the scripts that carry them

| script | question | answer |
|---|---|---|
| `diag_noisy_arm_fits.R` | how does the engine do when the emission is correct? | **bias -0.51 deg, 189 km, coverage 0.996.** The coverage failure does not reproduce at all. So ~1 deg of the real-data bias and essentially all of the coverage problem is emission fidelity to real light |
| `diag_fisher_sloppy.R` | is the bias a gauge freedom? | Yes. Condition number **3.8e4**; latitude is essentially THE soft direction (0.881 in the softest eigenvector, 0.005 in the stiffest), trading against `z50` at 1.4-4.4 deg/deg with the exchange rate **changing sign across the year**. This is why axis-aligned sweeps returned nulls |
| `diag_gauge_breaking.R` | does anything identify the calibration? | Yes. Gauge inflation **1.81** (5-day window) -> **1.08** (240-day) -> **1.01** (6 tags pooled), `se(z50)` falling as 1/sqrt(N) to 0.045 deg. The seasonal contrast over a full deployment identifies it |
| `diag_profile_recovery.R` | can `log_z` recover the calibration? | **Yes, exactly.** Pooled evidence maximiser lands on the truth (+0.00) while per-tag maximisers scatter sd 1.03 and one runs to the boundary. Accuracy optimal at the same point (189 km, vs 480 at -2 deg) |
| `diag_area_diffusion.R` | is the cell-area correction the problem? | No, the opposite — it prevents a +5 deg poleward error on the noisy arm |
| `diag_zenith_reference.R` | is the solar zenith right? | Yes, to 0.003 deg rms against an independent NOAA/Meeus implementation. **No atmospheric refraction anywhere**: +0.22 deg in the twilight band, absorbed by the fitted `z50` |
| `diag_aniso_and_filter.R` | does the blue filter distort the twilight curve? | Not measurably. Seasonal swing at fixed zenith is 0.014 of range in the twilight band, so the filter acts as a stable function of zenith and empirical response fitting absorbs it |
| `diag_axis_scales.R` | what is the animals' real movement scale? | 12 h step: NS 17.2 km, EW 22.4. Both axes directed; NS saturates near 78 km/sqrt(day), EW keeps growing to 192. The shipped isotropic 110 is too wide in latitude and too narrow in longitude |
| `diag_mode_vs_mean.R` | is the bias a readout artefact? | No. The mean is more biased than the mode; the ordering follows the measured skewness |
| `diag_filtered_vs_smoothed.R` | does the bias appear in the forward filter? | Yes, same shape, no growth in knot index. Needs `grid_posterior(fit, filtered = TRUE)` |

`LOG_calibration_window.md` carries the same table for the earlier calibration-window
work, which established the recipe still in use: **pooled haul-out geometry with a
per-tag 15-day intensity scale**, and the tangent tied to `max_light`.

## The outstanding item

**Replace haul-out calibration with pooled evidence profiling.** `fit_light_response()`
estimates `z50` on a 15-day haul-out window spanning under a degree of declination —
where the Fisher analysis says the parameter is worst identified — and freezes it for
a deployment spanning 47 degrees. The identifying information is in the data and is
discarded by construction. Profiling `log_z` (already returned by the engine) over
`z50`, summed across tags, is validated above and needs no engine change.

Calibration error costs **accuracy, not bias**: ~125 km per degree of `z50` against
only 0.07 deg of latitude. Given a recorded 5.6 deg spread in haul-out `z50` estimates
across tags, that is a large accuracy cost currently being paid silently.

## Naming convention

`diag_*` diagnostics · `fit_*` fitting runs · `anal_*` analysis of fit output ·
`show_*` reporting · `sim_*` simulation · `ingest_*`/`probe_*`/`inventory_*`/`fix_*`
data handling · `nes_common.R` shared helpers · `build_light_battery.R` the synthetic
apparatus.

Scripts written before 2026-08-10 are documented in `LOG_calibration_window.md` and
`SESSION_STATE.md`; several were superseded by later work and a few were harness
failures kept because the failure was instructive. Where this README and an older
document disagree, this README and `notes/latitude_bias_investigation.md` §0 win.

## Response family: settled (2026-08-18)

`fit_family_control.R` -- pipeline recipe verbatim, three response arms, real light,
6 tags, scored with `argos_at`.

| arm | median km | mean bias | coverage |
|---|---|---|---|
| tangent (shipped) | **257** | -1.466 | **0.593** |
| darkness regime | 317 | +0.76 | 0.338 |
| gompertz (as lookup table) | 839 | **-0.076** | 0.077 |
| logistic | 624 | -3.825 | 0.269 |

**The trade is monotone.** Every step that corrects the response's twilight level and
floor reduces bias and costs distance and coverage: bias 1.47 -> 0.76 -> 0.08,
distance 257 -> 317 -> 839. The darkness regime is the EFFICIENT point -- half the
bias reduction for 60 km, where gompertz spends a further 520 km for the rest.

**Gompertz is NOT a shipping candidate**, despite fitting the clear-sky envelope 6.8x
better than the tangent in the twilight band at no gauge cost (1.01 vs 1.39). The
envelope prediction did not survive the fit-level test. The mechanism is visible: the
tangent descends linearly to zero at z ~ 106 and so still discriminates latitude at
high zenith, while gompertz asymptotes to its (physically real, measured ~0.10 of
range) floor and goes flat there. **The tangent's wrongness is load-bearing.** The
floor belongs in the LIKELIHOOD's darkness regime, not in the response.

**Nothing tested fixes coverage** (0.59 / 0.34 / 0.08, all far below 0.95, and worse
as bias improves). That is emission SPREAD, not response shape -- the spike is 2.4-4x
tighter than the assumed scale and the slab is ~5x too heavy (`diag_real_residual.R`).

### Two traps this run exposed

1. **`all29_tags.csv` is the PRE-area-correction baseline** ("242 was two bugs
   agreeing" -- SESSION_STATE). Comparing a faithful harness against it showed ratios
   up to 1.58 and looked like a broken harness. Use `rescore29_results.csv` arm
   `tangent_areaON` as the corrected-kernel reference.
2. **TWO DEFINITIONS OF GROUND TRUTH are in use.** `argos_at()` (nes_common.R, used by
   `fit_all29.R`) returns NA when the bracketing Argos gap exceeds 24 h;
   `truth_at()` (defined locally in `fit_rescore29.R`, `fit_lambda_stage2.R`,
   `fit_dark_sweep.R`, `fit_area_arm.R`) is `approx(rule = 2)` and never returns NA,
   interpolating across arbitrary gaps and extrapolating past both track ends. They
   differ on **4.5% of knots** -- and those are the worst-supported ones. Numbers from
   different scripts are therefore NOT directly comparable. `argos_at` is the
   defensible one; scoring against truth interpolated across a >24 h gap adds the
   interpolator's error to the model's, preferentially where Argos coverage is poor.
