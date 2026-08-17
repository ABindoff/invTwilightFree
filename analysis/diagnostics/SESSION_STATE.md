# Session state — resume point

Written 2026-08-10. Everything below is on disk; nothing here depends on a live
session. If a restart has wiped the R library, reinstall first: terra, sf,
rnaturalearth, ncdf4, xml2, testthat, then `R CMD INSTALL . --no-multiarch`
(note `R` is aliased to `Invoke-History` in this PowerShell, so use the full path
to `R.exe`).

---

## Where the work stands

**Committed** (branch `fix/light-response-autodetect`, 6 commits):

| commit | what |
|---|---|
| 177d934 | `overcount` composite weight — negative result, later reverted |
| 2f28746 | light-response calibration: `scale_light`, `pool_light_responses()`, tangent tied to `max_light`, new test file |
| 2bc32e0 | analysis rewrite, calibration-window comparison section, render path fix, cache-key versioning |
| b137008 | clock acceptance rule (offset AND drift, bound 20 min) |
| 0c6ab2e | clock off by default; `hier_light_forced` batch demonstrating the cost |
| 8433e48 | species-range search domain (bit-identical to Argos-derived) |
| 13fb074 | revert of `overcount` |

**Uncommitted in the working tree**: anisotropic diffusion (`diffusion_lon`) in
`src/rust/src/lib.rs` (both smoother passes, shared `log_move_kernel`),
`R/TwilightFreeGrid.R`, `R/extendr-wrappers.R`, and
`tests/testthat/test-anisotropic.R`. Built and passing (168 tests, 0 failures).
**Not adopted** — see the negative result below. Decide whether to keep it as a
documented, defaulted-off argument or revert it as `overcount` was.

**Manuscript**: `manuscripts/supplementary_calibration.md` (~4,200 words, S1–S7),
plus `gelman2006`, `bindoff2026flat`, `bindoff2026rbloo` added to the .bib with a
header note on updating the two preprints on acceptance.

---

## The data

29 double-tagged deployments, three seasons, after ingesting the F. Vilches
delivery (`fvilches/`, gitignored and .Rbuildignore'd).

- 2021 (10): `analysis/cache/nes/archives_v3_30min.rds`, from `data/nes_untracked`
- 2022 (9) + 2023 (10): `analysis/cache/nes/fvilches_archives_v1_30min.rds`
  and `fvilches_argos_v1.rds`; manifest in `fvilches_manifest.csv`
- 7 further Argos-only deployments exist (no light): 2021024, 2021029, 2022047,
  2022049, 2022050, 2023031, 2023033. Useless for geolocation, ideal for an
  occupancy prior with zero circularity.

**Four traps in the new delivery, all fixed in `fix_fvilches_argos.R`:**
1. PTTs are redeployed across seasons, so each RawArgos file holds every fix that
   transmitter ever produced — three animals were placed in Scotland, 8172 km
   out, and two given 613-day trips. Clip to each deployment's own archive span.
2. Three date formats (`17-May-2023`, `24-May-23`, and a mis-parse that looked
   like `yy-mm-dd`). Now tries five formats and fails loudly if none parses.
3. A class-0 junk fix at the head of a record poisons the speed filter, which
   anchors on the first fix. Now anchors on the first class 1–3 fix.
4. The `18A` tag family names the temperature column `External Temperature`, not
   `External Temp`; an exact header match silently dropped 4 deployments.

Tag serials are reused across seasons: 11 serials cover 25 deployments, three
appear in all three years. **Deployments are nested within physical tags** — this
is the untested opportunity (below).

---

## Headline results, 29 deployments (`all29_tags.csv`, `all29_knots.csv`)

| season | n | median km | rmse_lat | cover_lat |
|---|---|---|---|---|
| 2021 | 10 | 210 | 2.25 | 0.77 |
| 2022 | 9 | 242 | 2.48 | 0.73 |
| 2023 | 10 | 279 | 2.86 | 0.68 |
| all | 29 | 244 | 2.53 | 0.72 |

2021 reproduces the earlier analysis (210 vs 208–209 km), so the rewritten
pipeline is consistent. Own haul-out geometry 233 km vs borrowed 268 km
(p = 0.216). **The new seasons are genuinely harder and it is not explained by
borrowed calibration** — 2023 animals with their own haul-out still sit at 258 km
against 2021's 198. Unresolved: could be the different Argos processing
(RawArgos vs Locations), the different endpoint derivation, or real inter-annual
variation. Separate these before anything goes in the paper.

---

## THE CENTRAL FINDING

Latitude error is **not information-limited**. Across 315 windows:

| | Spearman ρ with the even-harmonic content E2 |
|---|---|
| model's **reported** latitude sd | **−0.926** |
| **actual** latitude error | −0.148 |

Mean absolute latitude error by |declination| band: 2.02, 2.69, 2.49, 2.04,
2.00° — essentially **flat across a 5.7-fold change in available day-length
information**. The model already widens its intervals correctly when the light
goes uninformative; the error simply does not follow. The residual is a
bias term of roughly constant size that more signal does not reduce and nothing
in the light identifies.

This explains every failure below and should be the Discussion's central claim.

---

## Everything tried, and what it was worth

| approach | effect | verdict |
|---|---|---|
| light-response calibration (window + pooling + tangent) | 473 → 208 km | **adopted, the whole game** |
| decimation by maximum, timestamped at the retained sample | 3650 → 852 km, then a 4° longitude bias removed | **adopted** |
| hard bathymetric constraint | no accuracy change, −31% runtime | **adopted for speed** |
| species-range search domain | bit-identical | adopted for defensibility |
| clock correction | +863 km, 0/10 improved | **excluded** |
| SST fusion | degrades latitude | **excluded** |
| leave-one-out occupancy prior | 208 → 180 km, p = 0.375 | unproven |
| posterior mean vs mode | −5.5 km, p = 0.176 (was p=0.055 at n=8) | unproven |
| latent state-switching | no change; coverage 0.77 → 0.65 | rejected |
| anisotropic movement prior | coverage balances, accuracy −12%, log_z worse | rejected |
| per-tag λ by marginal likelihood | 209 → 283 km, 0/10 interior optima | rejected |
| composite-likelihood weight | point estimate moves 5.7° | **removed** |
| dive effort as movement covariate | 1.33× on σ | too small |
| response drift over deployment | real (1.82°) but r = −0.07 with bias | not the mechanism |
| depth attenuation | 96% of windows reach <10 m | nothing to model |
| ALAN in 2021032 | darkest tag at night | excluded |
| declination sign-flip | 7/9 tags swing the wrong way | **falsified twice** |
| harmonic interval scaling | model already tracks E2 at ρ=−0.93 | **inverted my proposal** |

---

## THE ACTIONABLE BUG (2026-08-10) — `calibration_logistic` is on the wrong scale

**This is the highest-value item and it is a one-line fix.**

`fit_light_response()` returns `calibration_logistic = c(floor, amp, z50, scale)`
where `floor` is the RAW dark level. The engine's length-4 branch computes
`floor + amp/(1 + exp((z - z50)/s))` and compares it against
`pmax(0, light - baseline)`. The observations are baselined; the floor is not.

Measured mismatch across the 2021 ten: **mean 18.2 units, 11.4% of the response
range**, up to 25.2 on individual tags.

So the logistic expectation sits ~11% too high at every zenith. That is almost
certainly why the logistic tested at 833 km against the tangent's 354 and was
rejected in favour of the tangent: it was never a fair comparison. The tangent
won because its floor is 0 by construction — wrong, but consistently wrong.

FIX: express `calibration_logistic`'s floor relative to `baseline`, then re-run
logistic vs tangent under the corrected calibration.

Note: 2021032, the tag that has defeated every diagnostic, has a mismatch of
only 1.6 units against the panel's 18.2.

## WHY IT MATTERS — the measured response vs what the package can express

`diag_residuals.R` and `diag_empirical_response.R`. Observed light against TRUE
Argos zenith, upper envelope per 1-degree bin, all 29 deployments:

| zenith | measured | tangent (what the engine gets) | gap |
|---|---|---|---|
| 40-80 | 143-157 | 149-160 | -3 to -13 (shading, as assumed) |
| 85 | 134 | 120 | +14 |
| **92** | — | — | **+31.7 (peak, 20% of range)** |
| 95 | 90 | 62 | +28 |
| **>110 (night)** | **39-44** | **0** | **+41** |

TWO distinct misspecifications:
1. **A non-zero floor** (~40 = 25% of range) that the tangent cannot express at
   all, but a correctly-scaled logistic CAN. Fix the bug above.
2. **A twilight bump** peaking at zenith 92, which NO monotone sigmoid can make.
   The Mk9 carries an optical filter maximising sensitivity around twilight (for
   template fitting), and twilight light is blue-shifted by Chappuis absorption,
   so a blue-weighted channel reads relatively more there. The 18A family shows
   the same bump (+20.9 vs Mk9 +22.7), so it is not one instrument's quirk.

Residual by zenith band, showing the sign flip that drives the panel bias:
day -22/-34/-14 (shading, 1-2% brighter than model), then **86-90: +3.7 with 77%
BRIGHTER than the model; 90-94: +9.9 with 88% brighter**, then night +12.9.

The spike's upper arm is 2*lambda — observations ABOVE the curve are penalised
twice as hard — so through twilight, where latitude is read, the model is
penalising the large majority of observations at the harsh rate. It relieves that
by moving toward the subsolar point: north in northern summer, south in winter.
**That is the panel-norm swing (7/9 tags).**

Twilight residual does NOT predict per-tag bias (r = -0.127, p = 0.51), so this
explains the COMMON bias, not the 5x per-tag spread.

## RUNNING when this was written

`sim_perfect_light.R` — 87 fits, three arms, results append to
`sim_perfect_results.csv` / `sim_perfect_knots.csv`:
- **A_perfect**: light generated from the fitted response at TRUE positions, no
  noise, fitted with the same calibration. The model is EXACTLY correct for the
  data, so any residual bias is the ESTIMATOR (grid discretisation, movement
  prior, endpoint bridge, or the non-linearity of timing -> latitude near the
  equinox where dphi/dH -> infinity).
- **B_noise**: A plus the model's own spike-and-slab draws.
- **C_observed_shading**: A plus each tag's OWN observed-minus-theoretical
  residuals, resampled within zenith bins.

Real light for comparison: 244 km, rmse_lat 2.53, cover_lat 0.72.
If A is unbiased and C is not, C's residuals are the target — and they are now
characterised above.

## ARM A RESULT (16 of 29 tags) — the estimator is exonerated

Perfect light, model exactly correct for the data:

| | arm A (perfect) | real light, same 16 tags |
|---|---|---|
| median km | 52 | 235 |
| rmse_lat | 0.65 | 2.47 |
| cover_lat | **1.000** | 0.71 |
| reported sd / actual sd | **1.97** | 0.58 |

The estimator accounts for **7% of the latitude error variance**. More
decisively, the coverage failure **inverts**: with correct data the intervals
are twice too WIDE. So the interval machinery is not systematically narrow, and
the real-data under-coverage is entirely extra error the model does not know
about — i.e. observational misspecification, which is where the two
misspecifications above live.

There IS a bias: −0.334 deg, negative in 16/16 tags (p = 1.8e-10). It is
**shrinkage toward the anchored endpoints**, not a defect:

- by true latitude: −0.03 (35-40), −0.37 (40-45), −0.44 (45-50), −0.46 (50-55).
  Both endpoints are at Ano Nuevo ~37 N, so the pull grows with distance north.
- by |declination|: **−0.93 at equinox**, −0.67, −0.36, −0.23, −0.15 at solstice.
  A 6-fold swing: where the light stops identifying latitude, the prior takes
  over. This is the dphi/dH effect, and it is the ONLY thing found all session
  whose size tracks the available information.

Fitted latitudes land on cell centres (fractional part 0.500, 97% in the middle
two bins) — the noiseless posterior is near-degenerate, so this is NOT a
half-cell indexing offset. Discretisation is not the mechanism.

Per-tag simulated bias does not predict per-tag real bias (r = +0.098, p = 0.72),
so the estimator explains none of the 5x per-tag spread either.

**Conclusion: stop looking at the estimator. Arm C is the experiment that
matters.**

Arm A is now complete at 29/29 and every number held: 53 km, bias −0.345,
negative in **29/29** (p = 2.2e-19), cover_lat 0.999. The declination gradient
survived the full panel: −0.965, −0.720, −0.348, −0.239, −0.155.

## ARM B RESULT (18 of 29) — lambda is 3x too loose

| | A perfect | **B + model's own noise** | real light |
|---|---|---|---|
| median km | 52 | **299** | 242 |
| rmse_lat | 0.65 | **3.62** | 2.52 |
| bias_lat | −0.33 | **−2.43** | −0.12 |
| cover_lat | 1.00 | **0.73** | 0.70 |

**The model's own error model is worse than reality.** Data drawn from the
assumed spike-and-slab, fitted with the same likelihood that generated it,
lands 23% further out than the real light does. That can only mean the assumed
noise is far larger than the tag's.

Measured against 200k residuals at true positions, in units of the response
range:

| | mean abs residual / range |
|---|---|
| assumed by `lambda = 1/(0.5*max_light)` | 0.375 |
| **actual, all zeniths** | **0.121** (3.1x too loose) |
| **actual, 86-94 deg — where latitude is read** | **0.086-0.103 (~4x)** |

So the likelihood is too flat and the fit is discarding information it has.
`likelihood_params[1]` is a one-number change, testable on the same 29
deployments. This is a THIRD actionable defect alongside the floor and the bump.

CAVEAT before acting: a too-loose lambda widens the posterior, which on its own
would OVER-cover, yet the real data under-covers at 0.70. Tightening lambda will
improve accuracy but will make coverage worse unless the response
misspecification is fixed first. **Fix the floor and the bump, then retune
lambda** — not the other way round.

Arm B's own −2.43 bias is the asymmetric spike doing what it is designed to do
at an absurd scale: the lower arm has twice the mean of the upper, so simulated
light is systematically dimmer, days look shorter, latitude runs south. At
lambda 3x too loose that asymmetry is worth 2.4 degrees.

## ARM C RESULT — THE DIAGNOSIS CLOSES (COMPLETE, all 87 fits, 29/29/29)

Perfect geometry, plus each tag's OWN observed-minus-theoretical residuals
resampled within 5-degree zenith bins.

| | median km | rmse_lat | cover_lat | rmse_lon |
|---|---|---|---|---|
| A perfect | 53 | 0.67 | 1.00 | 0.39 |
| B model's own noise | 314 | 3.58 | 0.73 | 1.54 |
| **C real shading** | **218** | **2.30** | **0.79** | **0.71** |
| **REAL LIGHT** | **244** | **2.53** | **0.72** | **0.91** |

C reproduces **90% of the km error and 91% of rmse_lat**, and the bias lands at
−0.304 against the real −0.166.

The result STRENGTHENED with n, which is the opposite of this session's pattern
(mean-vs-mode decayed p=0.055 -> 0.176; anisotropy's log_z gain went +600 -> −31).
Tracked across the run: km reproduction 87% -> 88% -> 90%, bias rho 0.777 ->
0.775 -> 0.801. Nothing here is a small-sample artefact.

**And it predicts WHICH tags are bad** — the 5x per-tag spread that nothing all
session has explained:

| | Spearman rho, C vs real (n = 29) |
|---|---|
| median km | **+0.922** (p < 0.001) |
| rmse_lat | **+0.925** (p < 0.001) |
| **bias_lat** | **+0.801** (p < 0.001) |

Transplanting a tag's residuals onto perfect geometry reproduces that tag's
real-world failure. The latitude error is observational, it lives in the light
residual, and it is tag-specific. **This is the paper's mechanism result.**

### The part that changes the scattering plan

C resamples residuals **within zenith bins**, which destroys all time structure
while preserving the zenith-conditional marginal. So C *is* the shuffled control
for burstiness. It still reproduces 90% of the error. Therefore:

- the **zenith-conditional residual distribution** drives accuracy and the
  per-tag spread,
- the **temporal burst structure** is worth at most the remaining ~13%.

But look at coverage: C sits at 0.79 while the real light sits at 0.72, and C is
the only arm that fails to reproduce the under-coverage. Correlated residuals
inflate the effective sample size, which is exactly how temporal structure would
show up. **So the burst structure shows up in COVERAGE, not in accuracy.**

That gives scattering (or any correlated-error model, or simply a
variance-inflation factor) a precise and much narrower job: fix the intervals,
not the point estimate. It is a downgrade of the earlier enthusiasm and should
be recorded as one.

Longitude: C 0.71 vs real 0.91, so C reproduces only 78% of the longitude error.
Something in longitude is still unaccounted for. Timing is the obvious suspect,
and the clock estimator is already excluded, so this needs its own look. Note
this is the ONE place where accuracy is under-reproduced, and it is the same
place time structure would act.

## IS THE RESIDUAL A DEFORMATION PROCESS? (`diag_scattering.R`, 27 deployments)

Asked because a scattering-transform likelihood is only worth building if the
residual has time structure. Hand-rolled a-trous first- and second-order
scattering, against a within-tag shuffled control (same marginal, no time
structure). Per-tag coefficients in `scattering_probe.csv`.

**Variance split:** 52% of the residual is deterministic in zenith (the floor
and the bump, removable with a better response), 48% stochastic. So the two
parametric fixes have a hard ceiling of halving the residual variance.

| scale | S1 real | S1 white | S2 real | S2 white | S2 ratio |
|---|---|---|---|---|---|
| 0.5 h | 0.373 | 0.499 | 0.199 | 0.065 | **3.0** |
| 1 h | 0.187 | 0.250 | 0.161 | 0.050 | **3.2** |
| 2 h | 0.078 | 0.126 | 0.106 | 0.039 | **2.7** |
| 4 h | 0.054 | 0.063 | 0.078 | 0.039 | 2.0 |
| 8 h | 0.057 | 0.032 | 0.103 | 0.055 | 1.9 |
| 16 h | 0.061 | 0.018 | 0.160 | 0.144 | 1.1 |

The control's S1 halves exactly per octave, the white-noise signature, which
validates the transform. The real S1 flattens then rises, carrying 3.4x the
white energy at 16 h. S2 excess of ~3x at 0.5-2 h is genuine amplitude
modulation, i.e. bursts. **The residual is a deformation process**, so
scattering is the right family of tool, but it is third in line behind the
floor and the bump.

**Unexpected: a diel term.** ACF of the stochastic residual runs 0.255, 0.288,
0.329, 0.266 over the first two hours (peaks near 1.5 h, NOT AR(1)), decays to
0.05 by 12 h, then returns to **0.258 at 24 h**. Zenith-detrending cannot remove
this because zenith is ambiguous between morning and evening, so it is a real
dawn-versus-dusk asymmetry. Diel dive behaviour is the obvious candidate. A
time-of-day covariate may be much cheaper than any transform. **Untested.**

Do not quote the script's own "decorrelation time 0.5 h": the 1/e rule is
inappropriate when the ACF never starts near 1. The honest description is a
large independent component plus a ~30%-amplitude correlated component
persisting one to four hours.

## THE SCALE BUG IS FIXED (2026-08-11)

`R/light_response.R`: `calibration_logistic[1]` is now `max(floor_l - q05, 0)`
in `fit_light_response()` and `max(r$floor - r$baseline, 0)` in
`pool_light_responses()`. Clamped at zero because the observations are and the
spike normaliser only supports [0, max_light]. Roxygen updated; five new
assertions in `test-light_response.R` (38 pass), including a pooling check so
the raw scale cannot creep back in.

The old value was **worse than the 18.2 estimate**: the engine was being handed
a mean floor of **74.1 units, 46% of the response range**. The corrected value
is 22.5 units, 14.0%. That fully accounts for the logistic's historical 833 km
against the tangent's 354, and means the four-parameter form has never actually
been tested. `fit_logistic_vs_tangent.R` is running, 58 fits.

### THE LOGISTIC IS DECISIVELY REJECTED, AND WHY IS THE INTERESTING PART

Paired on 25 of 29 tags:

| | median km | rmse_lat | **rmse_lon** | bias_lat | cover_lat |
|---|---|---|---|---|---|
| tangent | 246 | 2.58 | **0.98** | −0.06 | 0.71 |
| logistic | **634** | **6.67** | **0.94** | **−2.21** | **0.29** |

Better on 3 of 25 tags, Wilcoxon paired p = 4.6e-06, bias negative on **25 of
25**. The scale fix was correct and necessary and did NOT rescue the
four-parameter form. The "unfair comparison" hypothesis is dead.

**Longitude is untouched** (0.94 against 0.98, ratio 0.96) while latitude
triples. The response shape damages latitude only, exactly as the phase/duty
decomposition says it should: longitude is the carrier phase and does not care
what the amplitude curve looks like.

### HYPOTHESIS: the tangent wins because it is FLAT, not because it fits

This is the day's synthesis and it ties to ADB's threshold rationale above.

The tangent is clamped, so it is constant at `max_light` below z50−2*scale
(~76 deg) and constant at 0 above z50+2*scale (~108). Flat means insensitive to
position, so those observations cancel from the likelihood ratio. **The tangent
is an implicit threshold method**: it uses only a ~33 degree window and ignores
the rest, without anyone having chosen to.

The logistic is never flat. Its daytime limb slopes about 0.5 units/deg (143 at
40 deg falling to 122 at 80), so the 73% of observations in the flat regions now
carry positional leverage. Those are exactly the observations that are noisiest
(daytime residual sd 0.169 of range against 0.094 at night) and where the model
is most wrong. Latitude collapses.

Note this beats the alternatives on the evidence: it is not curve fidelity,
because the tangent describes the measured response terribly (rmse 25.1, against
the table's 5.1) and still wins by 2.6x.

**FALSIFIABLE PREDICTION, make it before running the table arm.** The measured
table is also not flat: 160 at 40 deg, 154 at 50, 152 at 60, 149 at 70, and a
night floor that wobbles 29/38/32 at 110/120/130, which is noise. If the
flatness hypothesis is right, **the plain table will fail too**, and a
CLAMPED table — measured shape through the transition, forced flat outside it —
will beat both. That is the measured-response version of ADB's threshold, with
no tuning constant to defend.

Two arms discriminate: (a) clamped table vs plain table, (b) the LOGISTIC with
the likelihood restricted to an 80-105 window. If (b) rescues the logistic, the
mechanism is confirmed and it is about flatness, not about shape.

### THE TABLE ARMS: MY PREDICTION FAILED (complete, 29/29/29/29)

All four responses, same 29 deployments, everything else identical:

| arm | median km | rmse_lat | rmse_lon | bias_lat | cover_lat |
|---|---|---|---|---|---|
| **tangent** | **242** | **2.51** | 1.16 | **−0.14** | **0.73** |
| table_flat | 451 | 5.65 | 1.09 | +1.96 | 0.36 |
| table_plain | 481 | 5.73 | 1.12 | +1.71 | 0.36 |
| logistic | 615 | 6.56 | 1.14 | −2.24 | 0.31 |

I predicted `table_flat` would beat all three. It does not: +210 km against the
tangent, better on only 6 of 29, p = 7.7e-05. **Recorded as a failed
prediction.**

**The flatness effect is real but small.** `table_flat` beats `table_plain` on
**26 of 29** tags, p = 3e-06 — the direction was right and it is about as clean
a paired result as this data produces. But it is worth only 30 km of the 210 km
gap. Flatness is a real second-order effect, not the explanation.

**Longitude is identical across all four** (1.09-1.16). Four wildly different
response curves, latitude ranging over 2.5x, and longitude does not move. The
response acts on latitude alone, which is the phase/duty decomposition showing
up as cleanly as it ever will. This is a paper figure.

**The tangent is the only form whose bias is near zero** (−0.14 against +1.96,
+1.71, −2.24). Note the three alternatives are badly biased in BOTH directions,
so this is not one shared mechanism with a sign.

### WHERE THE SPIKE ACTUALLY FIRES (`diag_penalty_budget.R`, 322,693 obs)

ADB's question: how much light is being penalised that should not be penalised
at all? Answer: most of it.

The tangent's penalty budget, by zenith band:

| band | % of obs | expects | tag reads | % on upper arm | **% of budget** |
|---|---|---|---|---|---|
| day <80 | 41.5 | 159.4 | 131.3 | 1.7 | 0.4 |
| 80-86 | 5.5 | 131.0 | 116.6 | 18.4 | 1.3 |
| 86-94 | 7.1 | 90.6 | 97.4 | 82.3 | 13.1 |
| 94-100 | 5.4 | 50.4 | 51.4 | 58.0 | 5.3 |
| 100-106 | 5.4 | 15.8 | 15.5 | 39.8 | 4.7 |
| **night >106** | **35.2** | **0.0** | **12.7** | **61.1** | **75.2** |

**Three quarters of the spike's entire budget is spent at night, on a tag whose
dark reading is simply not zero.** Under a correctly-floored table the same
observations cost 0.027 nats instead of 0.316, a 12-fold drop, and the total
budget falls 91% (47,729 nats to 4,495). **69% of the tangent's whole penalty
budget is spurious night dark-current.** Per 12-hour knot the tangent charges
3.55 nats, 2.67 of it at night.

A Kuo-Mallick spike is meant to fire on ANOMALOUS light. It is firing on the
sensor's noise floor, on 35% of the record.

### AND YET THE TANGENT WINS. The mis-specification is LOAD-BEARING.

242 km against the correctly-floored table's 451. The spurious penalty is not
harmless, it is doing the work. Mechanism:

The penalty is not constant across candidate positions. Where the model predicts
night, `mu = 0` and the observation costs 0.316 nats; where it predicts twilight,
`mu > 0` and it costs less. **That differential IS the latitude discriminant.**
Expecting zero at night turns "should it be dark here?" into a sharp binary test
with a large lever arm (12.7 units times lam_hi), and day length is exactly what
encodes latitude. Give the response its true floor and the lever collapses:
night costs 0.027 nats wherever you stand, night stops carrying information, and
latitude degrades to 451 km.

So the engine's accuracy rests on an accident. The zero-floor error is
accidentally implementing a day-length likelihood, which is the right idea
arrived at by the wrong route, and it explains everything the last two days have
turned up: why the tangent beats better-fitting curves, why flatness helps only
30 km of 210, and why longitude never moves (phase does not care about darkness;
duty cycle is darkness).

CAVEAT: the lever-arm argument is inferred from residuals at TRUE positions plus
the table result. It has NOT been measured directly. Do that by computing how
`mu` and the penalty vary across candidate latitudes at fixed time, under both
responses. Cheap, and it would turn the mechanism from well-supported into shown.

**THE REFORMULATION THIS POINTS TO** — and it is close to ADB's original design.
Separate the two jobs the tangent is currently doing at once:
1. a correctly-floored response for the CONTINUOUS part, so shading is modelled
   honestly rather than being confounded with darkness;
2. an explicit Kuo-Mallick indicator on "light where darkness is expected",
   carrying the day-length information deliberately and with a stated lever arm,
   instead of smuggling it in through a floor that is wrong by 40 units.

That is a principled model rather than a lucky one, and it is defensible in a
way "our expected-light curve is wrong at night and that is why it works" is not.

### NEW HYPOTHESIS: `expected` is a location parameter, not the clear-sky curve

The tangent fits the measured response WORST (rmse 25.1 against the table's 5.1)
and wins by 2.5x. Flatness explains 30 km of 210. So something else is doing the
work, and the likelihood's asymmetry is the obvious suspect.

The spike is asymmetric: below `expected` the rate is `lam_lo`, above it
`lam_hi = 2*lam_lo`. Calibration currently targets the **95th-percentile
envelope**, so by construction ~95% of observations fall on the CHEAP shading
arm, where the log-density's sensitivity to `mu` is the smaller rate. The
tangent sits far below the envelope at night (0 against a measured 39-44), which
puts night observations on the STEEP arm where they discriminate.

So the curve the engine wants may not be the physical clear-sky maximum at all,
but whatever location makes the asymmetric density most informative.

**TEST, cheap and decisive:** build tables at several envelope quantiles
(`env_q` = 0.5, 0.75, 0.95) and score them. If a lower-quantile table beats the
95th, the hypothesis is confirmed and there is a real improvement available,
because the table can sit wherever we want while keeping the measured SHAPE.
This also reframes `env_q` from a nuisance argument into the main tuning surface.

### A real but MINOR bug found on the way (do not quote it as the cause)

The logistic branch of `expected_light()` was unclamped while the linear branch
was not. `floor + amp` is `top - baseline` where `top` is the envelope's
brightest bin, but `max_light` is `q95 - q05` over the whole record, and for a
diving animal most of the record is dark. So the logistic curve exceeds
`max_light` on **10 of 10** tags checked, by 3.4% on average, and
`spike_normaliser`'s upper mass term goes negative. Now clamped in both the Rust
and the R twin. At a 3.4% overshoot the total normaliser is still positive and
far from the 1e-12 floor, so this is tidiness, NOT the explanation for the
collapse. My first reading of it as a broken likelihood was wrong.

### Early result: the corrected logistic is still far worse (n = 1)

Tangent arm complete, 29 fits, and it reproduces `all29_tags.csv` once pooling
was fixed to be WITHIN FAMILY (my first version pooled across all 29 and gave
141 km against the recorded 149 on 2021023; family-wise gives 148).

First logistic fit, tag 2021023: **845 km, bias_lat −7.23, cover 0.17**, against
the same tag's tangent at 148 km. The historical logistic scored 833 km. So the
scale fix, while correct and necessary, has NOT rescued the four-parameter form.
The "unfair comparison" hypothesis looks wrong. Likely cause: even corrected,
the logistic's transition is far too wide (width_deg ~35) against a measured
collapse of ~15 degrees, and it attains neither 0 nor max_light, which the spike
normaliser is documented to dislike. n = 1, so let it finish.

NOTE ON BUILD STATE: the running job loaded the package at 07:44 and that build
HAS the scale fix and does NOT have the table work (verified). A later
`R CMD INSTALL` failed with a locked DLL and restored the previous version, so
the system library is the pre-table build. The new build with
`light_response_table()` is installed in the scratch `testlib`. Reinstall to the
system library once no job is running.

### A trap: the floor mismatch does NOT explain the season effect

Seductive, and wrong. By season the alignment is perfect and monotone:

| season | n | floor mismatch | median km | rmse_lat | cover |
|---|---|---|---|---|---|
| 2021 | 10 | 11.0% | 210 | 2.25 | 0.77 |
| 2022 | 9 | 13.9% | 242 | 2.48 | 0.73 |
| 2023 | 10 | 17.0% | 279 | 2.86 | 0.68 |

But at the TAG level there is nothing: rho = −0.027 (p = 0.89) against median km,
−0.035 against rmse_lat, +0.002 against coverage. Within season, −0.467, −0.017,
+0.030, none significant. Three season means agreeing in order happens by chance
one time in six.

This is an ecological correlation, the same shape as the dive-effort Simpson's
paradox earlier in the session. **Do not use that table.** The season effect
remains unexplained and confounded with Argos processing and endpoint
derivation.

## THE NON-PARAMETRIC RESPONSE IS IMPLEMENTED (2026-08-11)

**First, a correction to the earlier diagnosis.** The "twilight bump" is NOT a
local maximum. The measured envelope is monotone (largest rise +2.0 in 90 steps
for Mk9). The +31.7 was the gap to the TANGENT, not a rise in the response.
What the channel actually does is: a daytime plateau near 155, a collapse across
roughly 15 degrees, and a floor at 39-44. The tangent's transition is too wide
and mis-centred, and it has no floor. So "no monotone sigmoid can express it"
was wrong. The right statement is that no TWO-parameter clamped line can place a
floor and a narrow transition at once.

### `calibration` now accepts a lookup table

`c(-1, z_min, dz, y_0, ..., y_{n-1})`, linearly interpolated on a uniform zenith
grid, held constant beyond either end, clamped to [0, max_light]. A negative
leading element is the sentinel: unambiguous because the linear form's intercept
is `slope * zero_at` and the logistic floor is now clamped at zero.

- `src/rust/src/lib.rs`: the branch in `expected_light()`, plus
  `reject_table_calibration()`.
- `R/light_response.R`: `light_response_table()` (exported), the R twin in
  `.tf_expected_light()`, and `envelope$measured` — the envelope BEFORE the
  `cummin`, because the monotone constraint is wrong at night (the dark reading
  drifts up with zenith and `cummin` pins it to the running minimum, biasing the
  fitted floor low).
- `tests/testthat/test-light_response.R`: full suite 189 pass, 0 fail.

**LIMITATION, deliberate and enforced.** The block samplers (`run_block_track`,
`run_block_hier`, so `TwilightFreeHier`) carry the calibration in a fixed
`[f64; 4]` and cannot hold a table. They now panic with a clear message rather
than truncate, which would have left `[-1, z_min, dz, y_0]` to be read as a
logistic with a floor of −1. Also fixed a latent panic there: `cal_stride > 4`
indexed past the array. Grid and particle engines take the table fine.

### It recovers the measured shape, without truth

`diag_table_shape.R`, 20 Mk9 haul-out fits pooled, against the envelope measured
at true Argos positions:

| | rmse to measured envelope |
|---|---|
| **table** | **5.1** |
| tangent | 25.1 |
| logistic with floor 0 | 25.9 |

Twilight 84-100: table 4.7, tangent 23.6. Night >110: table 8.1, tangent 41.4.

CAUTION: this is a statement about the CURVE, not about position error, and the
two have come apart before — the tangent beat the better-fitting logistic
because fitting the curve is not the objective. `fit_table_arm.R` is written and
scores it against Argos. **Do not run it until `fit_logistic_vs_tangent.R`
finishes**: the DLL is locked while a job runs, and concurrent fits contend.

F18A has no haul-out geometry in this pipeline, so its 4 deployments fall back
to the deployment-start windows. Worth checking why `departure_time()` finds
nothing there.

## DESIGN RATIONALE — why the original TwilightFree thresholded (for the Discussion)

ADB's original method used a threshold and modelled shading as a distribution
around the twilights, rather than using the whole light curve. Today's evidence
supports that choice, and sharpens why.

**73% of observations sit where the fitted response is flat**: 41% below 80
degrees, 31% above 110, leaving only ~27% in the 80-110 band where expected
light varies with zenith. And the flat regions are NOISIER, not quieter.
Residual sd as a fraction of range: 0.121 below 60, **0.169 across 60-80**,
falling to 0.095 at 94-98 and 0.094 at night.

The mechanism is not merely low signal-to-noise. If expected light were flat
across every candidate position those observations would cancel from the
likelihood ratio and cost only time. But the tangent forces expected light to
ZERO at night while the tag reads 39-44, so every night observation carries a
+40 residual on the upper arm at 2*lambda, about one nat each. That is constant
only where the model predicts night at all candidate positions. Near the edges
of the domain it does not, so a position predicting more darkness is penalised
more than one predicting twilight. **The least informative part of the record
emits a systematic pull toward shorter nights, i.e. toward the summer pole** —
the direction of the panel bias, with 73% of the data behind it.

The table fix should reproduce the threshold's benefit endogenously: predict 40
at night, observe 40, and the residual is near zero AND near constant across
positions, so the night data neutralises itself instead of being excluded. That
is the better story for a paper, since a threshold is a tuning choice a reviewer
will ask you to defend and a measured response is not.

TESTABLE PREDICTION, two cheap arms on the same 29 deployments: with a corrected
floor, restricting the likelihood to a zenith window should buy almost nothing;
with the tangent it should buy a lot. That settles whether the continuous
likelihood earns its keep or has been getting away with it because the
calibration work masked the cost.

## PROCESS TRAP: install AFTER the last source edit, not before

The table arms were launched on 2026-08-11 and died 20 seconds in with
`unused argument (flat_outside = win)`, then sat idle overnight. Cause: I ran
`R CMD INSTALL`, THEN added the `flat_outside` argument, tested it with
`pkgload::load_all()` (which reads source, so it passed), and launched a job that
loads the INSTALLED package. Sixteen hours lost.

`load_all()` passing is not evidence that a background job will run. Before
launching anything long, assert the feature against the installed build:

    Rscript -e "library(invTwilightFree); stopifnot('flat_outside' %in% names(formals(light_response_table)))"

Relaunched 2026-08-12 07:30, verified first.

## DARKNESS-REGIME SWEEP, ROUND 1 (48 fits, 8 tags, complete)

Subset baselines: tangent 257 km / cover 0.68; table_flat 477 / 0.31.

| ratio | pslab | km | rmse_lat | cover | beats table_flat | beats tangent |
|---|---|---|---|---|---|---|
| 50 | 0.35 | **303** | 4.37 | 0.44 | **8/8** | 3/8 |
| 50 | 0.15 | 314 | 4.58 | 0.42 | 8/8 | 3/8 |
| 20 | 0.15 | 339 | 4.37 | 0.39 | 8/8 | 2/8 |
| 20 | 0.35 | 339 | 4.50 | 0.40 | 8/8 | 2/8 |
| 8 | 0.35 | 395 | 4.97 | 0.35 | 8/8 | 2/8 |
| 8 | 0.15 | 399 | 4.90 | 0.33 | 8/8 | 1/8 |

**The indicator works, and closes 79% of the gap.** Every cell beats the
corrected-floor response on 8 of 8 tags; the best is 303 against 477,
p = 0.0078. It does NOT yet beat the tangent (303 against 257, better on 3/8,
p = 0.55), so the accidental lever is still ahead overall.

**ROUND 2 COMPLETE (72 fits). The optimum is now bracketed and it SATURATES**:

| ratio | 8 | 20 | 50 | 100 | 200 | 400 |
|---|---|---|---|---|---|---|
| median km | 397 | 339 | 308 | **291** | 290 | 289 |
| cover_lat | 0.34 | 0.39 | 0.43 | 0.47 | 0.48 | 0.49 |

Flat from ratio 100 onward, exactly as the slab cap −log(pslab/max_light)
predicts: a sharper arm cannot cost more than the contamination floor, so it
saturates instead of running away. **Use ratio ~100; anything beyond is free but
pointless.** Best cell 289 km against table_flat's 477 and the tangent's 257,
which closes **85%** of the gap. Still 3/8 against the tangent (p = 0.64).
`prob_slab_dark` barely matters (0.15 against 0.35) EXCEPT at ratio 50, where
0.15 gives 259 km on the contaminated tag and 0.35 gives 174 — a sharp arm needs
a bigger slab, exactly as designed.

### The tag the mechanism was built for

2023037, which carries the most bright-when-dark contamination (8.7% against a
median of 1.8%):

| | median km | cover_lat |
|---|---|---|
| tangent | 461 | 0.26 |
| table_flat | 274 | 0.55 |
| **darkness regime (ratio 20)** | **154** | **0.79** |

**Three times better than the tangent, and coverage from 0.26 to 0.79.** This is
the strongest single-tag result of the session and it is on the tag predicted in
advance to benefit most. Whatever happens to the panel mean, this is the
existence proof that the mechanism is real and correctly identified.

But at ratio 100/200/400 this tag DEGRADES: 182, 200, 217 km. Its optimum is
ratio 20 while the panel's is 100+. **The two disagree, and the reason is
structural.** Once the arm is sharp enough, every excursion above the dark
envelope saturates at the slab cap, so a moderate genuine excursion and gross
contamination cost the same and the gradation between them is lost. A clean tag
does not care, because it has few excursions; a contaminated tag needs the
gradation. That argues for estimating the contamination rate PER TAG rather than
fixing one weight for the panel, which is what the particle filter already does
for `prob_slab` and what the fraction of night observations above the envelope
gives directly (median 1.8%, worst 8.7%).

### The +1.4 deg north bias is NOT the pooled floor (`diag_floor_bias.R`)

Mechanism tested: if the table hands a tag a dark envelope BELOW what that tag
measures, its ordinary night readings look like "not dark", the night looks
short, and the fit runs poleward. Correlation of (table floor − own floor)
against bias: rho = −0.209, p = 0.277. Right sign, not significant, n = 29.

The test is also weaker than it looks and should not be quoted as a clean null.
The comparison baseline is each tag's `dark_level`, which the `cummin` biases
LOW (that is why the table floor averages 32.1 against `dark_level`'s 22.5). So
this compared the table against a known-biased yardstick. **The bias remains
unexplained**, and it is now the dominant residual problem: it is the reason the
tangent still wins on the mean despite the indicator recovering 85% of the gap.

## THE MOVEMENT KERNEL WAS BIASED POLEWARD. FOUND, QUANTIFIED, FIXED.

**This is the session's most important result and it is a genuine package
defect, not a modelling choice.**

`run_grid_hmm` builds its transition weight as `exp(-dist^2/var2)` times a
CONSTANT, with no per-cell normaliser and no cell area. The kernel is a density
on the sphere evaluated on a grid uniform in DEGREES, so a cell's probability
needs the `cos(latitude)` Jacobian. Without it every cell counts equally, and
because high-latitude cells are physically smaller they collect weight they have
not earned.

Isolated with no data at all (`diag_kernel_bias.R`), on the 100x50 degree grid
the analysis uses with a 110 km step:

| destination lat | relative incoming mass | 1/cos(lat) |
|---|---|---|
| 22 N | 0.788 | 0.802 |
| 42 N | 1.000 | 1.000 |
| 50 N | 1.156 | 1.156 |
| 66 N | 1.828 | 1.827 |

**r = 0.9999.** Mass at 68 N is 2.32x that at 22 N. Between 38 N and 50 N it is
**0.204 nats PER STEP** poleward, about 98 nats over a 480-knot track unopposed.

### It explains what nothing else did

The pull between what the emission prefers (latitude profiles) and where the
track lands was **positive in all four responses** (+2.78, +1.29, +1.14, +0.53)
whichever way their emissions leaned. And it explains why tempering the light
drove the fitted latitude past the domain centre to 51.2 rather than stopping at
45: weakening the emission hands the answer to a prior that leans north.

### The fix, and the trap in it

`area_correction = TRUE` (default; `FALSE` reproduces archived results).
Adds `ln(cos(lat))` for the DESTINATION cell, plus the initial distribution,
which had the same defect.

**The forward pass's destination is `i`; the backward pass's destination is `j`,
because it runs from knot k to k+1.** Applying the area to `i` in both would be
silently wrong and no isotropic test would catch it. Second time this session
the two smoother passes have differed in exactly this way.

### Measured, against predictions made in advance

Validation: `areaOFF` reproduces the archive on **8 of 8 tags exactly**.

| response | bias before | bias after | shift | predicted |
|---|---|---|---|---|
| tangent | +0.053 | **−1.375** | **−1.427** | −1.14 |
| dark | +1.350 | **+0.350** | **−1.000** | −1.29 |

Bias moved south on **16 of 16 tags**. The predicted shifts came from the
profile accounting and landed within 25%.

Full 8-tag result (32 fits, complete):

| response | area | km | bias | rmse_lat | cover |
|---|---|---|---|---|---|
| tangent | OFF | **257** | +0.05 | 2.75 | 0.68 |
| tangent | ON | 289 | −1.37 | 3.52 | 0.64 |
| dark | OFF | 291 | +1.35 | 4.37 | 0.47 |
| dark | ON | **286** | **+0.35** | **3.20** | 0.54 |

### THE CONSEQUENCE: the tangent's advantage was two errors cancelling

With the kernel fixed the tangent's track bias (−1.375) now almost equals its
EMISSION bias (−1.28) — the accounting closes, which is the check that the
defect really was the prior. But its accuracy gets WORSE, 257 to 289 km, because
its near-zero bias was a southward emission error cancelling a northward prior
error. Remove one and the other is exposed.

The 34 km gap between tangent and dark closes to a **tie** once the kernel is
fixed: 289 against 286. Head to head with the corrected kernel, dark is better
on 4 of 8 tags, mean −4 km, Wilcoxon p = 1.00. Per-tag differences run −317 to
+199, so 8 tags cannot separate them. Dark carries less bias (mean |bias| 0.87
against 1.37, smaller on 5 of 8) and better rmse_lat (3.20 against 3.52) but
worse coverage (0.54 against 0.64). None of that is significant at n = 8.

**Honest statement: with the defect fixed, the tangent's advantage disappears
and neither response wins.** Every ranking on disk was produced with the
poleward kernel and is no longer valid.

## THE RE-BASELINED RESULT, n = 29, CORRECTED KERNEL (58 fits, complete)

| response | km | bias | mean abs bias | rmse_lat | rmse_lon | cover |
|---|---|---|---|---|---|---|
| tangent | 276 | **−1.543** | 1.54 | 3.45 | 1.15 | **0.66** |
| dark | **272** | **+0.458** | **0.73** | **3.04** | 1.25 | 0.57 |

archived tangent on the BIASED kernel: 242 km, bias −0.14, cover 0.73

### The headline accuracy gets WORSE, and that is the correct outcome

The tangent goes 242 to 276 km, better on only 4 of 29 tags, and its bias goes
−0.14 to −1.54. The old number was flattered by a southward emission error
cancelling a northward prior error. **276 km is the honest figure for this
dataset; 242 was two bugs agreeing.** Any manuscript number must come from the
corrected engine.

### Head to head on a level field

| comparison | result | p |
|---|---|---|
| median km | dark better on 16/29, mean −4 km | **1.00** |
| rmse_lat | dark better on 16/29 | 0.325 |
| **mean abs bias** | **dark smaller on 23/29** | **0.0002** |
| coverage | tangent better (0.66 vs 0.57) | — |

**Distance is a genuine tie.** The one decisive result is bias: the darkness
regime is half as biased and beats the tangent on 23 of 29 tags at p = 0.0002.
That is the first strongly significant response comparison of the session, and
it only became visible once the kernel stopped leaning.

So the choice is now explicit rather than accidental: **dark for unbiasedness,
tangent for coverage.** Neither is better outright, and the honest framing for a
paper is that the clamped-linear response carries a systematic −1.5 degree
latitude bias that was previously masked.

### By season

| resp | 2021 | 2022 | 2023 |
|---|---|---|---|
| dark | 335 km / cover 0.41 | **213 / 0.74** | 262 / 0.57 |
| tangent | 253 / 0.71 | 271 / 0.68 | 304 / 0.58 |

Dark is much better on the two fvilches seasons and much worse on 2021, while
the tangent is uniform. Unexplained and worth chasing: 2021 is the original
delivery with its own Argos processing and endpoint derivation, which is the
same confound flagged for the season effect.

### Outstanding

1. **Coverage for dark (0.57).** It is less biased AND tighter, so its intervals
   are too narrow for its own errors. This is where the effective-sample-size
   idea belongs, but NOT as uniform tempering, which failed.
2. The tangent's −1.54 is the emission bias localised earlier to twilight and
   night. Now unmasked and measurable.
3. `table_flat` and `logistic` have NOT been rescored. Their rejections rest on
   the biased kernel.
4. **The block samplers were not touched by this fix.** `TwilightFreeHier` may
   carry the same poleward lean; the notes already record a "surrogate mesh past
   the pole" defect, which is the same family. CHECK BEFORE PUBLISHING ANY
   HIERARCHICAL RESULT.

(The rescore script threw a cosmetic error in its own summary block, a `dcast`
expecting an area=FALSE column that this run does not produce. All 58 result
rows are written and complete.)

## THE SPIKE AND THE SLAB DO NOT EAT EACH OTHER (`diag_identifiability.R`)

ADB's concern: the shading arm explains observations BELOW expected and the slab
explains those ABOVE, so between them they could explain anything, leaving the
likelihood flat and position unidentified while every diagnostic still looks
healthy. Real risk, since lambda is already 3.1x looser than the measured
residual scale.

Measured per 12-hour knot as logL(true) − logL(true + 3 deg latitude), on six
deployments, corrected-floor table plus darkness regime at shade_ratio_dark 20:

SIGNAL, nats/knot (rows lam_mult, columns prob_slab_dark):

| lam_mult | 0.05 | 0.15 | 0.35 | 0.60 |
|---|---|---|---|---|
| 0.5 | 0.30 | 0.28 | 0.24 | 0.18 |
| 1.0 | 0.42 | 0.38 | 0.32 | 0.25 |
| 2.0 | 0.59 | 0.54 | 0.47 | 0.38 |
| 3.0 | 0.74 | 0.69 | 0.61 | 0.51 |
| 4.0 | 0.90 | 0.84 | 0.75 | 0.65 |

**No collapse anywhere.** Signal falls gently with the slab and never
approaches zero; the most generous cell (0.18) is still 1.5x the TANGENT's 0.12.
The feared corner does not exist in this range.

Better, the trade is favourable in both directions at once. Tightening lambda
raises signal 3-fold AND lowers the cost of an injected false-light reading
(3.50 to 2.84 at pslab 0.05), because a tighter spike makes the slab relatively
more attractive for genuine outliers. lam_mult 4 with pslab 0.35 gives signal
0.75 and false-light cost 1.03, against the current setting's ~0.40 and ~2.8.
Pareto-better on both axes. And it puts a number on the lambda finding: being
3.1x too loose is costing more than half the available day-length signal.

### BUT THE METRIC FAILS ITS OWN CONTROL — do not use the map to choose

I proposed this as a way to pick parameters before the sweep finishes. It cannot
do that, and I should not have implied it could. The control:

| arm | Argos median km | snr (signal/sd) | wrong-sign rate |
|---|---|---|---|
| tangent | **242** | 0.302 | **36%** |
| table_flat | **451** | 0.366 | **40%** |

`table_flat` has HIGHER mean signal and HIGHER signal-to-spread than the tangent
and fits nearly twice as badly. So per-window discriminating power at the truth
does not predict track accuracy, and the map cannot tell us whether the darkness
regime will help. Only the sweep can.

Likely reason: the metric is a LOCAL gradient evaluated AT the truth. The HMM
integrates the whole grid, so what matters is the global shape and any
systematic pull (table_flat carries +1.96 deg of bias, which a metric centred on
the truth cannot see).

The one summary that does order the two known arms correctly is the **wrong-sign
rate**, and by that measure every cell in the grid (37-42%) is slightly worse
than the tangent's 36%. That is n = 2 and far too weak to act on, but it is a
mild negative for the sweep and worth remembering when the results land.

**Keep the map for the manuscript as an identifiability result** — it answers a
real objection about the two mixture components — but NOT as a model-selection
criterion.

## STEP 1: THE LATITUDE PROFILE (`diag_latitude_profile.R`) — AND IT KILLS STEP 2

11,636 profiles: 29 tags, 1,160 windows, 4 responses, 4 zenith bands. For each
12-hour window the log-likelihood is evaluated across latitude at the TRUE
longitude, so this is the emission's entire opinion with no prior, no endpoints
and no smoother in it.

### The emission bias does NOT explain the track bias

| arm | emission bias | track bias | track km |
|---|---|---|---|
| dark | **+0.11** | +1.40 | 289 |
| logistic | −2.77 | −2.24 | 615 |
| table_flat | **−0.82** | **+1.96** | 451 |
| tangent | **−1.28** | **−0.14** | 242 |

`table_flat` REVERSES SIGN between emission and track. `tangent` has ten times
more emission bias than track bias. Only `logistic` roughly matches.

**So calibrating the response to zero the emission bias — my step 2 — is not
well posed.** Zeroing it would not zero the track bias. Recorded as a proposal
killed by its own gating test, before any fitting time was spent on it.

Why it fails is instructive: profile sd is **8.0 to 9.5 degrees**. A single
window is nearly uninformative about latitude, and the track's 2.5 degree
accuracy is built by multiplying ~480 such profiles. The bias of a product of
broad skewed densities is not the mean of their biases.

### Where the emission bias IS made: twilight and night, never day

| arm | day | twilight | night |
|---|---|---|---|
| dark | +0.75 | −0.70 | −0.64 |
| logistic | +0.70 | −3.79 | −4.19 |
| table_flat | +0.75 | −2.50 | −1.14 |
| tangent | +0.15 | −1.31 | −1.52 |

Daytime is small and positive for all four. The southward pull is manufactured
entirely in twilight and night. Profile sd by band confirms where the signal is:
day 12.7-13.6 (uninformative), **twilight 4.4-8.0 (the informative band)**,
night 11.1-12.8.

### The darkness regime wins on EVERY emission metric and still fits worse

`dark` has the smallest emission bias (+0.11 against −0.82/−1.28/−2.77), the
tightest overall profile (8.02), and a twilight profile of **4.4 against the
tangent's 6.6**, a 33% sharpening of exactly the band that carries latitude.
And it fits at 289 km against the tangent's 242.

**Therefore the emission is not the binding constraint.** The tangent's advantage
is downstream of the likelihood. That is the first hard evidence for where NOT
to look, after ten mechanisms died in the response.

### THE REDIRECT: the windows are combined wrongly, not evaluated wrongly

| arm | profile says sd | means actually scatter | ratio |
|---|---|---|---|
| tangent | 8.79 | 4.96 | **1.77** |
| table_flat | 9.47 | 5.33 | **1.78** |
| dark | 8.02 | 5.23 | 1.54 |
| logistic | 9.28 | 6.58 | 1.41 |

Per window the emission is CONSERVATIVE: it claims 1.5-1.8x less knowledge than
it has. At track level the posterior is OVER-confident (coverage 0.72 against a
nominal 0.95). Conservative inputs, over-confident output. **The reversal can
only happen in the combination**, where the HMM treats each window's emission as
conditionally independent given position.

They are not independent. Measured earlier (`diag_scattering.R`): residual
autocorrelation 0.26-0.33 across 1-4 hours, and **0.26 at 24 hours**, a diel
dawn-versus-dusk asymmetry that zenith-detrending cannot remove. Correlated
residuals mean the effective number of independent observations is far below the
count, so the product of profiles is too sharp, and any correlation shared
across days displaces it systematically.

**NEXT, and it is one parameter:** temper the per-window emission by an effective
sample size, `logL_window / n_eff`. Cheap to implement, cheap to sweep, and it
predicts both a coverage fix and a bias reduction if the bias comes from
over-confident windows being systematically wrong. This is the
variance-inflation item from the original list, but now with a mechanism and a
place to apply it.

**ALSO CHEAP, from data already on disk:** split windows by whether their
twilight observations are dawn or dusk and compare emission bias. The 24-hour ACF
says the asymmetry is there; if it maps onto latitude bias, that is the diel
term made concrete.

## Untested and worth doing

1. **Tag reuse across seasons.** 11 serials, 25 deployments, 3 in all three
   years. Does a tag calibrated in 2021 calibrate its 2023 deployment? If yes the
   recipe becomes calibrate-once-reuse. This is the biggest open opportunity and
   the data now support it.
2. **Why 2023 is worse than 2021.** Argos processing, endpoint derivation, or the
   animals. Confounded at present.
3. **Variance-inflation factor for honest intervals.** Reported latitude sd 1.42
   against actual 2.44, ratio ~1.7. One number, estimable from a double-tagged
   subset, and the only practical answer to the coverage problem given the
   central finding above.
4. **Family-wise pooling.** Mk9 pooled z50 91.92 (19 haul-outs), 18A 91.74 (1).
   Only 0.18° apart — the split may be unnecessary, but 18A rests on one fit.
5. Re-run the mean-vs-mode and occupancy-prior tests at n = 29.

---

## Reproducing

Every script in this directory has a header stating the question it answers.
`nes_common.R` carries the shared crosswalk, Argos reader and scoring.
Run from the package root. `fit_all29.R` is the gateway: it produces
`all29_knots.csv` and `all29_tags.csv`, which most downstream tests read.

Two process rules learned the hard way: **persist results to disk, never to the
terminal** (a `Select-Object -Last 5` in a launcher truncated a five-hour run's
entire output), and **append per unit of work** so an interrupted run resumes
instead of restarting.

---

## Multi-agent review, and the code changes that followed

Four reviews: adversarial (attack the kernel fix), audit (other engines), constructive
(path to manuscript), mediator (adjudicate + frame choices).

### What the review overturned in my own reporting
- **"16 of 16 tags by 1.0-1.4 deg" was wrong twice.** It is 8 tags scored under 2 arms,
  the arms correlate at r = 0.665 (effective n ~ 8), and 1.0/1.4 are the two ARM MEANS.
  Actual per-fit shifts: -0.385 to -2.177.
- **The empirical work cannot identify the coefficient.** A 0-3x sweep through the engine
  is monotone with no plateau at 1. ANY weight in ~(0.1, 2) reproduces 16/16 and passed
  the entire test suite; a half-Jacobian passed all four assertions. Only the derivation
  licenses the coefficient. A weight sensitivity analysis belongs in the paper.
- **The fix was not new.** `area_correction` is bit-for-bit `area_prior()` +
  `identity_rule()` -- same log(cos(lat)), same 1e-6 clamp, same three insertion points.
  Confirmed independently in the test suite (identical lat, lon and log_z).
- **`notes/topology/sbc_design.md` recorded the identical operator as "negligible".**
  Reconciled, not refuted: that was 41 daily knots, BOTH endpoints anchored, -50 deg in
  strong austral winter light. Ours is 273-494 knots at 12 h with a free end near
  equinox. The effect scales with knot count, endpoint anchoring and likelihood
  flatness. **A null from a sharply identified design says nothing about a weakly
  identified one** -- and that generalisation is exactly how the term got retired.
- **Flipping the default silently invalidated the fit-kernel SBC**, whose generator says
  in a comment "no cos-lat area weighting, so generator and fit share a movement model"
  and called TwilightFreeGrid() without the flag. Not previously on the outstanding list.

### What survived, and the evidence that actually carries the claim
Row sums constant to 4.5e-6 corrected vs 0.20376 nats uncorrected (analytic 0.2037);
boundary truncation 2.6e-4 and latitude-flat; `areaOFF` reproduces the archive exactly on
8/8. The mechanism evidence is **not** the tags: it is a hemisphere test built during
review (emission muted, symmetric -60..60 grid) where the shift is exactly antisymmetric
across the equator to 3 dp and ~0 at the equator. A generic "moves tracks south"
alternative predicts southward everywhere and dies on that table.

The darkness regime also survived the in-sample tuning attack: tuned 8 tags 5/8
p = 0.3125, **held-out 21 tags 18/21 p = 0.0002** -- stronger out of sample. Tuning was on
median km while bias_lat RISES with the tuned `ratio`, so it selected against the claim.

### Code changes made (all tested; installed build asserted before any job launched)
1. `area_prior()` DEPRECATED (warns), gains attr `tf_area_prior`. Its docs previously
   stated the "far too weak" verdict as fact -- rewritten to state the conditionality.
2. `TwilightFreeGrid()` now ERRORS (not warns) if `area_correction = TRUE` and an
   `area_prior()` term are both supplied. This double-count was live and silent, behind
   advice published in README:72 and the manuscript.
3. **Block sampler Jacobian fixed.** `area_delta_deg()` added to all three MH sites
   (`single_site`, `rb_polish`, `block_update`), gated on `ctx.spherical`. `metric="flat"`
   correctly owes nothing -- NOT because the prior cancels (it cancels in the defective
   spherical branch too) but because the degrees->km map uses a FIXED reference latitude,
   an affine map with constant Jacobian. Verified not inert: spherical sits equatorward of
   flat at +55/+40/-40/-55, sign flipping across the equator.
4. `log_z` docstring: valid to difference only same-grid AND same-`area_correction`.
   Measured shift -0.32 nats/knot vs ln cos(45) = -0.347 -- the whole shift is the
   normaliser.
5. Tests: removed the tautology (`expect_identical(run(FALSE,..), off)` compared a call to
   itself); relabelled the pure-R test as a derivation check, honestly; added the
   equivalence test that PINS the coefficient, the double-count error test, and a
   hemisphere-flip regression test for the block sampler.
6. `inst/sbc/sbc_grid_hmm.R`: generator now READS `AREA` instead of hard-coding "no
   cos-lat weighting", so it cannot silently desynchronise from the engine default again.
   Design caveat recorded in the script itself.

### Still open
- SBC re-run in flight (both arms, 100 reps, like-for-like with the recorded baseline;
  pre-correction figures preserved as `sbc_*_pre_area.png`). Prediction: the SPHERE arm
  should improve, since the fit now carries the sphere's measure. If it calibrates
  cleanly the sphere arm becomes the substantive claim and fit-kernel the control.
- Block-sampler fix changes 3.4/3.5 numbers; nothing rescored yet.
- Coastline sub-stochasticity: `makeGrid(mask="sea")` drops land cells and the kernel is
  unnormalised, so every coastline is an edge. Not triggered by the NES fits (no mask).
  Unquantified for coastal species.
- `mask_matrix`/`mask_extent` in the SMC path are accepted and documented but never
  referenced in the body. The mask is INERT.

### Shared blind spots the reviews were briefed into (highest value remaining)
1. **Argos was never questioned.** Tens of km of error and, worse, non-uniform sampling
   in latitude and behaviour (haul-out and surface-active periods over-represented). A
   latitude BIAS of 0.5-1.5 deg measured against that is exactly where a truth-set
   artefact hides. Nobody checked whether Argos coverage correlates with latitude within
   tracks. This is the most dangerous assumption in the campaign.
2. **We refused n=16 and then accepted n=29.** One species, one colony, one strategy,
   overlapping years -- and tag-reuse is on the cut list, i.e. it is KNOWN some tags
   recur. p = 0.0002 and p = 3.7e-9 assume an independence the cut list concedes is
   absent. The results will likely survive; the stated p-values will not.
3. **Kernel and emission are not separable.** 242 km was two errors cancelling; nobody has
   shown that fixing them sequentially converges where fixing them jointly would.

## SBC re-run: both arms PASS; the last "known latitude limitation" was the harness

Ran both arms on the corrected engine (100 reps, like-for-like with the recorded
baseline; pre-correction figures kept as `sbc_*_pre_area.png`).

First pass:  sphere PASS/PASS,  fit-kernel FAIL(lat)/PASS.

I proposed the spike normaliser as the next target and **was wrong** -- it is
already implemented and correct (`spike_normaliser`, matching the generator's
`m_lo + m_hi` exactly, added precisely because SBC caught the mu-dependence). The
SBC header caveat and the manuscript's "the spike density is not renormalised"
sentence are both STALE.

Second hypothesis, row-stochasticity, also **largely refuted by measurement**: row
sums are flat to 3.7e-8 across the interior of the SBC grid. Only the outermost
1 deg row deviates, by 0.067 nats, and symmetrically at BOTH edges -- it pulls
toward the domain centre and cannot produce a directional latitude skew.

**Actual cause: the ranking.** The fit-kernel generator samples CELLS, so truth sat
exactly on a cell centre while draws were jittered by U(+/- CELL/2); when the
posterior concentrated on the truth's own cell the rank was deterministically 0.5.
Signature: a near-vertical step at rank 0.5 with a deficit in both tails, present in
every fit-kernel figure before AND after the area correction. Latitude only,
because at -50 deg a longitude cell is 71 km against sigma = 80 (spreads over
several cells) while a latitude cell is 111 km (can contain the posterior). The
sphere arm never showed it -- `sphere_step()` returns a continuous position, which
is exactly why that arm passed.

Fix: jitter the RANKED mid-knot within its cell only. Chain, endpoints and emission
stay on cell centres, so movement remains bit-exactly the engine's.

**Result: fit-kernel latitude FAIL -> PASS, with no engine change. Both arms now
pass both coordinates.**

Also added `SBC_ARMS` env control so one arm can be re-run without paying for both.

### What this changes
- The manuscript's §3.5 "mild latitude miscalibration ... traces to a known
  approximation, the spike density is not renormalised" is wrong twice: the spike
  IS renormalised, and the residual was a harness artefact. That whole caveat goes.
- SBC now cleanly validates the grid HMM location posterior on both arms, including
  the physical-movement arm, which is the one that tests CORRESPONDENCE rather than
  self-consistency. That is a stronger calibration claim than the paper currently
  makes.
- It does NOT explain the real-data coverage shortfall (0.57-0.66). SBC passing and
  real coverage failing is consistent and is the interesting result: the model is
  calibrated for data generated from itself, so the shortfall is misspecification
  against real light residuals, not an implementation error. The temporal-
  correlation attribution (ACF 0.26-0.33) remains the live hypothesis.
- Still not evidence about the NES latitude BIAS, which is a separate quantity.

## Attacking the tangent's -1.54 deg: it is not the response, it is SHADING

Emission-only profiling (no HMM, no movement prior): simulate light at a known
position, profile the light log-likelihood over candidate latitudes at the true
longitude, take the displacement of the peak. Every run carries a matched-response
CONTROL, and the control earned its keep three times.

**My bug first.** `eval_logpk_grid`'s 7th argument is `shade_ratio`, NOT a cell
count (cells = `length(lon)`). Passing `length(grid)` made the upper arm ~80x
steeper than the lower one and drove every profile poleward -- the control peaked
+4 to +13 deg off truth and exposed it. Existing scripts are CLEAN
(`diag_latitude_profile.R`, `diag_identifiability.R`, `test-dark-regime.R` all pass
2, the correct value). The roxygen block documented six parameters and omitted
`shade_ratio` entirely; now documented, with the failure mode named.

### Results (all with the corrected call)
1. **Clean data, matched response: unbiased.** Control 0 to 1.25 deg across
   25-65 deg, both hemispheres.
2. **The tangent response does NOT explain -1.54.** On clean spike-and-slab data it
   biases POLEWARD at 30-50 N (+5 to +13 deg), mirror-symmetric across the equator,
   shrinking with |lat| (cor(|lat|,|bias|) = -0.90). Swept the zero-crossing over
   z0 = 98,100,104,108,112,116: **no column goes equatorward.** Both of my stated
   predictions (equatorward, growing with latitude) are REFUTED.
3. **Shading is the mechanism, and it has the right sign.** With a PERFECTLY
   MATCHED response, dive-bout shading alone drives a strong equatorward bias
   (p_shade 0.3 -> -3 to -6.5 deg; 0.6+ saturates a +/-40 deg window) while the
   unshaded control stays ~0. **So this is not a response-form problem and no
   choice of response fixes it.**
   Algebra: the spike charges a fully shaded observation `lam*(mu - x)`, i.e. in
   proportion to the EXPECTED light. Daily expected light grows with latitude in
   summer, so the model is charged for being far north.
4. **Tightening lambda makes it WORSE**, not better: at p=0.3, lat 45, bias -3.0 at
   the default 1/50 vs -16.0 at 1/8. The long-deferred "lambda is 3.1x too loose"
   retune is therefore NOT a bias fix. It may still matter for coverage; it is not
   this.
5. **Ratio-scale shading arm: REFUTED by its own control.** Replacing the lower arm
   with `lam_r*(1 - x/mu)` (scale-free, so shading costs the same fraction of light
   wherever you are) destroys the unshaded signal: control +40/+35 deg at lat 45/55
   where the additive form gives +2.0/+0.5. Not the fix, at least not naively.

### Where this leaves the bias
The target moved. It is not the tangent, not the zero-crossing, and not lambda. It
is that the shading arm's cost scales with expected light, which is a latitude
proxy. A fix has to decouple shading cost from expected light WITHOUT flattening
the latitude signal that the same arm carries -- and those two are not obviously
separable, since day length is exactly what both are made of. That tension is the
real problem and it is now stated sharply.

Magnitudes do not yet match reality: p_shade = 0.3 gives -3 to -4 deg against the
observed -1.54, so the effective shading fraction is below 0.3 and should be
calibrated against the tags' actual duty cycle rather than guessed.

## Depth: tested as a covariate, never as a CORRECTION

Question raised: did we ever estimate clear-sky light from observed light at recorded
depth? Answer: no. `diag_depth_bias.R` tested depth as an EXPLANATORY COVARIATE for
the bias and was dismissed on "96% of windows reach <10 m". The correction itself --
inverting Beer-Lambert to recover surface irradiance -- has never been tried.

**The depth channel exists on all 29 tags** (`depth_min`, `depth_max`, `temp_surf`,
per 30-min bin, no missingness; 2021 archive + fvilches archive).

Light is decimated by MAXIMUM, so the retained sample sits at `depth_min`, and
`exp(-k*depth_min)` is the surface light that survived into it. Pooled `depth_min`
median 2.0 m, 96.3% under 10 m -- which reproduces the "96%" behind the original
dismissal. **But <10 m is not "no attenuation": at k = 0.07 /m it is a 50% light
loss.** Between-tag spread is large: median `depth_min` 0.0 to 8.5 m, surviving
fraction 0.97 down to 0.55.

### Dose-response test: NULL at n = 29
cor(surviving light fraction, bias_lat), area-corrected fits:
- tangent: **+0.009, p = 0.962** (batch-adjusted slope +0.03, p = 0.977)
- dark:    -0.287, p = 0.131
- controls clean throughout (rmse_lon p = 0.41/0.57, median_km p = 0.30/0.43), so
  the covariate is not proxying data quality.

An n = 10 pilot on the 2021 archive alone gave +0.572, p = 0.084 and looked
promising. **It was noise** -- and the dark arm's wrong sign at n = 10 was the
warning that should have been heeded before reporting it.

### What the null does and does not establish
It refutes: per-tag VARIATION in measured attenuation driving per-tag VARIATION in
bias.

It does NOT refute the shading mechanism, which was demonstrated in simulation with
a matched-response control. **Every seal dives more or less constantly, so the
shading penalty is close to COMMON-MODE: a shared offset with little between-tag
variance.** A correlation across tags has no power against a common-mode effect --
which is precisely the blind spot that made the original month-to-month test null
too. Two different correlational tests have now failed for the same structural
reason, and neither is evidence about the mean.

**The only test that can settle it is direct: apply the correction and refit.**
`light_corrected = light_observed * exp(k * depth_min)`, k swept over the
open-ocean range 0.04-0.10 /m, and see whether the mean bias moves off -1.54
toward zero. That is a rescore of 29 tags, not a correlation.

## Logger reuse: effective n is 15, not 29 -- and the headline results survive it

The review flagged that paired tests across 29 deployments assume an independence
that the deferred "tag reuse" item concedes is absent. Now measured.

**The two deliveries do NOT overlap** (2021 archive 10, fvilches archive 19,
intersection 0), so no animal is double-counted. The fvilches RAW directory does
re-export every 2021 TOPP id, but those are not in its archive cache.

**But 29 deployments were produced by only 15 distinct loggers.** 25 of 29 (86%)
reuse one; three loggers were deployed three times:
  2190018 -> 2021023, 2022041, 2023030
  2190072 -> 2021032, 2022046, 2023035
  2190068 -> 2021033, 2022042, 2023036
Animals differ, so animal-level effects keep n = 29. Every claim in this campaign
is about the SENSOR -- response shape, dark current, shading, twilight -- and two
deployments sharing a logger share all of it. **Between-logger variance is 42% of
the tangent bias variance**, so this is not a technicality.

### Recomputed with the logger as the unit
| test | by deployment (n=29) | by logger (n=15) |
|---|---|---|
| tangent bias negative | 29/29, p = 3.7e-9 | **15/15, p = 6.1e-5** |
| dark beats tangent on abs(bias) | p = 1.9e-4 | **p = 1.2e-4**, 14/15 |
| distance | tie, p = 1.00 | tie, p = 0.56 (median -73 km) |
| rmse_lat | p = 0.325 | p = 0.083 |

Both headline claims survive, and the bias advantage STRENGTHENS at logger level.
**Quote the clustered figures.** 15/15 at p = 6.1e-5 is weaker than 29/29 at
p = 3.7e-9 and is the defensible one.

## Depth pairing: light and depth share a clock, but our BINNING does not

Raw Mk9 archives log light and depth on the SAME 4 s clock, same rows, identical
0.11% missingness -- so there is no logger sampling-rate mismatch. The 30-min cache
however stores light decimated by MAXIMUM and depth as depth_min over the bin,
computed separately, and the brightest instant need not be the shallowest one.

Measured within-bin discrepancy (depth at light-max minus depth_min): median
0.5-1.0 m, q90 1.5 m, >5 m on ~3% of bins. Implied Beer-Lambert error factor at
k = 0.07: median 1.036-1.073, q90 1.11, **max ~1e13**. Negligible in the median,
unbounded in the tail. The raw files can supply the correct pairing (depth AT the
retained light maximum) since both channels share a row.

Combined with the noise-amplification finding, the Beer-Lambert route is only
trustworthy where the correction is also small, which is a real limit on the idea
rather than a tuning choice.

## The band decomposition: twilight is where the signal AND the bias both live

Emission bias by zenith band, profiled at true Argos positions (deg latitude; sd is
the profile width, i.e. how sharply that band identifies latitude):

| arm | day | twilight | night |
|---|---|---|---|
| tangent    | +0.148 (sd 12.69) | -1.313 (sd 6.60) | -1.516 (sd 12.22) |
| logistic   | +0.704 (sd 13.15) | -3.785 (sd 8.04) | -4.194 (sd 11.12) |
| table_flat | +0.751 (sd 13.64) | -2.499 (sd 6.70) | -1.145 (sd 12.76) |
| dark       | +0.750 (sd 13.64) | **-0.695 (sd 4.41)** | -0.645 (sd 12.45) |

1. **Daylight is unbiased in every arm** (+0.15 to +0.75). None of the bias is there.
2. All the negative bias sits in twilight and night.
3. **Twilight is 2-3x sharper than either other band** (sd 4.4-8.0 vs 11-14). It is
   the only band that is both informative and wrong. Day is unbiased but nearly
   flat; night is biased and flat.
4. **This explains the darkness regime.** `dark` halves the twilight bias AND
   tightens the twilight profile (4.41 vs 6.60). It wins on bias because it fixes
   the band that carries the latitude. It also explains the rejections:
   table_flat -2.50 and logistic -3.79 are WORSE than the tangent in that band.

### Twilight diving: not a factor in this dataset
Median depth_min is 2.5 m in day, twilight AND night, with identical across-tag
ranges. No diel signal at the 30-min bin level. All four dose-response correlations
null (twilight attenuation, day attenuation, twilight depth, day depth; p = 0.42 to
0.97, deployment- and logger-clustered alike). So the twilight bias is NOT caused by
twilight-specific diving; it is the emission model's shape in the low-light regime.

### Consequence
The target is the twilight response shape, not shading and not depth. The supervised
"learn light-at-depth from expected light at Argos positions" route is also
deprioritised on its own terms: the band that needs help is the one where clear-sky
light is hardest to predict, and arm C already established that real residuals
reproduce 90% of the error.

## Depth-correction harness FAILED ITS CONTROL -- numbers uninterpretable
Uncorrected arm gave 712 km / bias +3.29 against the known rescore29 baseline of
276 km / -1.54 on the same engine and tags. Two causes, both mine:
1. **No within-family response pooling.** fit_rescore29.R pools responses within
   logger family; that step is what took the campaign 473 -> 208 km. I used raw
   per-tag fits.
2. **The HEADROOM factor was applied inside the likelihood parameters**, so
   max_light*1.4 also made lambda 1.4x looser and diluted the slab by 1/1.4. Given
   this session's finding that lambda strongly controls the shading-induced bias,
   that is not a neutral change.
The 725 km separation between arms is therefore the harness, not the shading arm.
Fixable, but deprioritised: the band decomposition says depth is not the lead.

## Re-tuning the darkness regime on TWILIGHT bias (emission-only, no HMM fits)

The three constants were tuned on median km, and the earlier sweep showed bias_lat
RISING with `ratio` -- so that tuning selected against the quantity that matters.
`dark_frac` had never been swept. Retuned against twilight-band emission bias, which
the band decomposition identifies as the only band that is both sharp and wrong.

Three rounds, each pushed until the optimum stopped moving. Admissibility required
BOTH that the twilight profile not widen AND that the DAY band not move.

**Round 2's apparent winner was killed by the day control.** dark_frac = 0.80 looked
best (twilight -0.199) but moved day bias +0.741 -> +0.285: dark_frac is a fraction
of max_light, so at 0.80 the "darkness" regime has swallowed most of the daylight
curve and is no longer a darkness regime. Disqualified. With dark_frac held at its
current 0.35 the day band does not move at all (spread 0.005 across 32 configs).

### Result (logger-clustered, 15 loggers / 29 deployments)
| config | twilight bias | twilight sd | night | day |
|---|---|---|---|---|
| tangent (no regime)        | -1.252 | 6.44 | -1.681 | +0.356 |
| dark, CURRENT f0.35_r100_p0.35 | -0.773 | 4.37 | -0.813 | +0.741 |
| retuned f0.35_r400_p0.03   | **-0.221** | **3.28** | -0.366 | +0.741 |
| raw optimum f0.35_r1600_p0.01 | -0.188 | 3.20 | -0.183 | +0.742 |

Converged, not edge-limited: ratio 400 -> 800 -> 1600 buys only -0.221 -> -0.207 ->
-0.200 at pslab 0.03. That is the slab saturation predicted in advance -- a sharp arm
cannot cost more than -log(pslab/max_light), so it stops improving rather than
running away.

### Recommended setting: dark_frac 0.35, ratio 400, prob_slab_dark 0.03
NOT the raw optimum. prob_slab_dark = 0.01 is below the MEASURED night contamination
rate (median 1.8% per tag, worst 8.7%), so it would under-absorb on most tags and on
the worst one badly; it wins here only because the profile diagnostic scores bias at
the true position and never has to survive a contaminated night. 0.03 keeps ~95% of
the gain and stays above the measured median. dark_frac stays at 0.35 because that is
where the day control is clean.

Versus the tangent this is an 82% reduction in twilight bias and a 49% tighter
twilight profile; versus the shipped darkness setting, 71% and 25%.

### Not yet established
This is EMISSION-ONLY. It says the likelihood's twilight bias falls; it does NOT say
the fitted tracks improve. The HMM rescore at n = 29 is the test, and the campaign has
already seen an emission-level prediction fail at the fit level. Also unresolved: 20
warnings recurred in every round and were never inspected, and twilight rests on 223
windows across 24 tags (vs 840 for day), so it is the sharpest band and the thinnest.

## Grid search: calibration vs shading. It is SHADING (lambda), not calibration.

Two knobs, both tunable against Argos, emission-only (no HMM), 45 configs x 29 tags
x 25 knot windows:
  delta     rigid shift of the tangent along the zenith axis, c(a + b*delta, b).
            This is literally "what light level do we call civil twilight", and it
            leaves the daytime slope untouched so it cannot be confounded with gain.
  lam_mult  multiplier on the spike's lower-arm rate, i.e. the price of observing
            light below the expected curve. shade_ratio held at 2.

CONTROL 1 PASSES: delta 0, lam 1 reproduces the shipped tangent exactly --
day +0.357 (known +0.36), twilight -1.215 (-1.25), night -1.658 (-1.68).

### Result: lambda owns the bias
bias range across delta 0.76 deg; across lambda **1.76 deg**. At delta = 0:

| lam_mult | 0.25 | 0.5 | 1 (shipped) | 2 | 4 |
|---|---|---|---|---|---|
| bias | -1.55 | -1.73 | -1.31 | -0.60 | **+0.03** |
| profile sd | 12.1 | 10.7 | 8.77 | 6.73 | **5.04** |

Zero crossing at lambda ~3.7x TIGHTER than default. CONTROL 2 passes convincingly:
the profile SHARPENS as lambda tightens (8.77 -> 5.04), so the bias reduction is not
bought by flattening the likelihood.

**Convergence with an independent measurement:** the campaign separately measured the
default lambda as 3.1x too loose against the residual scale. This grid, from Argos
truth and knowing nothing about that, puts the zero at ~3.7x. Two unrelated lines of
evidence landing together is the strongest result of the campaign.

delta's effect is weak, monotone in the UNHELPFUL direction (later twilight -> more
southward), so no calibration choice within +/-6 deg can fix the bias.

NOTE this CONTRADICTS the earlier synthetic which said tightening lambda makes the
shading-induced bias worse. The synthetic used an invented shading process; this uses
the real records against real truth. Prefer the grid.

### Leave-one-logger-out (15 folds; deployment-level CV would leak via shared sensor)
lam_mult = 4 chosen in **all 15 folds** -- a stable, real optimum. delta's choice was
scattered (-1.5 to 6), i.e. noise, consistent with its weak surface.
  in-sample  |bias| 0.026
  held-out   |bias| 1.264  vs shipped 1.481, p = 0.277
The gap 0.026 -> 1.264 is the overfitting, and it has a clear cause: lambda removes
the MEAN bias but per-logger bias remains scattered -2.6 to +3.5. **A single global
lambda cannot fix per-logger variation**, which is where the remaining error lives.

### Stage 2 running: 4 fit-level points, 116 fits
lam 1.0 (reproduction control, must return 276 km / -1.543), 2.0 and 4.0 (predicted
better), and **0.5 (predicted WORSE)**. The 0.5 point is the one that matters: if the
fits also worsen there, the emission ranking transfers; if they improve, the ranking
is inverted and stage 1 says nothing about fits.

## The retuned darkness regime FAILED at fit level -- emission tuning does not transfer
Reproduction control passed exactly (max |km diff| 0.0, max |bias diff| 0.0000).

| arm | km (mean) | bias | abs(bias) | rmse_lat | cover |
|---|---|---|---|---|---|
| tangent | 276 | -1.543 | 1.543 | 3.449 | 0.656 |
| dark (shipped) | 272 | +0.458 | 0.734 | 3.041 | 0.568 |
| dark_retuned | - | +1.032 | 1.231 | 4.582 | 0.483 |

Retuned is SIGNIFICANTLY WORSE: abs(bias) better on 9/29 p = 0.006; rmse_lat
p = 0.0013; coverage p = 0.017. Logger-clustered, better on 3/15, p = 0.030.

**Why it failed, and it matters:** at the shipped setting the emission profile says
twilight bias -0.773 (southward) while the FIT says +0.458 (northward). They disagree
in SIGN. Reducing the emission's negative bias therefore pushed the fit further
positive and overshot. The band decomposition remains valid as DESCRIPTION; it is not
usable as a tuning objective.

## CORRECTION: "the kernel fix made accuracy worse" depends on the summary statistic
242 and 276 are MEANS of per-tag medians (fit_rescore29.R's convention).
| summary | archived (uncorrected kernel) | corrected |
|---|---|---|
| mean of per-tag medians | 242 | 276 (worse) |
| median of per-tag medians | 253 | **241 (better)** |
The direction REVERSES under the robust summary, so the regression is driven by a few
large-error tags, not a general degradation. "276 is the honest figure" was one
statistic's answer stated as the answer, and the write-up must report both.

## There is NO good default lambda -- it is genuinely per-tag

Per-logger lambda that zeroes that logger's own latitude bias (delta = 0, emission
grid, 15 loggers):
- **only 5 of 15 have an optimum inside the swept 0.25x-4x window.** 8 need >4x,
  2 need <0.25x. So the per-logger optimum spans at least 16x, probably more.
- among the 5 interior ones the spread is still 3.9x (0.62 to 2.39), sd of
  log2(lam_opt) = 0.82.
- best single default (lam = 2) leaves mean abs(bias) 1.114, worst logger -2.87, and
  only **6 of 15 loggers within 1 degree**. lam = 4 manages 8 of 15.

So a default cannot work and the problem becomes estimating lambda WITHOUT truth.

### Why the obvious criterion already failed, and what it constrains
Recorded earlier: "per-tag lambda by marginal likelihood: 209 -> 283 km, 0/10
interior optima, rejected". Same pathology as above -- log_z is MONOTONE in lambda.
Cause: the light residuals are temporally autocorrelated (ACF 0.26-0.33 over 1-4 h),
so the likelihood over-counts independent information and always rewards a sharper
emission. Pseudo-replication.

**This rules out any criterion that treats observations as independent** -- log_z,
AIC, naive one-step-ahead scoring will all be monotone. A working criterion must
respect the autocorrelation.

### Candidate truth-free criteria, ranked
1. **Dawn/dusk self-consistency (preferred).** Day length is measured twice daily by
   nearly independent halves of the record. Fit latitude from dawn-side observations
   only, then dusk-side only, form z = (lat_dawn - lat_dusk)/sd_claimed. Correct
   lambda gives sd(z) = 1; too tight gives >1, too loose <1. Guaranteed interior
   optimum, no truth, no template, no held-out data, and it targets exactly what
   lambda is -- a calibration parameter. VALIDATABLE here against per-logger
   Argos-optimal lambda.
2. **Block cross-validation / filtered-vs-smoothed.** Hold out a contiguous day,
   predict it from the rest. This is the implementable form of the user's
   "forward-backward agreement" idea and it respects autocorrelation.
3. **Movement realism.** Weakest as an OBJECTIVE: minimising implausible movement
   risks re-expressing the movement prior, i.e. tuning the likelihood until the
   posterior matches the prior. Stronger in the user's species-specific form
   (elephant seal out-and-back, shearwater coast-hugging north then sprint south)
   because that is external structure the isotropic Brownian model does not contain.
   Best used as a CHECK on the answer, not as the thing optimised.

### Prerequisite now running
Validating any truth-free criterion needs per-logger Argos-optimal lambda for all 15
loggers, not 5. Wide sweep launched: 19 lambdas over 0.0625x-32x, delta fixed at 0
(settled as noise), 29 tags.

## Per-tag lambda spans 256x. No default can work.

Wide sweep, 19 lambdas over 0.0625x-32x, delta = 0, 29 deployments / 15 loggers.
Per-logger lambda that zeroes that logger's own latitude bias:
  **interior optima 11 of 15 | median 0.20 | IQR 0.14-2.24 | range 0.07-18.41
  | fold-spread 256x | sd of log2(lam_opt) = 2.80**
Four loggers still optimise outside even this range.

Cost of the best global default (lambda = 2): mean abs(bias) 1.114, worst logger
-2.87, and only 6 of 15 loggers within 1 degree. Across the WHOLE ladder no global
value does better than ~1.11, and the best coverage any value achieves is 8/15
within 1 degree (lambda 2.83-4). The mean does have an interior optimum near
lambda = 2, so a default is not arbitrary -- it is simply far too coarse.

**The easy truth-free handle fails**: profile sd at the default does NOT predict a
logger's optimal lambda (r = -0.14, p = 0.68, n = 11), so the width of the likelihood
carries no usable information about how much it should be sharpened.

## Dawn/dusk criterion built (fit_dawndusk.R + anal_dawndusk.R)
Ground-truth-free. Each day contains two transitions; estimate latitude from the
rising half and the setting half separately, then
    z = (lat_dawn - lat_dusk) / sqrt(sd_dawn^2 + sd_dusk^2)
Calibrated lambda gives sd(z) = 1; too tight over-confident (>1), too loose (<1). It
has an interior optimum BY CONSTRUCTION, which is exactly what log_z lacked, and it
compares two halves of the SAME day so shared slow structure is differenced out
rather than counted twice -- the autocorrelation that made log_z monotone.

Truth-free throughout: longitude is estimated per window from the light alone (2-D
coarse pass, once at lambda = 1, reused across lambda since the carrier phase should
not depend on the shading rate -- the spread is reported as a check), and the
dawn/dusk split is the sign of d(zenith)/dt at the estimated position.

mean(z) is a SEPARATE diagnostic: a systematic dawn-vs-dusk offset is an asymmetry
(wrong response on one limb, or a clock error), not a dispersion problem, and no
lambda can fix it.

STATED CAVEAT: the two halves of one day share weather and behaviour, so they are not
fully independent. That makes them agree more than chance, deflating sd(z) and biasing
the chosen lambda TIGHT. Hence the criterion must be VALIDATED against the
Argos-optimal lambda per logger (anal_dawndusk.R does this), not trusted alone.

## Coarse-then-refine grid search (R/refine_grid.R) -- built, and its first design FAILED

Motivation: longitude is a carrier phase and tightly identified. Measured over 2184
knots on 5 tags, the smallest contiguous longitude band holding 1-eps of the marginal
posterior, out of 100 columns:
  eps 1e-2: median 8, max 11 | 1e-4: 11/16 | 1e-6: 15/21 | 1e-9: 18/25
So most of a generous domain never holds mass.

### The first design failed its acceptance test, and the reason is worth keeping
Naive coarse-then-refine gave 44-56x speedup but moved the answer by 2.49 deg of
latitude -- an approximation, not an optimisation. Diagnosis:

| latitude envelope (1-1e-6 of mass) | |
|---|---|
| fine (1 deg)   | 24.5 to 50.5 |
| coarse (4 deg) | **38.0 to 38.0** |

**A coarse grid does NOT blur the posterior, it SHARPENS it.** Four degrees of
latitude is a large enough change in day length that one cell beats its neighbours
outright, so the coarse posterior collapsed onto a single cell and its MAP latitude
was constant at 38.0 while the fine MAP ranged 36.5-42.5. A coarse pass therefore
cannot bound a fine posterior without being deliberately broadened.

### The fix, both parts principled rather than tuned
1. `coarse_inflate = 3` -- the coarse pass runs with its movement scale inflated,
   which broadens its posterior back into an honest envelope.
2. `cell_margin = coarse.size` -- a fit cannot localise better than its own cell
   size, so a box from a 4 deg fit is padded 4 deg on top of the 5-sigma movement
   margin.
Both are documented in the function WITH the measured failure, so they are not
mistaken for over-caution later.

### Result: acceptance test passes, and the honest speedup is ~2.6x
max diff lat 0.50 deg, lon 0.68 deg (both inside one cell) at every eps.
  full domain 100x50 at 1 deg: 328 s, 5000 cells
  eps 1e-4: 122 s (2.7x) | 1e-6: 127 s (2.6x) | 1e-9: 153 s (2.1x)
The 44-56x was buying speed by clipping. 2.6x matches the 3-4x predicted from first
principles beforehand (emission is cells x obs, forward-backward is cells x
neighbours and is the larger term, and the band must exceed the posterior by the
5-sigma movement radius or the kernel itself gets truncated).

Box remains conservative: lat 22-65 vs a fine envelope of 24.5-50.5, so
`coarse_inflate = 2` would likely recover speed -- but that must be tested on several
tracks, not tuned on the one that motivated it.

Also added: `TwilightFreeGrid()` now returns `diffusion`, `step_hours` and
`area_correction`, so a margin can be sized in km from a fit alone.
`log_z` is NOT comparable across different domains (a sum over cells) -- documented
alongside the same warning for `area_correction`.

Suite: 302 pass, 0 fail, 7 skip (raster/SGAT absent).

## lambda_scale built: per-knot shading rate (Rust + R), and a simulation study

Engine: `run_grid_hmm` gained `lambda_scale`, a per-knot multiplier on lambda. Empty
= constant and bit-identical to before. `LightMix` was built ONCE per fit, so the
shading rate was necessarily constant across a whole track; it is now built per knot
when a schedule is supplied. Validated length and positivity in both R and Rust.

R: `TwilightFreeGrid(lambda_scale = NULL)`, plus `declination_lambda_scale()` and
`solar_declination()`. The schedule is normalised to geometric mean 1, so it
REDISTRIBUTES the rate through time without changing its level -- otherwise a change
of schedule would silently change the overall rate too and the two could not be told
apart. Default slope -0.081 log2 per degree of |declination| (the measured value);
`slope = 0` recovers constant exactly, which makes it a clean control.

Tests (24, all passing): NULL == constant 1 bit-identical; a constant scale c is
EXACTLY equivalent to scaling likelihood_params[1] by c (this fixes the argument's
MEANING -- without it the multiplier could be doing something adjacent); a varying
schedule changes the fit; input validation; the schedule is tighter at equinox than
solstice; slope 0 is exactly flat; max_ratio caps; solar_declination signs.

**A test I wrote wrong and had to fix.** The "is it inert?" test first asserted the
MAP track changes, and it FAILED against a working implementation: the MAP is snapped
to cell centres, so on a well-identified synthetic it can be bit-identical while the
likelihood underneath has changed a lot. Assert on `log_z` and the posterior instead.

Suite: 326 pass, 0 fail, 7 skip.

### Simulation study launched (sim_lambda_schedule.R)
The real data cannot test this: all 29 trips share one phenology, so time of year,
latitude and behaviour are confounded. Simulation breaks that -- **the shading process
is CONSTANT in time by construction**, so if the best lambda still varies with
declination it can only be the information geometry of day length.

PART A (must pass first): does simulated data with constant shading reproduce
lambda_opt falling with |declination|? If not, PART B is untestable.
PART B: scheduled vs constant at MATCHED LEVELS, best-of-each compared -- a schedule
must beat the best CONSTANT lambda, not merely the shipped one.

Scenarios: north_span (solstice->equinox, as the seals do), south_span (mirror),
equator (crossing, least day-length information), and **solstice_only as the control
-- declination barely moves there, so the schedule must do NOTHING**. A schedule that
helps even there is not doing what it claims.

NOTE: install is blocked (the stage-2 job holds the DLL), so the simulation runs
against the SOURCE build via load_all. Verified the feature is present in that build
before launching -- the failure mode recorded earlier this session was launching a
long job against a stale INSTALLED binary.

## STAGE 2 COMPLETE: lambda is a bias-vs-coverage dial, and the shipped value is not optimal

29 tags x 4 arms. Reproduction control exact (max |km diff| 0.0, max |bias diff| 0.0).

| lambda | km | bias | abs(bias) | rmse_lat | coverage |
|---|---|---|---|---|---|
| 0.5  | **266** | -1.691 | 1.691 | **3.302** | **0.779** |
| 1.0 (shipped) | 276 | -1.543 | 1.543 | 3.449 | 0.656 |
| 2.0  | 298 | -1.188 | 1.239 | 3.609 | 0.521 |
| 4.0  | 328 | -0.646 | **0.927** | 3.983 | 0.421 |

**cor(emission bias, fit bias) = +0.979** over all four arms. The emission surface
DOES transfer for lambda -- unlike the darkness-regime parameters, where emission and
fit disagreed in sign. So emission-only tuning is parameter-specific, not generally
invalid, and the two-stage design was the right call.

The lam0.50 falsification point behaved exactly as predicted: worse on bias for 27/29
tags, p < 0.0001. Predictions held at both ends.

### The finding that matters
**lambda = 0.5 beats the shipped setting on BOTH accuracy (266 vs 276 km) and
coverage (0.779 vs 0.656)**, costing 0.15 deg of extra bias. The shipped lambda = 1 is
optimal for NOTHING -- it sits between the accuracy/coverage optimum and the bias
optimum. Tightening buys bias and pays in km, rmse_lat AND coverage, monotonically
across the whole range.

**Coverage 0.779 is the best figure this campaign has produced.** The under-coverage
that has dogged the manuscript (0.57-0.66 vs nominal 0.95) is substantially a
CONSEQUENCE OF LAMBDA BEING TOO TIGHT, not an unexplained miscalibration. That turns
a stated limitation into a tuning choice, and it should change what section 3.5 says.

NOTE this also retires the long-standing "lambda is 3.1x too loose" recommendation in
the opposite direction: the residual-scale measurement was right that the SPIKE is
looser than the residuals, but acting on it costs accuracy and coverage. On the
evidence, lambda should go LOOSER than shipped, not tighter, unless bias is the only
thing being optimised.

### Queued
lambda 0.25 and 0.125: coverage is still climbing at the loose end and accuracy has
not turned over, so the optimum is outside the swept range on the loose side -- the
same edge-of-grid failure this campaign has hit three times. Run after the simulation
frees the CPU.

## SIMULATION RESULT: mechanism confirmed (Part A), schedule NOT justified yet (Part B)

### PART A -- the mechanism is real and geometric. CONFIRMED.
Simulated data with a shading process CONSTANT IN TIME by construction still shows
lambda_opt falling with |declination|:
  simulated r = -0.335, p = 2.4e-07, slope -0.104 per degree
  real data  r = -0.31,  p = 5e-5,   slope -0.081 per degree
No behaviour, no phenology, no latitude confound -- so the relationship is the
information geometry of day length, which is exactly what the 29 real deployments
could never establish (they share one phenology). **This is a genuine finding.**

### PART B -- no significant fit-level benefit. NULL.
Paired within track AND level (schedule minus constant):
| scenario | median d abs(bias) | p |
|---|---|---|
| equator | -0.123 | 0.071 |
| north_span | -0.057 | 0.86 |
| south_span | +0.081 | 0.44 |
| solstice_only (CONTROL) | +0.008 | 0.44 |
| spanning pooled | -0.085, better 53/96 | **0.58** |

The control is correctly null, so the schedule is not spuriously helping everywhere.
It simply is not helping. The best-of-each table (4.228 -> 4.189, 2.333 -> 2.101,
2.858 -> 2.630) looked favourable but compares minima of aggregate curves; paired,
there is nothing. Do not quote best-of-each -- it flatters the schedule.

### Why the null is not decisive: the simulator is not calibrated to reality
At level 1, constant: simulated abs(bias) 3.82, rmse 6.66, coverage 0.24
                      real       abs(bias) 1.54, rmse 3.45, coverage 0.66
Two to three times harder on every metric. A schedule redistributing a factor of 3.7
in lambda has little to work with when the fit is that poorly determined. Likely
culprits: P_SHADE = 0.45 with full attenuation, and DIFF = 90 km/sqrt(day) over 210
days, which wanders far more than a real seal.

**Next: calibrate the simulator to reproduce the real error scale, then re-run Part B.**
Until then the schedule is unproven, not disproven.

### Recommendation
KEEP `lambda_scale` -- it is a capability, defaults to NULL, bit-identical unused, and
the engine could not previously express a time-varying rate at all. Do NOT make
`declination_lambda_scale()` a default. The mechanism is established; the benefit is
not.
