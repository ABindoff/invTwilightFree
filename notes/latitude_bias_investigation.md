# The residual latitude bias: what it is, what it is not, and what to do next

Status as of 2026-08-18. Branch `fix/light-response-autodetect`.
**Read sections 0, 0a, 0b and 0c first, in that order.** 0a and 0b correct section 0;
0c is the closest this dataset comes to an answer.

This note exists so that the manuscript can be edited from a written record rather
than from recollection. Everything here is reproducible from scripts in
`scratch/nes_calibration/` (gitignored, but present on the working machine) and
from the numbers quoted below.

---

## 0. CURRENT STATE (2026-08-17, evening) — read this first

This note grew by accretion and contains retractions. Everything below is superseded
by this section where they conflict.

### The engine is sound on correctly-specified data

The decisive control is the battery's **noisy** arm: clear-sky light plus a draw from
the engine's OWN spike-and-slab mixture (KS-verified against `LightMix`), so the
emission is correct by construction, while the TRACK is a real Argos track rather
than a draw from the movement prior — the situation on real data.

| | honest synthetic | real data |
|---|---|---|
| latitude bias | **-0.51 deg** (6/6 negative, p = 0.031, one track gives half) | -1.54 deg |
| median error | 189 km | ~244 km |
| latitude coverage | **0.996** | **0.64** |

**The coverage failure does not reproduce at all** when the emission is right. So
roughly a degree of the real-data bias, and essentially all of the coverage problem,
is the emission not matching REAL light — not the inference, geometry, or prior.

### What was tested and does NOT help

Censoring the emission, the Ito drift correction, tightening or splitting the
movement prior per axis, changing the point-estimate readout, unclamping the
expected curve. All measured; all neutral or worse. Two were falsified only after
being implemented. **Do not re-run these.**

### The gauge is real, and your endpoints already suppress it

Fisher analysis of the emission (553 one-day windows): condition number **3.8e4**,
and latitude is essentially THE soft direction (weight 0.881 in the softest
eigenvector, 0.005 in the stiffest), trading against the threshold `z50` at 1.4-4.4
deg of latitude per degree, with the exchange rate **changing sign across the year**.
That explains why axis-aligned hyperparameter sweeps returned null or produced
compensating shifts: they were moves along a level set.

But the fitted sensitivity is only **0.07 deg/deg**, ~20x smaller, because both
endpoints are pinned to truth. The pinned deployment and recovery positions are
already acting as closure constraints and suppressing the gauge mode.

### THE DELIVERABLE: identify the calibration by pooled evidence, not at the haul-out

Gauge inflation (cost in latitude of not knowing the calibration) falls from **1.81**
on a 5-day window to **1.08** over a full 240-day deployment, and **1.01** pooling six
tags, with `se(z50)` shrinking as 1/sqrt(N) to 0.045 deg. The seasonal contrast over
a full deployment identifies the calibration.

`fit_light_response()` estimates `z50` on a **15-day haul-out window** — under 1 deg
of declination span, where inflation is 1.7 — and freezes it for a deployment
spanning 47 deg. The identifying information is in the data and is discarded.

**Validated fix, no engine change required.** Profile `log_z` (already returned) over
`z50`, summed across tags, and take the maximiser:

| dz50 | -2.0 | -1.0 | -0.5 | **0.0** | +0.5 | +1.0 | +2.0 |
|---|---|---|---|---|---|---|---|
| rel log_z | -48.0 | -13.7 | -4.3 | **0.0** | -0.9 | -7.1 | -35.5 |
| median km | 480 | 314 | 245 | **189** | 194 | 210 | 350 |

Pooled maximiser lands **exactly on the truth (+0.00)** while per-tag maximisers
scatter with sd 1.03 and one tag (2023037) runs to the boundary and identifies
nothing alone. Accuracy is optimised at the same point.

**Calibration error costs ACCURACY, not bias**: ~125 km per degree of `z50`, but only
0.07 deg of latitude. Given the recorded 5.6 deg spread in haul-out `z50` estimates
across tags, this is a large accuracy cost currently being paid silently.

### What remains unexplained

About 1 deg of real-data latitude bias and the whole of the coverage failure, both
localised to emission fidelity against real light. The battery cannot investigate
this: it generates from the model. Next step is characterising the real residual
distribution (observed light minus fitted expectation at Argos-known positions)
against the spike-and-slab the model assumes.

## 0a. THE PROBLEM SPLITS IN TWO (2026-08-18) — depth explains the seasonality, not the bias

Independent review (two Fable agents) plus verification found that `depth_min` — the
shallowest depth in each 30-minute decimation window, present in the cached panel and
never regressed against the emission residual at observation level — is a large,
structured term in the emission.

**Verified independently.** In the graded zenith band (82-100 deg, where the tangent
is not clamped; it IS clamped for 94% of the record):

| depth_min | share | median residual | slab-scale fraction |
|---|---|---|---|
| <2.5 m | 58.1% | +0.028 | 0.002 |
| 2.5-10 m | 37.7% | +0.022 | 0.000 |
| 10-25 m | 0.3% | **-0.103** | 0.017 |
| 25-100 m | 1.0% | **-0.246** | 0.333 |
| >100 m | 2.9% | **-0.492** | 0.884 |

It is a GATE, not a gradient — flat to 10 m, then a cliff. **95.1% of slab-scale
events (residual < -0.3 of max_light) occur below 10 m.** Depth explains **51.3% of
within-tag residual variance**. And 41.9% of graded windows never come shallower than
2.5 m, so the 30-minute MAXIMUM is not recovering surface light.

**That refutes the recorded "depth attenuation: dead" verdict**, whose stated reason
was that 96% of windows reach within 10 m of the surface. Reaching 10 m is not
reaching the surface, and the gate sits exactly at 10 m. The two earlier depth tests
also had specific defects: `diag_depth_bias.R` aggregated to MONTH level, destroying
a 30-minute signal; `fit_depth_corrected.R` inverted Beer-Lambert in LINEAR units on
a LOG-scale sensor.

### What gating at 10 m does, and does not do

**DOES remove the seasonal drift and the slab** (`diag_depth_gate.R`):

| | ungated | gated |
|---|---|---|
| NH winter | -0.373 | +0.633 |
| equinox | +0.053 | +0.594 |
| NH summer | +0.537 | +0.557 |
| **seasonal swing** | **0.910** | **0.076** |

A 92% collapse, in zenith-equivalent degrees. Twilight-band tail mass falls from
0.030/0.026/0.009 to 0.001/0.000/0.000.

**DOES NOT remove the latitude drift.** Pre-registered requirement was to at least
halve it; the within-tag slope goes -0.0850 -> -0.0707 zenith-equivalent deg per deg
of latitude, a **17% reduction. PREDICTION FAILED.**

### So the problem is two problems

1. **Seasonal emission drift = diving behaviour.** Removable with the tag's own depth
   channel. Explains the slab entirely.
2. **Latitude drift = still unexplained**, and it is the one that maps onto the bias:
   -0.085 zenith-equivalent deg per deg latitude, times the measured exchange rate of
   1.4-4.4, gives **-0.12 to -0.37 deg of latitude bias per degree of northward
   excursion**. Over a 15 deg excursion that is -1.8 to -5.6 deg at the northern
   extreme — the right order, and structured as the bias is.

### Immediately actionable

**`prob_slab` should be near zero on depth-gated light**, not 0.10. Measured tail
after gating is 0.000-0.001: a ~100x overstatement of contamination, and unlike
`lambda` it is not compensating for anything.

### The open question, sharpened

**What makes twilight light drift dimmer, relative to a zenith-only model, as the
animal moves north — after depth is controlled for?**

## 0b. A correction to section 0's framing

Section 0 concluded "the coverage failure does not reproduce on the honest synthetic,
so the gap is emission fidelity to real light". That is sound for COVERAGE and was
wrongly extended to BIAS. Adversarial review showed the bias is substantially
accounted for by the cos(latitude) area tilt acting on a WIDE PER-KNOT LIKELIHOOD,
with magnitude scaling as prior force / effective latitude precision — the honest
synthetic does not escape it, it simply has ~3.3x more precision.

- Battery's own per-tag law `bias = 0.03 - 0.147 * var(err_lat)`, extrapolated to the
  real error variance, predicts **-1.37 deg** against an observed **-1.54**.
- Real data, area OFF vs ON: bias -0.166 -> -1.543, i.e. the correction moves bias
  **south by 1.38 deg on 28/29 tags**, and the per-tag shift scales with posterior
  latitude variance (r = -0.65, p = 1.3e-4).

**So the elimination row "the area correction is not a bias source" is WRONG.** Both
are true: the uncorrected kernel drifts poleward AND the corrected one carries an
equatorward tilt. Keeping the correction is still right; calling it inert was not.

This also resolves the "80x too small" dismissal recorded earlier. The tilt shift is
`-sigma^2 tan(lat)`, and the right sigma is the **per-knot LIKELIHOOD** width (6-9
deg, in the Fisher output), not the smoothed posterior width (1.1 deg). That gives
~0.85-1 deg, and the arithmetic closes.

## 0c. THE ANSWER, as far as this dataset can give one (2026-08-18)

**Real twilight light along these tracks is not a pure function of solar zenith.**
Relative to a fitted zenith-only response it gets systematically dimmer as the animal
moves north, at the surface, independent of diving:

    within-tag slope = -0.0707 zenith-equivalent degrees per degree of latitude
                       (se 0.0032, t = -21.9), depth_min <= 10 m

("Zenith-equivalent" = residual divided by the fitted response slope, so +1 means the
light is one degree of zenith brighter than the model expects.)

### The control that makes this a finding and not an artefact

The identical analysis on the synthetic battery — light GENERATED from a pure zenith
function, so it has no atmosphere, no cloud and no latitude dependence whatsoever:

| light | latitude slope |
|---|---|
| synthetic, noiseless | **+0.0009** (se 0.0013, t 0.7) |
| synthetic, model-correct noise | **+0.0010** (se 0.0094, t 0.1) |
| **real, gated at 10 m** | **-0.0707** (se 0.0032, t -21.9) |

Zero in synthetic, 22 sigma in real. Response fitting, the tangent, the clamping and
the zenith-equivalent conversion create NO slope on light that has none.

### It is not depth. Eliminated three ways

| gate | kept | slope |
|---|---|---|
| none | 100% | -0.0850 |
| 10 m | 96% | -0.0707 |
| 5 m | 75% | -0.0746 |
| 2.5 m | 58% | -0.0711 |
| **1 m** | **37%** | **-0.0816** |

Keeping only observations where the animal came within ONE METRE of the surface
leaves the slope as large as ungated. Depth inside the gate does not vary
monotonically with latitude (1.45 / 3.06 / 3.08 / 1.82 m), and controlling for it
linearly removes 11%.

### It is the right size to be the whole bias

-0.0707 over a typical 15 degree excursion is -1.06 zenith-equivalent degrees; at the
measured Fisher exchange rate of 1.4-4.4 deg latitude per degree of z50 that is
**-1.49 to -4.67 degrees**, bracketing the observed **-1.54**.

### What it is NOT: Rayleigh scattering in air

Rayleigh acts through AIR MASS, a function of zenith alone — the optical path at a
given z is identical at 37 N and 52 N — so an empirically fitted response absorbs it
entirely. The same argument disposes of twilight scattering geometry in its simple
form: shadow height above the observer depends on how far the sun is below the
horizon, not on where the observer stands. Consistent with the earlier measurement
that the curve at fixed zenith is stable across seasons to 1-2%.

### What it might be, and why this dataset cannot say

Cloud climatology (the North Pacific storm track puts more cloud north), sea state,
water clarity, or solar-path geometry (the terminator's orientation at fixed z does
differ with latitude). **Not separable here.** One colony and one phenology means the
animals go north IN SUMMER, so latitude and solar geometry are collinear at
r = -0.665: single-predictor R2 is 0.1417 (latitude) against 0.1395 (dz/dt), the joint
model reaches only 0.1442, and the latitude-by-dz/dt cross-tabulation is non-monotone
in both directions. Separating them needs a second colony at a different latitude, or
tags whose excursion and season decouple.

### What this means for the manuscript

This is a **characterised, quantified limitation with a clean control**, not an
unexplained residual bias. The honest statement is that archival light on a
wide-ranging animal is not a pure function of solar zenith over a 15-degree
excursion; the departure is measurable at -0.07 zenith-equivalent degrees per degree
of latitude; and it accounts for the observed latitude bias in full. It generalises
to every light-geolocation study of a wide-ranging animal, which is the audience.

Reproduced end to end by `analysis/diagnostics/diag_latitude_drift.R`.

## 1. The problem

Fitted against Argos truth on 29 double-tagged northern elephant seal
deployments, the grid HMM has a mean **latitude bias of -1.54 degrees** (southward)
at the shipped settings, with interval coverage 0.66 against a nominal 0.95.
Longitude is essentially unbiased (+0.2 to +0.3 deg) and well calibrated (0.99).

A long sequence of candidate explanations was tested and rejected: the movement
model (CRW), the movement scale, Rao-Blackwellisation, the spike normaliser, the
response floor, the shading-arm rate, the composite-likelihood weight, anisotropic
diffusion, response drift, depth attenuation, artificial light, and a declination
sign flip. Each is recorded in `scratch/nes_calibration/SESSION_STATE.md`.

## 2. The null gate — the experiment that changed the diagnosis

`scratch/nes_calibration/build_light_battery.R` builds a battery of **synthetic
light records with exactly known truth**: real Argos movement, real measured
per-observation attenuation, and only the DATE varied (whole-day offsets, which
preserve local time of day so the diel structure of diving stays aligned with the
sun). Six tracks x four offsets = 24 records.

`fit_battery_null.R` then fits that light with **the same logistic response that
generated it**, endpoints pinned to truth. If the estimator cannot recover a known
track from noiseless, correctly-modelled light, then everything downstream of it is
measuring something else.

**It cannot. 48 of 48 fits:**

| level | latitude bias (by track) | rmse | median error | coverage |
|---|---|---|---|---|
| perfect (noiseless, generating response) | **-0.643 deg** | 1.03 | 73 km | 0.999 |
| noisy (clear sky + the engine's own spike-and-slab) | -0.455 deg | 1.92 | 208 km | 0.993 |

All 6 tracks negative (-0.46 to -0.81), p = 0.031 (the floor for n=6 signed-rank).
All 4 date offsets negative (-0.45 to -1.02), so it is **not seasonal**.

**So roughly 40% of the real-data bias is present with perfect, correctly-specified
light.** It is not shading, and it is not response misspecification.

## 3. What was ruled out, and how

| candidate | verdict | evidence |
|---|---|---|
| Argos error in the truth | **excluded by construction** | the battery generates light AT the Argos positions and scores against those same positions, so the error is common to signal and reference and cancels exactly |
| solar zenith wrong | **excluded** | `diag_zenith_reference.R`: an independent NOAA/Meeus implementation agrees with `solar_zenith()` to bias -0.0003 deg, rms 0.003, max 0.009 over 200k samples. Also cancels by construction, since the battery generates with the same function the engine fits with |
| atmospheric refraction | real, but a REAL-DATA issue only | the package computes geometric zenith with no refraction; ours sits +0.221 deg above the sun's apparent position in the twilight band (85-96 deg), up to +0.578 at the horizon. Largely absorbed by the empirically fitted `z50`; the unabsorbed part is the SHAPE distortion across the band |
| readout: posterior mode | **excluded** | `diag_mode_vs_mean.R`: mode -0.639, marginal-lat mode -0.692, median -0.758, mean -0.797. The mean is MORE biased. The ordering is what the measured skewness (-0.145) predicts, so readouts differ only by skew; the posterior itself is displaced |
| grid quantisation | **excluded** | kernel drift identical at cell 1.0 and 0.5 deg to five decimal places |
| the cell-area correction | **exonerated** | removing it moves estimates NORTH (+0.254 deg, 6/6 tracks) — it does not null the drift, it flips the sign |

## 4. A real property of the kernel — but NOT the cause. See section 4a.

> **Read sections 4a and 4b before using anything in this section.** The drift
> described here is real and exactly characterised, but implementing its exact
> correction changed the fitted posterior by ~1e-4 degrees. It is not the mechanism
> behind the bias. Section 4b identifies what is: two per-observation terms in the
> EMISSION. Section 4b also corrects two claims made elsewhere in this note --
> the "monotone in declination" result (an artefact of taking absolute values) and
> the magnitude bookkeeping.

For a point at angular distance `d` and bearing `theta` from `(phi0, lambda0)`:

    sin(phi) = sin(phi0) cos(d) + cos(phi0) sin(d) cos(theta)

Averaging over uniform bearing kills the second term, so `E[sin phi] < sin(phi0)`:
an isotropic kernel on a sphere pulls the mean of sin(latitude) toward zero, i.e.
toward the equator, **before any likelihood is involved**. To second order, with
`u = phi - phi0`,

    u ~= d cos(theta) - (1/2) tan(phi0) d^2 + (1/2) tan(phi0) d^2 cos^2(theta)

and since `E[cos^2 theta] = 1/2` the last term returns half the second, giving

    E[delta phi] ~= -tan(phi0) * sigma^2 / (2 R^2)   per step.

**Measured and confirmed** (`diag_kernel_drift.R`, pure numerics reproducing the
engine's haversine distance, 5-sigma truncation and pre-filters):

- numeric vs analytic correlation **0.99994**
- sigma^2 scaling exact: D440/D110 = **16.05** against a predicted 16
- resolution-independent: cell 1.0 and 0.5 deg agree to 5 dp
- dropping the area factor gives the **mirror image** (35 N, D=110: -0.00298 with,
  +0.00298 without), so the two settings bracket zero symmetrically and **neither
  is drift-free**

Fit-level confirmation (`diag_area_diffusion.R`, 6 tracks, 2x2, perfect light;
reproduction control returned -0.567 exactly):

| arm | bias_mode | rmse | km |
|---|---|---|---|
| area OFF, D=110 | -0.313 | 0.636 | 58 |
| area ON,  D=110 (control) | -0.567 | 0.919 | 68 |
| area OFF, D=440 | -0.934 | 1.433 | 98 |
| area ON,  D=440 | -1.670 | 2.536 | 146 |

D440 - D110 = -1.103 deg, 6/6 tracks, p = 0.031. The ratio is 2.95 with the area
factor and 2.93 without — the same mechanism in both.

Magnitudes line up. At 45 N, D=110: -0.00426 deg/step x 480 knots = -2.04 deg
undamped, against an observed -0.567 (likelihood damping ~3.6x). At D=440: -32.8
undamped against -1.670 (damping ~19.6x). The damping grows as the prior loosens,
which is why a 16x per-step ratio appears as 2.95x at fit level.

**This is correct probability and wrong biology.** The engine faithfully
implements Brownian motion on a sphere, which relaxes to the uniform distribution
on the sphere and therefore genuinely drifts equatorward. As an animal movement
prior it asserts that a seal, absent data, prefers the equator — by about 2 degrees
over a 480-knot track at the shipped `diffusion = 110`.

It bites hardest here because the latent latitude is prior-dominated: the
Gelman-Pardoe pooling factor is `pi_lat = 0.87` against `pi_lon = 0.30`.

### Prior art in this repository

`R/sensor_prior.R` already named the term. Its documentation for `area_prior()`
noted that "isotropic Brownian motion on a sphere also carries an equatorward drift
in its transition (a `tan(lat)` term)", but judged it "far too weak" on SBC
evidence. That judgement was wrong, and for an instructive reason — see below.

## 4a. NEGATIVE RESULT: correcting the drift changes nothing

The Ito correction of section 6 was implemented (`drift_correction`, lib.rs
`ito_lat_shift`) and tested. **It is exact and it does not work.**

- **Exact**: the shift nulls the measured one-step drift to five decimal places
  at every latitude and diffusion tested (-0.00298 -> -0.00000 at 35 N, D = 110).
- **Live**: `log_z` moves (-53051.892 -> -53051.923 at D = 110), so it is applied.
- **Inert on the posterior**: on a 476-knot record at 1 degree with noiseless,
  correctly-modelled light, the posterior MEAN latitude moved by **0.00009 deg** at
  D = 110 and **0.0008 deg** at D = 440 — three to four orders of magnitude below
  the -0.586 and -2.025 it was meant to explain. **Not one of 476 posterior modes
  changed.**

### Why, and why this is the useful part

**This is a forward-backward SMOOTHER, not a filter.** A drift enters the forward
pass going forward in time and the backward pass going backward in time, so at
interior knots it largely cancels in the smoothed marginal. Correcting both passes
consistently cancels in exactly the same way. **A drift that cancels cannot have
been producing the bias in the first place.**

So the mechanism does not act through the transition as a drift. It must act on
the **smoothed marginal** — which is where the one term with a measured fit-level
effect also lives: the area factor multiplies each knot's marginal directly
(`log_cell_area(lat_i)` added to `alpha[k][i]`, once per knot, outside the sum over
sources), rather than entering as a drift that two-sided smoothing can cancel.

### The scaling agreement was a coincidence

The sigma^2 story looked strong because the D = 440 arm was ~3x worse and the
per-step drift scales as sigma^2 (16.05 measured against 16 predicted). But the
real mechanism also widens with the kernel, so that observation is consistent with
both stories and **discriminates neither**. The kernel diagnostic measured a
genuine property of the kernel; it was over-read as the property that reaches the
posterior.

### Disposition

`drift_correction` is kept and **defaults to FALSE** — bit-identical to every
number in this campaign when off, documented so the result is not re-derived, and
defensible on its own terms as the correct prior. Same disposition as
`diffusion_lon`; contrast `overcount`, which was reverted.

## 4b. CORRECTION (2026-08-17): the bias is not monotone in declination, and the
## mechanism is in the EMISSION

Two things asserted earlier in this note are wrong. Both were found by independent
review and then verified directly against `equinox_tilt.csv`.

### The "monotone in |declination|" result was an artefact of the absolute value

Keyed on SIGNED declination, the bias CHANGES SIGN:

| signed declination | mean bias | posterior sd |
|---|---|---|
| -22.3 (NH winter solstice) | -0.902 | 1.00 |
| -12.4 | -1.598 | 1.36 |
| **-7.4** | **-2.156** | 1.65 |
| 0.0 (equinox) | -1.672 | 1.68 |
| +7.5 | -0.311 | 1.26 |
| +12.5 | -0.027 | 1.04 |
| **+22.3 (NH summer solstice)** | **+0.101 POLEWARD** | 0.83 |

At the same |declination|, winter and summer are -2.156 and -0.311; averaging them
manufactured the monotone curve reported earlier. Signed declination fits better
than absolute (R2 0.338 vs 0.268). **The equinox is not the worst case** -- that is
declination -5 to -10.

So the correct description is a tilt pointing toward the **LONGER-DAY side**:
equatorward in northern winter, poleward in northern summer, its effect scaled by
posterior width. Not "the light goes blind at the equinox so a constant equatorward
force wins".

### The magnitude puzzle was my arithmetic

The note previously made much of the fitted slope being 48.9 where "1.0 would be
exact". That calibration is wrong. A PERSISTENT per-knot tilt is amplified by the
smoother's correlation length: bias = c * sum_j Cov(phi_k, phi_j), not c * sigma^2.
Even the pure area tilt would regress with slope ~5-10, never 1.0. The required
per-knot tilt is about **0.1 nats/deg**, not the 0.5-0.9 claimed.

### The mechanism: two per-observation terms in the emission

Measured per-knot tilt at truth on noiseless light (record 2021033/90):

| variant | tilt (nats/deg) |
|---|---|
| full engine emission | -0.116 |
| spike normaliser frozen | -0.063 |
| symmetric arms (`shade_ratio = 1`) | -0.071 |
| **both removed** | **-0.004** |
| pure geometric mismatch | -0.005 |

By zenith band: day -0.092, twilight -0.019, night -0.004 -- it is made in DAYLIGHT,
on the shoulder of the response, not at twilight.

1. **Spike-normaliser gradient** (`spike_normaliser`, lib.rs:65). `Z(mu)` shrinks as
   `mu` approaches the clamp at `max_light`, so cells predicting brighter daytime
   light get a data-independent bonus. Added once PER OBSERVATION.
2. **Arm asymmetry** (`spike_density`, lib.rs:238, `lam_hi = 2*lam_lo`). On noiseless
   data each daytime-shoulder observation charges poleward candidates (which predict
   dimmer than observed) at `2*lambda` and equatorward candidates at `lambda`.

Sufficiency: a 1-D latitude-only forward-backward chain on these exact emission
profiles reproduces the 2-D fit with knot-error correlation **0.974**, and removing
both components plus the area term takes the chain bias from -0.818 to **-0.045**.

### The `perfect` arm is NOT a clean null

Noiseless light sits at the spike's MODE, but `E[score] = 0` holds only for
model-correct data, and the spike's mean is not its mode. The perfect arm therefore
measures the score of one off-model dataset rather than isolating the estimator.
**The honest estimator number is the NOISY arm, -0.455**, where `E[score] = 0` kills
the arm-asymmetry term to first order. Section 2's -0.643 should be read with that
caveat.

### CONFIRMED by three pre-registered predictions (2026-08-17)

Predictions written down before running, chosen because the leading alternatives
predict different outcomes. Posterior-mean latitude bias:

| record | baseline | `shade_ratio = 1` | `lam_mult = 2` | SH mirror |
|---|---|---|---|---|
| 2021033/90 | -1.017 | **-1.291** | **-0.556** | **+0.375** |
| 2023032/0 | -0.692 | **-0.960** | **-0.356** | **+0.733** |
| *predicted* | | *-1.2, WORSE* | *-0.4, better* | *+0.6, SIGN FLIP* |

All three correct in direction on both records, and close in magnitude.

- The **mirror** is the decisive one. Light regenerated from the same response at
  `(-lat, lon)` on the same dates flips the sign of the bias, because the bias is
  EQUATORWARD and equatorward in the southern hemisphere means increasing latitude.
  No apparatus artefact -- grid indexing, interpolation, binning, truth handling --
  predicts a hemisphere flip.
- **`shade_ratio = 1` is the discriminating one.** Symmetrising the arms makes the
  bias WORSE, which this mechanism requires and a "the asymmetry is the culprit"
  reading forbids. Posterior sd behaves as required throughout (1.181 -> 1.365
  symmetrised, -> 0.890 at doubled lambda): since bias ~ tilt * sum_j Cov, a wider
  posterior amplifies the residual normaliser tilt even as the arm term is removed.

### Corroboration already on disk

The `shade_ratio` sweep of 2026-08-06, recorded at the time as a negative result,
is out-of-sample confirmation: ratio 1.00 -> bias -1.88, 2.00 -> -0.88, 4.00 ->
+1.10. Monotone, and *reducing* the ratio to 1 made it WORSE -- which is what this
mechanism predicts and what a "symmetrise the arms" reading does not.

## 5. Why SBC did not catch it

SBC simulates from the model's own prior and likelihood, so it verifies
**self-consistency**. A prior that is internally coherent but wrong about the world
is precisely the failure mode SBC cannot see: the equatorward drift is present in
both the generator and the sampler, so ranks come out uniform. Both grid-HMM arms
PASS (see `notes/topology/sbc_design.md`) while this bias is live.

Real-data interval coverage against independent truth is the complement, and the
two disagree. **Report both.** This is a substantive methodological point for the
Discussion, not a caveat.

## 6a. Next step, now that the mechanism is confirmed (2026-08-17)

**The normaliser cannot simply be removed.** It is what makes the spike a proper
density over `[0, max_light]`, and its absence is exactly what SBC caught before
(see `notes/topology/sbc_design.md`). So this is not a bug fix. The options are:

1. **Change the emission family** so its normaliser does not depend on the expected
   value. Anything whose `Z` is constant in `mu` carries no tilt. This is the
   principled route and it is a MODELLING decision with consequences for the
   manuscript's likelihood section.
2. **Remove the clamp interaction.** The tilt is largest where `mu` approaches
   `max_light` (lib.rs:221) -- it is made in DAYLIGHT (-0.092 nats/deg) not twilight
   (-0.019). A response that does not saturate, or a support that does not truncate
   at `max_light`, would reduce it without changing the family.
3. **Accept and report it.** A correctly-normalised asymmetric emission carries a
   latitude tilt; that is a property of the likelihood, not an error in it. Then the
   honest statement is the size of the resulting bias and its seasonal sign.

Do NOT tune `shade_ratio` to cancel it: symmetrising makes it worse (above), and
the ratio-4 setting that zeroes bias on real data does so by trading one tilt
against another, not by removing either.

**An implication worth pursuing separately.** The tilt is made in DAYLIGHT, on the
shoulder of the response. Much of this campaign -- the `lambda` tuning, the
darkness-regime work -- assumed the action is at twilight because that is where
latitude INFORMATION lives. The information is at twilight; the bias is not.

## 6. Superseded next step (kept for the record)

The Ito correction was the previous next step. It was implemented, tested, and
falsified by its own prediction — see section 4a. **Do not retry it.**

What section 4a establishes is where to look instead: the mechanism survives
two-sided smoothing, so it is **not** a drift in the transition. It acts on the
smoothed marginal. Two candidates, in order:

1. **The area factor as a per-knot marginal TILT, not a drift.**
   `log_cell_area(lat_i)` is added to `alpha[k][i]` once per knot, outside the sum
   over source cells, so it multiplies each knot's marginal by `cos(lat)` directly
   — and unlike a drift, a tilt does not cancel between the forward and backward
   passes. This is the only term with a *measured* fit-level effect (+0.254 deg,
   6/6 tracks, p = 0.031), so it is the one to characterise properly.
   A first-order estimate (`-sigma_post^2 * tan(lat)`, about 0.02 deg at
   `sigma_post = 1.1 deg`) is an order of magnitude too small to explain +0.254,
   so the accounting is not yet right and that gap is the thing to chase. Measure
   it directly on the smoothed marginal rather than deriving it.

2. **How the tilt at OTHER knots reaches knot k through the chain.** The
   single-knot estimate above treats each knot in isolation; the smoother couples
   them, and the coupling strengthens as the movement prior widens — which would
   also produce the observed growth with `diffusion` without any drift.

### Method note, learned the hard way

Both remaining candidates are statements about the **smoothed marginal**, so test
them there. The failure in section 4a came from measuring a property of the
one-step kernel in isolation and assuming it propagated. Any future mechanism
should be demonstrated to survive forward-backward smoothing *before* a correction
is designed for it.

### Scope caution (still applies)

`log_cell_area()` is shared between `run_grid_hmm` (initial, forward and backward
passes) and the block/hierarchical sampler's `emit` path (`log_area_element`,
around lib.rs:1641). Anything done to it hits both engines.

Do not "fix" the tilt by weakening `log_cell_area` to `0.5 * ln(cos lat)`: it would
null the one-step drift, but that drift is not the problem, and it would corrupt
the stationary measure to chase a number.

## 7. What remains unexplained

The null gate accounts for about 40% of the real-data -1.54 deg. After the Ito
correction, the remainder needs its own diagnosis. Current candidates, in order:

1. **Refraction shape distortion** (section 3) — geometry fitted at the haul-out is
   applied at other latitudes and seasons where the refraction-vs-zenith relation
   is sampled differently.
2. **The knot-stationarity approximation** — the engine holds position fixed at one
   cell for a whole 12 h knot while the animal moves continuously. Should be
   near-symmetric in time, so expect ~0.1 deg at most.
3. **Per-tag calibration residuals** — per-tag latitude bias still ranges -7.6 to
   +2.6 deg on real data.

## 8. Consequences for the manuscript

- Section 3.5's account of interval coverage should change. Under-coverage is
  partly a consequence of the shading rate `lambda` being tighter than optimal
  (0.5 beats the shipped 1.0 on accuracy AND coverage: 266 vs 276 km, 0.779 vs
  0.656) and partly this prior drift. It is not an unexplained miscalibration.
- The lambda tuning results should be re-read in this light: lambda has been tuned
  against a fixed -0.57 deg estimator offset, so the setting that zeroes bias is
  partly compensating for a defect rather than describing shading. **The
  lam0.25/lam0.125 arms were deliberately not resumed for this reason.**
- Add the SBC-vs-real-coverage argument (section 5) to the Discussion.
- The movement section should state that the grid HMM's prior is spherical Brownian
  motion, and say what that implies for latitude.
