# Latitude bias in light-based geolocation: current state, for independent review

Supersedes `EVIDENCE_DOSSIER.md`, whose framing was overturned twice.

## The system

`invTwilightFree` reconstructs an animal's track from an archival light logger. A
**grid HMM** over a 1-degree lon/lat lattice, knots every 12 h, exact
forward-backward marginals. Per observation the emission is a spike-and-slab:
an asymmetric Laplace around the expected light (rate `lambda` below the
expectation, `shade_ratio * lambda` above, truncated to `[0, max_light]`) mixed with
a uniform slab of weight `prob_slab`. ~24 observations per knot. Movement prior is
isotropic Gaussian in great-circle distance, `diffusion` km/sqrt(day). A cos(latitude)
cell-area factor is applied per knot. Expected light comes from a fitted response
curve against solar zenith.

Data: 29 northern elephant seal deployments, double-tagged with Argos, ~240 days
each, from Ano Nuevo. Deployment AND recovery positions are always known (an
archival tag must be physically recovered), and for a central-place forager both are
the same colony — a degenerate closure constraint.

## The headline problem

| | honest synthetic | real data |
|---|---|---|
| latitude bias | **-0.51 deg** | **-1.54 deg** |
| median error | 189 km | ~250 km |
| latitude coverage | **0.996** | **0.59** |

The "honest synthetic" is a battery where the emission is correct BY CONSTRUCTION:
clear-sky light plus a draw from the engine's own spike-and-slab (KS-verified), on a
real Argos track. **The coverage failure does not reproduce there at all.** So the
gap is emission fidelity to real light, not inference, geometry, or prior.

## Eliminated, with the method that eliminated each

| candidate | verdict |
|---|---|
| Argos error in truth | cancels by construction in the battery (light generated AT the Argos positions) |
| solar zenith algorithm | agrees with an independent NOAA/Meeus implementation to 0.003 deg rms over 200k samples |
| atmospheric refraction | real (+0.22 deg in the twilight band) but absorbed by the empirically fitted `z50` |
| point-estimate readout | posterior mean is MORE biased than the mode; ordering follows measured skewness |
| grid quantisation | drift identical at 1.0 and 0.5 deg cells |
| cell-area correction | not a bias source — it PREVENTS a +5 deg poleward error |
| spherical movement drift | derived, measured (correlation 0.99994 with theory), corrected exactly, and **inert** (posterior mean moved 9e-5 deg) |
| knot-grid alignment | invariant to 0/3/6/9 h phase shifts (bias moves 0.01 deg) |
| movement prior scale | tightening to measured movement (D=48 vs 110) gives bias -2.95, coverage 0.27 |
| per-axis anisotropic prior | no gain on bias |
| response family | gompertz fits the envelope 6.8x better in twilight and is 3.3x WORSE at fit level |
| calibration parameter uncertainty | marginalising over the z50 posterior widens intervals by **0.014%** |

## What is established about the structure

**The likelihood is sloppy and latitude is the soft direction.** Fisher analysis over
553 one-day windows: condition number 3.8e4; latitude carries weight 0.881 in the
softest eigenvector and 0.005 in the stiffest, trading against the threshold `z50` at
1.4-4.4 deg per degree, **with the exchange rate changing sign across the year**.
This is why axis-aligned hyperparameter sweeps returned nulls for weeks — they were
moves along a level set.

**But the endpoints suppress it ~20-fold.** The FITTED sensitivity to a z50 error is
only 0.07 deg/deg, against the 1.4-4.4 free-offset rate. Pinning both ends removes
most of the estimator's ability to move along the soft direction.

**The calibration IS identifiable from the seasonal contrast.** Gauge inflation falls
1.81 (5-day window) -> 1.08 (240-day) -> 1.01 (6 tags pooled), with `se(z50)`
shrinking as 1/sqrt(N) to 0.045 deg. Profiling `log_z` over `z50` pooled across tags
recovers a known truth exactly. **But the pipeline fits `z50` on a 15-day haul-out
window spanning under 1 degree of declination — where it is worst identified — and
freezes it for a 240-day deployment spanning 47 degrees.** Calibration error costs
~125 km per degree of z50 but only 0.07 deg of latitude.

**The emission's spread is measurably wrong**, measured against the pipeline's own
calibration with clamped bands excluded (the tangent response is pinned at 0 or
`max_light` for 94% of observations; only zenith 80-102 carries graded information):

| zenith band | implied shade_ratio | implied lambda multiple | tail mass vs prob_slab |
|---|---|---|---|
| 80-86 | 2.87 | 6.10x | 0.34 |
| 86-90 | 1.91 | 5.21x | 0.26 |
| 90-94 | 1.87 | 3.91x | 0.09 |
| 94-98 | 1.15 | 7.75x | 0.00 |
| 98-102 | 1.70 | 9.01x | 0.00 |

So: `lambda` is **4-9x too loose**, `prob_slab` (0.10) is **3x or more too heavy**,
and `shade_ratio` (2.0) is about right. The required ratio is NOT monotone in zenith,
so the engine's single darkness threshold cannot express it.

**And there is a direct contradiction.** The residuals say `lambda` should be 4-9x
TIGHTER. The fit-level sweep says the opposite: loosening improves coverage
monotonically (0.421 at lam 4 -> 0.656 at 1 -> 0.779 at 0.5 -> 0.857 at 0.25). The
reconciliation is correlation: thinning to 1-in-4 observations per knot costs no
accuracy, so effective sample size is ~1/4 of nominal, and one `lambda` is carrying
both the noise scale and the independence correction.

**Tempering does not resolve it.** `temper` (the reciprocal of effective sample size)
was predicted to separate the two. At lam 0.25 / temper 0.5: coverage 0.854 (vs 0.857
for lam 0.25 alone) and bias -1.786 (worse than -1.466 baseline). **Temper adds
nothing beyond lambda.** Also: coverage 0.854 is BIMODAL — four tags at 0.97-1.00 and
two at 0.55-0.62 — not a calibrated model.

## The open question

Why is real-data latitude bias -1.5 deg when the same engine is nearly unbiased on
correctly-specified data, and no emission parameter, response family, prior setting
or readout choice removes it?

## Premises that have NOT been examined

1. **That Argos is truth.** A diving seal transmits only at the surface, so Argos
   samples positions non-randomly in time and possibly in space, while the light
   estimate integrates over everything. Never tested.
2. **That the pinned endpoints are benign.** They suppress the soft direction 20x.
   Both pins are at the same colony.
3. **That the 29 deployments are exchangeable.** Coverage is bimodal across tags;
   per-tag latitude bias ranges -7.6 to +2.6 deg.
4. **That the decimation to 30-minute maxima preserves what matters.** Chosen because
   correlation with -zenith rose from 0.43 to 0.81, but never revisited.

## What a useful answer looks like

A DISCRIMINATING experiment: one whose outcome differs between the leading
explanations. Every hypothesis that has failed here fit the existing numbers first,
including matching sign, scaling and interaction terms. Fitting the numbers is not
evidence. A full fit is ~6-8 minutes; pure numerics are free.
