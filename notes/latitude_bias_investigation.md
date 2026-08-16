# The residual latitude bias: what it is, what it is not, and what to do next

Status as of 2026-08-16. Branch `fix/light-response-autodetect`.

This note exists so that the manuscript can be edited from a written record rather
than from recollection. Everything here is reproducible from scripts in
`scratch/nes_calibration/` (gitignored, but present on the working machine) and
from the numbers quoted below.

---

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

## 4. The cause: an isotropic movement prior on a sphere drifts equatorward

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

## 5. Why SBC did not catch it

SBC simulates from the model's own prior and likelihood, so it verifies
**self-consistency**. A prior that is internally coherent but wrong about the world
is precisely the failure mode SBC cannot see: the equatorward drift is present in
both the generator and the sampler, so ranks come out uniform. Both grid-HMM arms
PASS (see `notes/topology/sbc_design.md`) while this bias is live.

Real-data interval coverage against independent truth is the complement, and the
two disagree. **Report both.** This is a substantive methodological point for the
Discussion, not a caveat.

## 6. Next step

Implement an **Ito correction**: add `+tan(phi) * sigma^2 / (2 R^2)` to the
movement kernel so that latitude is a martingale under the prior. This removes the
coordinate drift without touching the emission, the area measure, or the movement
scale.

The falsifiable prediction, testable with the existing harness unchanged:

1. the null gate's -0.567 deg (area ON, D=110) should fall toward zero;
2. **more diagnostically**, the D=440 arm should stop being ~3x worse than D=110,
   because the σ²-scaled drift is what makes a looser prior worse.

If (1) improves but (2) does not, the correction is being absorbed rather than
fixing the mechanism.

### Scope caution

`log_cell_area()` is shared. It is used by `run_grid_hmm` (initial, forward and
backward passes) and the block/hierarchical sampler's `emit` path
(`log_area_element` around lib.rs:1641). A drift correction belongs in the
TRANSITION, not in the stationary area measure, so the change must be scoped to the
movement kernel deliberately rather than applied by editing `log_cell_area`.

A one-line empirical shadow of the correction — `log_cell_area` returning
`0.5 * ln(cos lat)`, i.e. halfway between the two bracketing settings — would also
null the drift, but it is a coincidence of the algebra rather than a defensible
model, and it would corrupt the stationary measure. Do not ship it.

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
