# Simulation-based calibration for invTwilightFree: design and caveats

Scaffold: `inst/sbc/sbc_grid_hmm.R` (+ vendored `inst/sbc/ecdf_bands.R`).
Method: Talts et al. 2018; simultaneous ECDF bands of Sailynoja, Burkner &
Vehtari 2022. Goal: show the per-knot location posteriors are calibrated (the
credible regions cover at the nominal rate), which is the honest-uncertainty
contribution flagged as the highest-value remaining work.

## What SBC tests here

Across many replicate tracks drawn from the model's own prior-predictive, the
rank of the true location among posterior draws must be uniform. A flat rank
ECDF (hugging the diagonal) means calibrated. A U-shape means the posteriors are
too narrow (overconfident, the dangerous direction); an n-shape means too wide
(conservative); a slope means directional bias.

## Design decisions (and why)

1. Fixed vs estimated. The grid HMM estimates only the latent track; calibration,
   diffusion and the likelihood parameters are fixed inputs. Per the SBC cardinal
   rule, a parameter must be drawn from its prior in generation only if the fit
   estimates it. So these are FIXED in both the generator and the fit (SBC option
   2), and need no prior draw. The tested "parameter" is the latent location.

2. The track prior is the movement model. We draw the true track from the EXACT
   transition kernel the HMM uses: a 2D isotropic Gaussian step with
   `sd = diffusion * sqrt(dt_days)` moved along a great circle. If the generator
   kernel did not match the fit kernel, SBC would fail for that reason alone, so
   `sphere_step()` mirrors the Rust move exactly.

3. Known endpoints. Both deploy and retrieve are fixed in the fit (their true
   values), matching real use where both are known. The interior is inferred.

4. One independent rank per replicate. Latent states within a track are
   correlated, so pooling ranks across knots within a replicate is not
   independent and would mis-calibrate the ECDF band. We therefore rank the true
   latitude and longitude at the MIDPOINT knot only (farthest from both fixed
   ends, widest posterior), giving one independent rank per replicate per
   coordinate. Testing additional fixed knot indices separately is fine; pooling
   across knots within a track is not.

5. Discrete posterior -> ranks. The grid posterior is a categorical over cells
   (exposed via `grid_posterior()`; the new `posterior` field on `run_grid_hmm`).
   We draw L cells from it and add uniform within-cell jitter (+/- cell/2) so the
   draws are continuous and grid ties do not bias the rank.

6. No convergence filter needed. The grid HMM forward-backward is exact and
   deterministic, so there are no divergences and the L draws are iid from the
   discrete posterior (ESS = L). SBC here tests the implementation and the grid
   resolution jointly; use a fine enough grid that discretisation is not the
   limiting factor, or a coarse-grid failure will look like miscalibration.

## The real caveat: the spike is not normalised

The engine's spike density `lambda*exp(-lambda*(mu-L))` (and the doubled-rate
branch above `mu`) is not normalised over the clamped light range `[0, max_light]`,
and that normaliser depends weakly on `mu` (hence on the cell), mainly for cells
where the sun is fully up or down. The HMM posterior is therefore not exactly the
posterior of a proper probability model.

For SBC we generate from the PROPER normalised mixture (a normalised spike plus
the uniform slab). Consequences:

- If SBC PASSES, the un-normalisation is immaterial at this grid/regime: good news
  and a clean calibration claim.
- If SBC shows a systematic U or n shape, part of it may be the spike
  normalisation mismatch rather than a forward-backward bug. The clean fix is to
  divide the spike by its per-`mu` normaliser in the engine, making the model
  proper; that is a small, well-scoped change worth doing if SBC flags it.

Localise before concluding (per the skill): re-run on the convergence-trivial set
(all of them here), stratify by regime (winter vs equinox), and if a coordinate
fails, check whether fixing the suspect (e.g. running a finer grid, or the
normalised spike) restores calibration.

## Regimes to report

The scaffold runs an austral-winter window where latitude is well identified.
Report that regime explicitly. Then run two more as separate tests: an equinox
window (the hemisphere-degenerate case, where you expect the bimodal posterior
and should check the class-aware path holds coverage) and a higher-diffusion
case (wider movement). SBC validates the algorithm under the priors and regime
you ran, not the package defaults in general; say so.

## First pilot result (100 reps, winter regime)

Longitude: calibrated (ECDF inside the band, hugging the diagonal).

Latitude: a one-signed positive hump peaking ~+0.2 near mid-rank, returning to 0
at both ends. That is a DIRECTIONAL BIAS (not a width problem, which would be
antisymmetric about 0.5), with an excess of low ranks, i.e. the true latitude
sits systematically low in the posterior. It is latitude-specific (longitude is
clean), which rules out a generic harness bug.

Leading explanation: the grid HMM transition kernel is a planar Gaussian over
lon/lat cells with a `-2 ln sigma` normaliser and a count-uniform initial prior;
it omits the spherical area element `cos(latitude)`. Over this grid `cos(lat)`
ranges ~0.41 to 0.83 (a factor of 2), so cells do not carry their true area
weight, which biases latitude and is neutral in longitude, matching the figure.

Diagnostic (now in the runner): run BOTH generators.
- "sphere" (a real, area-respecting animal) vs "fit_kernel" (the HMM's own
  discrete kernel, same planar approximation as the fit).
- sphere FAIL(lat) + fit_kernel PASS(lat)  =>  the planar/area approximation is
  the cause. Fix: weight cells by `cos(lat)` in the transition kernel and the
  initial prior (the spherical area measure). Longitude unaffected.
- both FAIL(lat)  =>  deeper bug; prime suspect is the un-normalised spike
  emission (its per-mu normaliser depends on latitude). Fix: normalise the spike.

This is the SBC payoff: a real latitude bias that point-accuracy benchmarks did
not surface, with a concrete, well-scoped fix once the generator contrast says
which mechanism it is.

## Diagnosis and fix (both generators failed latitude)

Both generators failed latitude; longitude passed in both. The fit-kernel run
matches the movement model to the engine exactly, so its only remaining
generator/fit difference is the emission. It still failed latitude, which pins
the cause on the emission: the un-normalised spike. Its normaliser `N(mu)`
depends on the expected light `mu`, hence on latitude through the solar geometry,
so the missing constant tilted the per-cell likelihood by latitude. The sphere
run failed with an extra one-signed (directional) component, indicating the
missing `cos(lat)` area term contributes a smaller, secondary latitude bias.

Fix applied: the engine now normalises the spike by
`N(mu) = (1 - e^{-lambda*mu}) + 0.5*(1 - e^{-2*lambda*(max_light - mu)})` at all
five call sites (`spike_density` takes `max_light`). The observation model is now
proper; `test-light_likelihood_norm.R` checks it integrates to 1. The SBC
generator already draws from the proper normalised mixture, so fit-kernel SBC
should now PASS. Re-run:

1. fit-kernel SBC -> expect latitude PASS (emission fixed; gen and fit agree).
2. sphere SBC -> if latitude still shows a directional bias, add `cos(lat)` area
   weighting to the transition kernel and initial prior, then re-run.

Two consequences to carry forward: (a) the likelihood changed, so the benchmark
numbers and any cached fits must be regenerated; (b) the manuscript Eq. 2.2 must
now include the normaliser `N(mu)` (the spike is `f_spike / N(mu)`).

## After the spike fix: fit-kernel PASS, sphere still fails latitude

Re-run result: fit-kernel SBC passes both coordinates (generator and engine now
share the same proper, normalised emission and the same movement kernel), which
validates the spike normalisation and the forward-backward implementation. The
sphere run still fails latitude, now with a clean antisymmetric shape (below the
band at low ranks, above at high ranks): the "too wide" / over-cover signature.

Cause: the grid kernel treats lon/lat cells as equally weighted and omits the
`cos(lat)` area element, so under a realistic area-respecting animal the latitude
prior is too flat and the posterior over-covers. Fix (no Rust change): add
`log(cos(lat))` per cell at every knot via the existing terms/aux mechanism,
which threads consistently through the forward and backward passes because
`logpk` feeds both. Implemented as `area_prior()` (pair with `identity_rule()`),
now included in the SBC fit.

With the area-weighted fit, the canonical test is the SPHERE generator (a real
animal); fit-kernel now intentionally mismatches (its generator omits area
weighting) and is expected to fail, retained only as the emission control.

If sphere latitude calibrates with `area_prior()`, the recommendation is to make
area weighting the default in `TwilightFreeGrid` (it is the physically correct
measure), folding it in with the benchmark regeneration that the spike fix
already requires.

## Correction: the residual is the planar-vs-spherical MOVEMENT model, not area

The `area_prior()` (log(cos lat) emission term) was applied correctly (a single
fit's posterior changed) but had negligible effect, and did not fix the sphere
SBC. Reason: an emission-side per-cell prior is only the STATIONARY part of the
spherical area measure. Isotropic Brownian motion on a sphere also has an
equatorward drift in its TRANSITION (a `tan(lat)` term; in colatitude the
generator carries `cot(theta)`), which the planar Gaussian kernel omits and which
a per-cell prior cannot reproduce. That drift, not a width issue, is the
gen/fit mismatch behind the sphere latitude failure.

Resolution and reporting:
- The calibration claim is the FIT-KERNEL SBC (data from the fitted model): it
  passes both coordinates. The spike normalisation was necessary and sufficient
  for that. This is standard SBC and what goes in the paper.
- The SPHERE SBC is a model-misspecification robustness probe. The planar
  movement kernel over-covers latitude under fully spherical movement at high
  latitude / long tracks. Report as a known limitation. A spherically-exact
  transition (area-weighted with per-source renormalisation, giving the correct
  drift) is future work; it is a Rust change to the forward/backward transition,
  not an emission term.
- `area_prior()` is retained as a weak optional prior with corrected docs; it is
  NOT a spherical-movement fix. It was removed from the SBC fit.

Net: SBC validation is in hand. The spike normalisation is a real correctness
fix (carry it into the benchmarks and Eq. 2.2). The cos(lat) area weighting line
in the manuscript movement section should be dropped or recast as future work.

## Status

Pilot at `SBC_REPS = 100`; scale to 300-1000 for the reported result. Requires
the package rebuilt with the `posterior` field (a `devtools::load_all()` /
reinstall), then `devtools::document()` for the new `grid_posterior()` export.
