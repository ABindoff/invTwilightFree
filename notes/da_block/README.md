# Surrogate-posterior block / delayed-acceptance sampler (single-track prototype)

A validated R prototype porting the delayed-acceptance / surrogate-posterior
**block** sampler (from the `da_gibbs_geolocation` DSpark-flavoured hierarchical
sampler) onto invTwilightFree's real twilight-free light likelihood, for a
**single track**. This is the first half of Option B: get the block kernel
exact and robust on one track before adding the multi-individual hierarchy.

The exact likelihood is the package's own `eval_logpk_grid()`. Every claim below
is checked against a **gold always-exact single-site sampler** run alongside —
the cheapest correctness check, and the one that surfaced every subtlety here.

## Files

Single track:
- `da_block_singletrack.R` — the sampler (gold single-site + surrogate-posterior
  block), the coarse grid-HMM smoother surrogate (pure R), the batched exact
  emitter, and the gold-vs-block comparison + plot. Driven by a `CFG` list.
- `validate_emit.R` — the M3 gate: proves the pure-R per-obs emitter matches
  `eval_logpk_grid()` to ~1e-14 (needed before batching can be trusted).
- `probe_likelihood.R` — plots the exact per-knot likelihood vs latitude; the
  diagnosis behind M2.

Hierarchy (multi-individual partial pooling):
- `sim_panel.R` — simulate a panel of N individuals with per-individual movement
  scale `σ_i`, forward zenith→light model (`simulate_panel()`; guarded, no
  side effects when sourced).
- `da_hier.R` — the hierarchical Gibbs sampler and the block-vs-gold check on the
  pooled parameters (see the H1/H2 section below). Source with `CFG$lib_only=TRUE`
  to reuse its functions without running (used by the SBC harness).

Calibration (SBC):
- `sbc_stageA.R` — simulation-based calibration of the single-individual movement
  variance `sig2` (block+polish kernel + conjugate update).
- `sbc_stageB.R` — SBC of the full hierarchy (`β` and `sig2_i` with pooling).
  Both reuse the package's `inst/sbc/ecdf_bands.R`.
- `figures/` — the figures referenced below.

## What was learned

### M1 — exactness on well-identified data (`figures/m1_m3_exactness.png`)
On `sim_short` the block reproduces the gold posterior: per-knot means coincide
to 0.02–0.05°, sds within ~20%, acceptance ~0.5. Per-sweep exact work is
identical to gold — the win is **mixing**: the block reaches the same posterior
in ~1/5–1/7th the sweeps, because it moves whole track segments while single-site
crawls.

### M2 — the naive port fails on equinox, and why
(`figures/m2_naive_surrogate_fails.png`, `figures/m2_likelihood_diagnosis.png`)

A naive per-knot moment-matched Gaussian surrogate (the direct analogue of the
reference) **stalls** on `sim_equinox`: cold-started from the northern fix the
block never goes south, while gold follows the true +20°→−40° migration.

Root cause, from probing the likelihood directly: **invTwilightFree latitude is
identified jointly, not per-knot.** A single knot's light peaks far from the
truth (true +18° → peak at −66°); the per-knot surrogate sd blows up to ~37°, the
block proposal degenerates to an uninformative random-walk bridge, and acceptance
collapses. The reference assumed per-node informative observations (`y_t ≈ x_t`);
twilight-free latitude violates that. Broad-enough-to-cover-both-branches is
**necessary but not sufficient** — a single broad Gaussian is too weak to mix.

### M4 — the fix (`figures/m4_coarse_hmm_fix.png`)
Three ingredients, each necessary:
1. **Coarse grid-HMM smoother surrogate** — per-knot `(μ_k, Σ_k)` from a
   forward–backward that fuses accumulated light + movement + the fixed start, so
   the surrogate is correctly *located* (ambiguous-knot sd 42.7° → 6.9°). Built in
   pure R (emissions from `eval_logpk_grid`, own RW transition).
2. **Surrogate-mean initialisation** — escapes a genuine delayed-acceptance
   pathology: from a cold start in the tight surrogate's tail, the
   `sum(ll_exact − ll_cheap)` correction *rejects* moves toward the surrogate
   mode. Init at the smoother mean puts the chain in the right basin.
3. **Random block phase per sweep** — removes residual under-dispersion at the
   steepest knots.

With all three, the equinox block follows the true descent and coincides with
gold (means within 0.26σ, sd ratio 0.81–1.36, acceptance 0.39).

### M3 — batching the exact evaluations (`figures/m1_m3_exactness.png`)
`emit_block_ll()` replicates the spike-and-slab light model in pure R (using the
exported `solar_zenith()` for the zenith) and evaluates all `2·B` exact
log-likelihoods of a block (proposal + current) in **one** `solar_zenith` call.
Validated to 5.7e-14 abs / 2e-16 rel vs `eval_logpk_grid` (`validate_emit.R`), so
exactness is preserved. On `sim_short`: likelihood FFI calls drop 28 → 7.2 per
sweep (3.9×), ~19× fewer in total once the mixing win is included. The systems
payoff scales with how costly the exact evaluation is relative to R glue.

### H1 / H2 — the hierarchy (`figures/h1_panel.png`)
The original goal: partial pooling of the movement scale across individuals,
where the validated per-individual block kernel is reused thousands of times.

- **Model:** `σ_i² ~ InvGamma(a_pop, β)`, `β ~ Gamma(g0, h0)`; Gibbs alternates
  a per-individual track update with conjugate draws of `σ_i²` and `β`.
- **Hyperprior must be vague** (`g0 = h0 = 1e-3`): `σ_i²` is a per-knot km²
  variance (~O(1000)), so an O(1) hyperprior pins the population level wrongly.
- **Block-only under-mixes** the high-frequency track wiggle that drives the
  movement-variance estimate (block proposals are smooth coarse-HMM bridges),
  giving systematically low `σ_i` — a 2.2σ gap vs gold at one individual.
- **Fix = hybrid:** block moves (mix the global track position) + a **batched
  red-black single-site polish** (mix the wiggle), both exact. Red-black =
  interior knots of one parity are conditionally independent given the other, so
  each half-sweep proposes and evaluates all of one parity in one `emit` call
  (2 emit calls per polish, not `K−1`).
- **Result:** block ≡ gold on every pooled parameter (`|Δμ|/σ` 0.01–0.33; pop
  scale 34.7 vs 36.2 km/day). Likelihood FFI calls: gold 1.008M `eval` vs block
  221k (mostly cheap `emit`) — **~4.6× fewer**. Wall-time is *slower* in this toy
  (the likelihood is cheap, so R glue dominates); the call-count reduction is the
  win when the exact evaluation is expensive, which is the real-geolocation case.

Reproduce: `Rscript -e 'CFG <- list(compare=TRUE, sweeps=4000L, gold_sweeps=6000L); source("da_hier.R")'`

### H3 — the pooling benefit (`figures/h3_pooling_benefit.png`)
Hierarchical (pooled) vs independent (unpooled, vague per-track prior) estimates
of `σ_i`, on the **same** panel with the **same** track kernel and initialisation
— only the `σ_i` prior differs (pooled shares the population level `β`;
independent uses a vague per-track `InvGamma`).

On N=12 short (7-day) tracks, where movement scale is weakly identified from
light:
- **Independent** estimates are unstable — 2–49 km/day for truths of 25–71,
  several collapsing near zero (the classic degeneracy of a vague variance prior
  on weak data: a smooth track → tiny increments → `σ→0`, with nothing to stop
  it).
- **Pooled** shrinks each individual toward the *learned* population level (~31),
  because individuals with more signal lift the shared `β` and regularise the rest.
- **RMSE vs truth 27% lower** (27.2 → 19.8); estimate SD 17.5 → 3.2.

Honest nuance: with data this weak, pooling approaches *complete* pooling
(over-shrinks the true heterogeneity). Movement scale is fundamentally hard to
identify from light alone, so neither method recovers it well; pooling's value is
**regularisation and stability** — turning unusable independent estimates into
stable ones. That is precisely the argument for a hierarchical fit on real
short-deployment tags.

Reproduce: `Rscript -e 'CFG <- list(benefit=TRUE, N=12L, days=7, sweeps=4000L); source("da_hier.R")'`

### H4 — partial pooling with heterogeneous track lengths (`figures/h4_partial_pooling.png`)
The richer, textbook regime. A panel of N=12 with track lengths spanning 3–21 days
(`days` may be a per-individual vector in `simulate_panel()`), so per-track
information — and therefore the right amount of shrinkage — varies across
individuals.

- **Short tracks (3–5 days):** independent `σ_i` estimates are wild (truth 47.5 →
  6.5; truth 24.8 → 2.2). Pooling shrinks them hard toward the population →
  **RMSE 63% lower**.
- **Long tracks (14–21 days):** independent estimates are already good (truth
  50.2 → 50.5; truth 70.8 → 68.9 — light *does* identify movement scale given
  enough data). Pooling barely touches them; where it does it slightly
  over-shrinks → **RMSE 34% worse**.
- **Net: 45% lower RMSE**, driven entirely by the short tracks. Mean `|shrinkage|`
  is 21 km/day (short) vs 8 (long) — the shrinkage **decreases monotonically with
  track length** (right panel), which is the defining signature of partial
  pooling: the model borrows strength adaptively, heavily where data is weak and
  hardly at all where it is strong.

This also refutes the earlier "~65% GLS bias" worry: that was a *weak-data*
effect. Long-track independent estimates are near-unbiased, so pooling correctly
leaves them alone.

Reproduce: `Rscript -e 'CFG <- list(benefit=TRUE, N=12L, days=c(3,3,5,5,7,7,10,10,14,14,21,21), sweeps=4000L); source("da_hier.R")'`

### SBC — absolute calibration (`figures/sbc_a_sig2.png`, `sbc_b_beta.png`, `sbc_b_sig2.png`)
The gold check is *relative* (block ≡ gold). It cannot validate two things: (a) the
conjugate `sig2_i`/`β` updates, which are **shared** by both samplers, and (b)
absolute calibration. Simulation-based calibration covers exactly that gap. Both
stages generate from the inference model *exactly* — draw parameters from proper
priors, a knot-level piecewise-constant RW track, and normalised spike-and-slab
light at fixed calibration — then fit and rank the true value among posterior
draws, checked against the Sailynoja et al. simultaneous ECDF band from
`inst/sbc/ecdf_bands.R`.

- **Stage A — single-individual `sig2`** (block+polish kernel + conjugate update):
  **PASS**, 150 replicates, L=100. The rank-ECDF sits inside the 95% band with no
  systematic skew. This confirms the polish fix (§H2) in an absolute sense — the
  movement variance is recovered without bias.
- **Stage B — full hierarchy** (`β` + `sig2_i`, pooling on): **PASS**, 70
  replicates (70 `β` ranks, 280 `sig2_i` ranks). Both parameters calibrate.

A useful incident: Stage B *initially failed* — not a bug (the conjugate updates
are provably correct) but **under-convergence** of the coupled `β`/`sig2` chain
from a cold, over-smooth start. Starting `β` at the prior mean and lengthening
burn-in fixed it. SBC catches convergence problems the gold check (converged by
construction) cannot, which is a second reason to run it before the Rust port.

Caveat: this is the **self-consistency** SBC (generator matches the sampler's
knot-level lon/lat metric), so it validates the *algorithm*. The spherical-metric
robustness probe belongs with the Rust port — mirroring the `sphere` vs
`fit_kernel` split already in `inst/sbc/sbc_grid_hmm.R`.

Reproduce: `Rscript -e 'REPS <- 150L; source("sbc_stageA.R")'` and
`Rscript -e 'REPS_B <- 70L; source("sbc_stageB.R")'` (run from `notes/da_block/`).

### Native Rust port (R1/R2) — wall-time A-B
The R prototype is glue-bound (Cholesky, `ll_cheap` loops, per-call FFI). Two
native kernels in `src/rust/src/lib.rs` port the hot loop (surrogate still built
in R and passed in as `mu_k`, `P_k`):
- `run_block_track` — the single-track block+polish sweep loop (hand-rolled
  Cholesky, DA correction, red-black polish, current-emit caching).
- `run_block_hier` — the whole hierarchy Gibbs natively (per-individual track
  update + conjugate `sig2_i` + conjugate `β`), no per-sweep FFI.

Drivers `ab_singletrack.R`, `ab_hier.R` run the same model both ways and compare.

- **Correctness:** the Rust posteriors coincide with the R prototype — single
  track mean gaps ≤0.06°/0.13σ, sd ratios ~1; hierarchy pop-scale 35.0 vs 35.1
  km/day and per-individual `sig2` `|Δμ|/σ` 0.01–0.14.
- **Wall time (matched sweeps): ~6×.** Single track 53s→8.6s (4× naive, 6× with
  emit caching); hierarchy (N=6, 4000 sweeps) **144s→22.5s**.
- **The ceiling is the likelihood.** The per-obs zenith (`acos`) + spike-slab
  (`exp`, `log`) dominates and is *already compiled* in both (R via `solar_zenith`,
  Rust natively). The port removes the orchestration, FFI, repeated ephemeris, and
  half the likelihood calls (caching) — not the math. Beyond ~6× would require a
  cheaper likelihood (fewer obs, a twilight approximation), orthogonal to the port.

The Rust kernels are internal (no `@export`, like `run_grid_hmm`); the drivers use
`devtools::load_all()`. Calibration transfers by construction: the Rust posteriors
equal the R ones, which pass SBC.

### Spherical movement metric (S1/S2) — `figures/sbc_spherical_metric.png`, `sbc_flat_metric.png`
The prototype used a flat tangent-plane RW at a fixed reference latitude
(`P_move = C⁻¹/sig2`, `C = diag(km_lon², km_lat²)`). That mis-scales longitude far
from the reference and doesn't match the grid engine's great-circle transition.
The Rust kernel now offers `metric = "spherical"` (default in `TwilightFreeHier`),
targeting the engine's model `p(x_k|x_{k-1}) ∝ exp(−gcdist²/2·sig2)`:

- The block **proposal** stays linear-Gaussian (so the tridiagonal Cholesky draw
  is unchanged). The **DA acceptance** gains a movement-prior correction,
  `−[(gc²−planar²)_prop − (gc²−planar²)_cur] / (2·sig2)` summed over each move's
  edges, upgrading the target prior from the linear approximation to great-circle.
- The conjugate `sig2` update uses `SS = Σ gcdist²` (great-circle) instead of the
  planar sum. Exactness is preserved; `metric = "flat"` reproduces the old path.

SBC (self-consistency: generate great-circle tracks, fit, rank `sig2_i`), tuned to
**wide, high-latitude tracks** (deploy −62°, ~11° latitude range) where the metric
matters most:
- **Spherical fit: calibrated** — the rank ECDF wobbles randomly within the band.
- **Flat fit: biased** — a systematic bowl dipping to the band edge (−0.08),
  i.e. `sig2` under-estimated, because the fixed deploy-latitude longitude scaling
  under-scales movement as the track heads equatorward. It scrapes inside the 95%
  band at 60 reps but the *shape* (systematic vs random) is the diagnosis; more
  reps or wider tracks push it out.

Takeaway: spherical is correct and is the right default; the flat metric is a fine
approximation for narrow-latitude tracks (which is why the earlier flat-based work
held up) but degrades for wide-ranging or high-latitude animals.
Reproduce: `Rscript -e 'REPS_S <- 60L; source("sbc_spherical.R")'`.

## The exactness invariant (do not break when extending)

The block correction is valid **only** because the block is drawn from the *same*
surrogate Gaussian that appears in `ll_cheap`: the movement prior cancels and the
Metropolis correction is exactly `sum(ll_exact − ll_cheap)` over the block. A
different surrogate, or a MAP-style draft, invalidates the acceptance ratio.

## Reproduce

```r
# from the package root, with the package installed/loaded
setwd("notes/da_block")

# M3 gate: the batched emitter matches eval_logpk_grid
source("validate_emit.R")            # expect max abs diff ~1e-14

# M1 + M3: exactness and FFI batching on a well-identified track
Rscript -e 'CFG <- list(dataset="sim_short", surrogate="coarse_hmm", coarse_res=1.0,
  diffusion=50, init="surrogate", sweeps=6000L, gold_sweeps=30000L, gold_burn=8000L,
  block_len=5L); OUT_PNG="figures/m1_m3_exactness.png"; source("da_block_singletrack.R")'

# M2: the naive surrogate fails on equinox (block stalls in the north)
Rscript -e 'CFG <- list(dataset="sim_equinox", surrogate="moment_match", diffusion=300,
  init="start", sweeps=8000L, gold_sweeps=20000L, gold_burn=6000L, mesh_lat_pad=30,
  block_len=6L); OUT_PNG="figures/m2_naive_surrogate_fails.png"; source("da_block_singletrack.R")'

# M2 diagnosis: per-knot likelihood vs latitude
source("probe_likelihood.R")

# M4: the coarse-HMM surrogate + surrogate init + random phase recovers equinox
Rscript -e 'CFG <- list(dataset="sim_equinox", surrogate="coarse_hmm", coarse_res=1.5,
  diffusion=300, init="surrogate", sweeps=8000L, gold_sweeps=20000L, gold_burn=6000L,
  mesh_lat_pad=30, block_len=6L); OUT_PNG="figures/m4_coarse_hmm_fix.png"; source("da_block_singletrack.R")'
```

## Next

- Per-tag **SST/bathymetry** demonstrated end-to-end from a real data column via
  `terms = function(id, df) ...` (the mechanism is wired; only a worked example
  with `sst_source()` remains).
- The genuinely-bimodal **equinox** case within the hierarchy (M4 machinery + the
  multimodal fallback), if any panel member crosses the equator near equinox.
- A finer aux grid (currently nearest-cell on the coarse surrogate mesh) if a
  sensor field needs sub-mesh resolution.

Note: the coarse-HMM surrogate is currently built in pure R because the installed
`run_grid_hmm` was a stale binary that segfaulted on fixed points (fixed by
reinstalling from source, commit 206f8ff). On a fresh build the surrogate could
instead use the real `run_grid_hmm` per-knot marginals.
