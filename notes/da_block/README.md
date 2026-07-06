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
  pooled parameters (see the H1/H2 section below).
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

- Richer **partial-pooling** regime (longer or heterogeneous track lengths), where
  well-identified individuals stay near their data while noisy ones shrink — a
  more nuanced picture than the near-complete pooling of the short-track panel.
- The genuinely-bimodal **equinox** case within the hierarchy (M4 machinery + the
  multimodal fallback), if any panel member crosses the equator near equinox.
- Port the validated kernel to the compiled engine if the wall-clock (not just
  FFI-call) win is needed at scale.

Note: the coarse-HMM surrogate is currently built in pure R because the installed
`run_grid_hmm` was a stale binary that segfaulted on fixed points (fixed by
reinstalling from source, commit 206f8ff). On a fresh build the surrogate could
instead use the real `run_grid_hmm` per-knot marginals.
