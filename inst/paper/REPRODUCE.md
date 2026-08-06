# Reproducing the invTwilightFree manuscript results

Every number and figure in the manuscript is produced by a seeded script in the
repository. Run all commands from the package root unless noted, with the
package built (`devtools::load_all(".")`) or installed (`devtools::install()`).
The comparator packages (SGAT, FLightR) and the legacy `TwilightFree` are needed
only for the benchmark tables.

Software used for the reported run: R 4.6.0, Rust 1.95, Windows 11, 24 logical
cores.

## Tables 1 and 2 — accuracy, speed, robustness

Driver: `scratch/run_p8_benchmark.R` (invTwilightFree is refit fresh at seed 42;
SGAT/FLightR are read from the cached fits in `vignettes/cache_*.rds`).

```r
Rscript scratch/run_p8_benchmark.R
```

- Geometric track: simulated in-script (adverse seed 123, ideal seed 42; 180 days,
  10-min cadence). invTwilightFree fit with the guided filter, 1000 particles, seed 42.
- Seal track: `scratch/simulate_seal_track.R` (seed 999) then
  `scratch/generate_seal_light.R` (seed 555) produce
  `scratch/simulated_seal_light_scenarios.rds` (40-day CRW, four light scenarios).
- The two 183-byte SGAT seal caches (`cache_sgat_seal_shaded.rds`,
  `cache_sgat_seal_alan.rds`) are failure placeholders: SGAT could not extract
  twilights under those noise levels, which is why Table 2 shows "—".

To regenerate the comparator caches from scratch (slow: minutes to hours each),
see `scratch/run_benchmark_*.R` and the `vignettes/comprehensive_benchmark.Rmd`
and `vignettes/benchmark_comparison2.Rmd` pipelines.

## Figure 1 — reconstructed tracks across seal scenarios

`vignettes/comprehensive_benchmark.Rmd` (set `NOT_CRAN=true`; requires the seal
scenario RDS and the comparator caches).

## Figure 2 — equinox hemisphere prior

```r
Rscript inst/paper/equinox_hemisphere_prior.R
```

Seed 20240320. Simulates a near-stationary animal at +5 deg latitude over a
five-day window straddling the March 2024 equinox, fits the grid HMM with free
endpoints, light-only versus light + `hemisphere_prior("N")`. Prints the
posterior mode, SD and southern-hemisphere mass and writes
`equinox_hemisphere_prior.png`. Reported: light only SD 3.9 deg / P(south) 0.12;
light + prior SD 0.4 deg / P(south) 0.00.

## Figure 3 — partial-pooling benefit

```r
cd notes/da_block
Rscript -e 'CFG <- list(benefit=TRUE, N=12L, days=c(3,3,5,5,7,7,10,10,14,14,21,21), sweeps=4000L, seed=1L); source("da_hier.R")'
```

Seed 1. Panel of 12 tags, track lengths 3-21 days. Reports RMSE of the recovered
per-tag movement scale, independent versus pooled. Reported: net 45% lower
(21.3 -> 11.6 km/day); short tracks 64% lower (28.6 -> 10.4); long tracks little
changed (9.4 -> 12.8); mean |shrinkage| 21.1 (short) vs 8.6 (long) km/day.

## Figure 4 — simulation-based calibration

- Grid HMM location posterior (panels c, d):
  ```r
  Rscript -e 'source("inst/sbc/sbc_grid_hmm.R")'
  ```
  Seed 1, 100 replicates each of two generators. Writes `sbc_sphere.png`
  (spherical-movement generator: latitude and longitude both PASS the 95%
  simultaneous ECDF band) and `sbc_fitkernel.png` (self-consistency generator:
  longitude PASS, latitude FAIL, the documented un-normalised-spike latitude
  approximation).

- Hierarchical movement variance (panel a), Stage A:
  ```r
  cd notes/da_block
  Rscript -e 'REPS <- 100L; source("sbc_stageA.R")'
  ```
  100 replicates: movement variance PASS. Writes `sbc_stageA_sig2.png`.

- Population scale + per-tag variance (panel b), Stage B:
  ```r
  cd notes/da_block
  Rscript -e 'REPS_B <- 70L; source("sbc_stageB.R")'
  ```
  Both parameters PASS.

- Correlated random walk, persistence and movement variance:
  ```r
  Rscript -e 'REPS <- 100L; source("inst/sbc/sbc_crw.R")'
  ```
  100 replicates of a 3-tag panel, ranks against the 95% simultaneous ECDF
  bands. Both parameters PASS (300 ranks: rho mean rank 0.475, deciles
  0.13/0.11; sigma^2 mean rank 0.504, deciles 0.10/0.11). Writes `sbc_crw.png`
  and `sbc_crw.rds`.

  Two settings matter and are easy to get wrong. The generating prior must
  equal the inference prior, so `rho_max` is passed to the sampler and used as
  the simulation range. And rho mixes slowly because it is coupled to the
  latent track: at 1500 sweeps thinned by 3 the effective sample size is about
  49 of 334 draws, which alone fails the bands. The defaults (6000 sweeps,
  thin 12) give near-independent draws.

- Spherical vs flat movement metric (panel e):
  ```r
  cd notes/da_block
  Rscript -e 'REPS_S <- 60L; source("sbc_spherical.R")'
  ```
  Spherical PASS; flat shows a systematic bias for wide-latitude tracks.

## Real-data application — ten northern elephant seals

```r
rmarkdown::render("analysis/nes_hierarchical_geolocation.Rmd")
```

Renders from any working directory: the setup chunk locates the package root by
walking up for `DESCRIPTION` and sets knitr's root there. Seed 42 throughout.

**Requires `data/nes_untracked/`, which is not in the repository and is not
publicly shareable** (co-author's data, gitignored and `.Rbuildignore`d). Without
it the report stops at the environment check rather than producing anything
misleading. Everything else it needs it fetches or caches: ERDDAP `ncdcOisst21Agg`
and `etopo360` on first run, then `analysis/cache/nes/`.

Seven fitting batches over 2200 tag-days, scored against Argos: hierarchical x
{light, soft fusion, light with the raw clock}, grid HMM x {light, hard fusion},
FFBS x {light, soft fusion} on three tags, plus the naive-versus-pooled
calibration comparison. Several hours on 24 cores; results are cached per batch,
so an interrupted run resumes.

**Cache keys carry a version suffix** (`fit_grid_light_v2` and so on). Bump it
when anything upstream of a fit changes, or a re-run will silently report the
previous numbers under the new text.

For a fast structural check of all 85 chunks, pass `max_tags` and `max_days`:

```r
rmarkdown::render("analysis/nes_hierarchical_geolocation.Rmd",
                  params = list(cache_dir = tempdir(), out_dir = tempdir(),
                                max_tags = 3, max_days = 40))
```

Point `cache_dir` somewhere disposable when doing this: a truncated run writes
fit caches under the same keys as the real one.

## Figure locations

Copies of the generated figures used in the manuscript are collected in
`manuscripts/figures/` (fig2_*, fig3_*, fig4a-e_*).
