# =====================================================================
# SBC Stage A: calibration of the per-individual movement variance sig2.
# Tests the block+polish track kernel together with the conjugate sig2
# update -- the machinery that is SHARED between the gold and block
# samplers and therefore NOT covered by the gold check.
#
# The generator matches the inference model EXACTLY (self-consistency SBC):
#   sig2 ~ InvGamma(a0, b0)                                  [proper prior]
#   knot-level RW: x_1 = deploy (fixed), x_k = x_{k-1} + N(0, sig2 * C^{-1})
#   location is piecewise-constant per knot (as the likelihood assumes)
#   light | x: normalised spike-and-slab at the knot's zenith
# Fit: block sampler (independent, same prior); rank true sig2 among draws.
# If the kernel + conjugate update are calibrated, ranks are uniform.
# =====================================================================
suppressPackageStartupMessages(library(invTwilightFree))

# ---- fixed SBC regime ----
DEPLOY_LON <- 150; DEPLOY_LAT <- -50
DATE   <- "2024-06-15"        # austral winter: latitude well identified
NDAYS  <- 12
STEP_H <- 12
CAL    <- c(558.5, 5.818)     # fixed c(intercept, slope)
LP     <- c(0.03, 64, 0.10)   # fixed c(lambda, max_light, prob_slab)
A0     <- 5; B0 <- 3200        # sig2 ~ InvGamma(A0, B0): mean 800 km^2 (~40 km/day)

# config for the reused hierarchy functions (fixed cal + metric for SBC validity)
CFG <- list(lib_only = TRUE, step_hours = STEP_H,
            calibration = CAL, likpar = LP, lat_ref = DEPLOY_LAT,
            a0 = A0, b0 = B0, coarse_res = 1.0, mesh_pad = 15.0, inflate = 1.5,
            surrogate_diffusion = 45.0, block_len = 5L, polish = TRUE)
source("da_hier.R")           # defines make_individual, update_block, ss_track, ... (no run)
source("../../inst/sbc/ecdf_bands.R")   # reuse the package's SBC ECDF-band machinery

# ---- precompute the fixed time grid, knots, obs->knot bins ----
times  <- seq(as.POSIXct(DATE, tz = "UTC"), by = "10 min", length.out = NDAYS * 24 * 6)
ut     <- as.numeric(times)
t0 <- ut[1]; t1 <- ut[length(ut)]
K  <- ceiling((t1 - t0) / (STEP_H * 3600)) + 1
tstep <- (t1 - t0) / (K - 1)
KNOTS <- t0 + (0:(K - 1)) * tstep
kidx <- findInterval(ut, KNOTS, left.open = TRUE) + 1L   # obs -> knot (t_prev, t_curr]
kidx <- pmin(pmax(kidx, 1L), K)
km_lat <- 111.0; km_lon <- 111.0 * cos(DEPLOY_LAT * pi / 180)

# ---- normalised spike-and-slab sampler (matches the engine density) ----
sample_spike <- function(mu) {
  lam <- LP[1]; maxl <- LP[2]
  m_lo <- 1 - exp(-lam * mu); m_hi <- 0.5 * (1 - exp(-2 * lam * (maxl - mu)))
  if (runif(1) < m_lo / (m_lo + m_hi)) mu - (-log(1 - runif(1) * m_lo) / lam)
  else                                 mu + (-log(1 - runif(1) * (2 * m_hi)) / (2 * lam))
}

# ---- generator: draw sig2 from the prior, a knot RW track, and light ----
gen_replicate <- function() {
  sig2 <- 1 / rgamma(1, shape = A0, rate = B0)         # InvGamma(A0, B0)
  s_lon <- sqrt(sig2) / km_lon; s_lat <- sqrt(sig2) / km_lat
  lon_k <- numeric(K); lat_k <- numeric(K); lon_k[1] <- DEPLOY_LON; lat_k[1] <- DEPLOY_LAT
  for (k in 2:K) { lon_k[k] <- lon_k[k-1] + rnorm(1, 0, s_lon); lat_k[k] <- lat_k[k-1] + rnorm(1, 0, s_lat) }
  # light: piecewise-constant location per knot
  lon_obs <- lon_k[kidx]; lat_obs <- lat_k[kidx]
  z  <- solar_zenith(ut, lon_obs, lat_obs)
  mu <- pmin(pmax(CAL[1] - CAL[2] * z, 0), LP[2])
  light <- vapply(seq_along(ut), function(j)
    if (runif(1) < LP[3]) runif(1, 0, LP[2]) else sample_spike(mu[j]), numeric(1))
  df <- data.frame(time = times, light = light, true_lat = lat_obs, true_lon = lon_obs)
  list(sig2 = sig2, df = df)
}

# ---- single-individual block+polish fit -> sig2 posterior draws ----
CTR$calls <- 0L; CTR$emit <- 0L
fit_sig2_draws <- function(E, sweeps, burn, thin) {
  x <- E$init(); sig2 <- ss_track(E, x) / (E$K - 1)
  out <- numeric(0)
  for (s in seq_len(sweeps)) {
    x <- update_block(E, x, sig2, cfg$block_len)
    SS <- ss_track(E, x)
    sig2 <- 1 / rgamma(1, shape = A0 + (E$K - 1), rate = B0 + 0.5 * SS)
    if (s > burn && (s - burn) %% thin == 0) out <- c(out, sig2)
  }
  out
}

# ---- harness ----
SBC_REPS <- if (exists("REPS")) REPS else 20L
SWEEPS <- 3000L; BURN <- 600L; THIN <- 24L    # L = (3000-600)/24 = 100 draws
ranks <- integer(0); Ldraws <- integer(0); t_start <- Sys.time()
for (r in seq_len(SBC_REPS)) {
  g <- gen_replicate()
  E <- make_individual(g$df, DEPLOY_LON, DEPLOY_LAT)
  draws <- fit_sig2_draws(E, SWEEPS, BURN, THIN)
  ranks <- c(ranks, sum(draws < g$sig2)); Ldraws <- c(Ldraws, length(draws))
  if (r %% 10L == 0L) cat(sprintf("  SBC %d/%d (%.0fs elapsed)\n", r, SBC_REPS,
                                  as.numeric(Sys.time() - t_start, units = "secs")))
}
L <- Ldraws[1]
cat(sprintf("\nStage A SBC: %d replicates, L=%d draws each\n", length(ranks), L))

df_ranks <- data.frame(rank_sig2 = ranks)
sbc_ecdf_report(df_ranks, c(rank_sig2 = "sigma^2 (movement variance)"), L,
                conf = 0.95, fig = "sbc_stageA_sig2.png", tag = "block+polish | sig2")
saveRDS(list(ranks = ranks, L = L), "sbc_stageA_ranks.rds")
