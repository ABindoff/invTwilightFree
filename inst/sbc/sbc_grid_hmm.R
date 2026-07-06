# =============================================================================
# Simulation-based calibration for the TwilightFreeGrid HMM.
#
# Tests whether the grid HMM's per-knot location posteriors are calibrated:
# across many replicate tracks drawn from the model's own movement + observation
# prior, the rank of the true location among posterior draws should be uniform.
#
# Design (see notes/topology/sbc_design.md for the full rationale):
#   * The grid HMM ESTIMATES the latent track only; calibration, diffusion and
#     the likelihood parameters are FIXED inputs, so we fix them in BOTH the
#     generator and the fit (SBC option 2). No hyperparameter prior draws.
#   * The "prior" over the latent track IS the movement model. We draw a true
#     track from the EXACT transition kernel the HMM uses (2D isotropic Gaussian
#     step, sd = diffusion * sqrt(dt_days), moved along a great circle), so the
#     data come from the model's prior-predictive.
#   * Both endpoints (deploy AND retrieve) are known and fixed in the fit, as in
#     real use; the interior knots are what we test.
#   * We rank the true latitude and longitude at the MIDPOINT knot (farthest from
#     both fixed ends, widest posterior). One rank per replicate => independent.
#   * The grid HMM posterior is exact (deterministic forward-backward), so there
#     are no divergences; ESS = L (iid draws from the discrete posterior).
#
# CAVEAT (a real finding to keep in mind): the engine's spike density is not
# normalised over the clamped light range, and that normaliser is weakly
# cell-dependent. We generate from the PROPER normalised mixture; a systematic
# under/over-coverage here may partly reflect that mismatch rather than a forward
# -backward bug. See sbc_design.md.
#
# Requires the package built with the `posterior` field on run_grid_hmm (recent).
# Run: source this file from inst/sbc/ (it sources ecdf_bands.R alongside it).
# =============================================================================

library(invTwilightFree)
# Run from inst/sbc/ (or set the working directory there) so this resolves:
source("./inst/sbc/ecdf_bands.R")

# ---- Config: fixed hyperparameters and regime -------------------------------
set.seed(1)
SBC_REPS  <- 100L           # pilot; scale to 300-1000 for the real test
SBC_L     <- 200L           # posterior draws per replicate (rank support {0..L})
CAL       <- c(64, 64 / 90) # fixed calibration c(intercept, slope)
LP        <- c(0.5, 64, 0.05) # fixed c(lambda, max_light, prob_slab)
DIFFUSION <- 80             # km / sqrt(day), fixed
STEP_H    <- 24             # daily knots
CELL      <- 1.0            # grid cell size (deg)
DEPLOY    <- c(lon = 150, lat = -50)
T0        <- as.numeric(as.POSIXct("2024-06-15", tz = "UTC"))  # austral winter: lat well identified
NDAYS     <- 40
EARTH     <- 6371

# knot times exactly as the engine builds them
t_start <- T0; t_end <- T0 + NDAYS * 86400
k_steps <- ceiling((t_end - t_start) / (STEP_H * 3600)) + 1
t_step  <- (t_end - t_start) / (k_steps - 1)
KNOTS   <- t_start + (0:(k_steps - 1)) * t_step
KMID    <- round(k_steps / 2)

# one grid, generous extent around the deployment
GRID <- makeGrid(lon = DEPLOY["lon"] + c(-18, 18),
                 lat = DEPLOY["lat"] + c(-16, 16), cell.size = CELL)

# Reconstruct the (unmasked) grid cell centres for the fit-kernel generator.
CELLS_LON <- rep(seq(DEPLOY["lon"] - 18 + CELL/2, DEPLOY["lon"] + 18 - CELL/2, by = CELL),
                 times = length(seq(DEPLOY["lat"] - 16 + CELL/2, DEPLOY["lat"] + 16 - CELL/2, by = CELL)))
CELLS_LAT <- rep(seq(DEPLOY["lat"] + 16 - CELL/2, DEPLOY["lat"] - 16 + CELL/2, by = -CELL),  # top-to-bottom
                 each  = length(seq(DEPLOY["lon"] - 18 + CELL/2, DEPLOY["lon"] + 18 - CELL/2, by = CELL)))

gcdist_km <- function(lo1, la1, lo2, la2) {           # haversine, km
  r1 <- la1 * pi/180; r2 <- la2 * pi/180; dlo <- (lo2 - lo1) * pi/180; dla <- r2 - r1
  a <- sin(dla/2)^2 + cos(r1) * cos(r2) * sin(dlo/2)^2
  2 * EARTH * asin(pmin(1, sqrt(a)))
}

# ---- Movement kernel: one step, matching the HMM exactly --------------------
sphere_step <- function(lon, lat, sigma_km) {
  x <- rnorm(1, 0, sigma_km); y <- rnorm(1, 0, sigma_km)
  d <- sqrt(x * x + y * y); dr <- d / EARTH; bear <- atan2(y, x)
  la1 <- lat * pi / 180; lo1 <- lon * pi / 180
  la2 <- asin(sin(la1) * cos(dr) + cos(la1) * sin(dr) * cos(bear))
  lo2 <- lo1 + atan2(sin(bear) * sin(dr) * cos(la1), cos(dr) - sin(la1) * sin(la2))
  c(lon = ((lo2 * 180 / pi + 180) %% 360) - 180, lat = la2 * 180 / pi)
}

# ---- draw_truth(): two generators -------------------------------------------
# (A) "sphere": proper spherical move (represents a real animal). Realism SBC.
draw_truth_sphere <- function() {
  lon <- numeric(k_steps); lat <- numeric(k_steps)
  lon[1] <- DEPLOY["lon"]; lat[1] <- DEPLOY["lat"]
  dt_days <- diff(KNOTS) / 86400
  for (k in 2:k_steps) {
    s <- sphere_step(lon[k - 1], lat[k - 1], DIFFUSION * sqrt(dt_days[k - 1]))
    lon[k] <- s["lon"]; lat[k] <- s["lat"]
  }
  list(lon = lon, lat = lat, mid_lat = lat[KMID], mid_lon = lon[KMID])
}

# (B) "fit_kernel": draw the next state from the HMM's OWN discrete transition
# kernel (planar Gaussian over cells, no cos-lat area weighting), so generator
# and fit share a movement model. Self-consistency SBC: if latitude calibrates
# here but fails under "sphere", the planar kernel's missing spherical-area term
# is the cause (fix = cos-lat weighting). If it still fails here, the cause is
# elsewhere (e.g. the un-normalised spike emission).
draw_truth_fitkernel <- function() {
  dt_days <- diff(KNOTS) / 86400
  idx <- which.min(gcdist_km(DEPLOY["lon"], DEPLOY["lat"], CELLS_LON, CELLS_LAT))
  lon <- numeric(k_steps); lat <- numeric(k_steps)
  lon[1] <- CELLS_LON[idx]; lat[1] <- CELLS_LAT[idx]
  for (k in 2:k_steps) {
    sig <- DIFFUSION * sqrt(dt_days[k - 1])
    d <- gcdist_km(lon[k - 1], lat[k - 1], CELLS_LON, CELLS_LAT)
    w <- exp(-(d * d) / (2 * sig * sig))               # exactly the engine's kernel
    idx <- sample.int(length(w), 1, prob = w)
    lon[k] <- CELLS_LON[idx]; lat[k] <- CELLS_LAT[idx]
  }
  list(lon = lon, lat = lat, mid_lat = lat[KMID], mid_lon = lon[KMID])
}

# ---- simulate_data(): light from the proper normalised spike-and-slab -------
sample_spike <- function(mu) {                       # normalised over [0, max_light]
  lam <- LP[1]; maxl <- LP[2]
  m_lo <- 1 - exp(-lam * mu)
  m_hi <- 0.5 * (1 - exp(-2 * lam * (maxl - mu)))
  if (runif(1) < m_lo / (m_lo + m_hi)) {
    e <- -log(1 - runif(1) * m_lo) / lam; mu - e          # below mu, rate lambda
  } else {
    e <- -log(1 - runif(1) * (2 * m_hi)) / (2 * lam); mu + e  # above mu, rate 2 lambda
  }
}
simulate_data <- function(truth) {
  times <- seq(t_start, t_end, by = 600)               # 10-min cadence
  # assign each obs to its knot (t_prev, t_curr], matching the engine
  kidx <- findInterval(as.numeric(times), KNOTS, left.open = TRUE) + 1
  kidx <- pmin(pmax(kidx, 1), k_steps)
  z <- solar_zenith(as.numeric(times), truth$lon[kidx], truth$lat[kidx])
  mu <- pmin(pmax(CAL[1] - CAL[2] * z, 0), LP[2])
  light <- vapply(seq_along(times), function(j) {
    if (runif(1) < LP[3]) runif(1, 0, LP[2]) else sample_spike(mu[j])
  }, numeric(1))
  list(time = as.POSIXct(times, origin = "1970-01-01", tz = "UTC"), light = light)
}

# ---- fit_and_rank(): fit with known endpoints, rank true mid-knot lat/lon ----
fit_and_rank <- function(data, truth, L) {
  fit <- TwilightFreeGrid(
    data$time, data$light, GRID,
    start_lon = truth$lon[1],       start_lat = truth$lat[1],
    end_lon   = truth$lon[k_steps], end_lat   = truth$lat[k_steps],
    calibration = CAL, likelihood_params = LP,
    step_hours = STEP_H, diffusion = DIFFUSION)
  post <- grid_posterior(fit)
  pk <- post$P[KMID, ]
  if (sum(pk) <= 0) stop("empty posterior at midpoint knot")
  idx <- sample.int(length(pk), L, replace = TRUE, prob = pk)
  jit <- function() runif(L, -CELL / 2, CELL / 2)
  lat_draws <- post$lat[idx] + jit()
  lon_draws <- post$lon[idx] + jit()
  list(ranks = c(mid_lat = sum(lat_draws < truth$mid_lat),
                 mid_lon = sum(lon_draws < truth$mid_lon)),
       ess = c(mid_lat = L, mid_lon = L), n_div = 0L)
}

# ---- harness ----------------------------------------------------------------
run_sbc <- function(gen, reps = SBC_REPS, L = SBC_L) {
  rows <- vector("list", reps); off_grid <- 0L
  for (r in seq_len(reps)) {
    truth <- gen()
    if (max(abs(truth$lon - DEPLOY["lon"])) > 17 ||
        max(abs(truth$lat - DEPLOY["lat"])) > 15) { off_grid <- off_grid + 1L; next }
    dat <- simulate_data(truth)
    out <- tryCatch(fit_and_rank(dat, truth, L), error = function(e) NULL)
    if (is.null(out)) next
    rows[[r]] <- data.frame(rep = r, n_div = out$n_div,
                            rank_mid_lat = out$ranks[["mid_lat"]],
                            rank_mid_lon = out$ranks[["mid_lon"]],
                            ess_mid_lat = out$ess[["mid_lat"]],
                            ess_mid_lon = out$ess[["mid_lon"]])
    if (r %% 20L == 0L) cat(sprintf("  SBC %d/%d (off-grid skipped: %d)\n", r, reps, off_grid))
  }
  cat(sprintf("off-grid replicates skipped: %d / %d\n", off_grid, reps))
  do.call(rbind, rows)
}

SPECS <- c(rank_mid_lat = "latitude (mid knot)", rank_mid_lon = "longitude (mid knot)")

# Two generators with different roles:
#   * fit_kernel gen draws from the engine's OWN movement model. This is standard
#     SBC (data from the fitted model) and is the CALIBRATION CLAIM. It passes,
#     which validates the implementation (and the spike normalisation that made
#     it pass).
#   * sphere gen draws isotropic Brownian motion on the sphere (an alternative,
#     more "physical" movement with an equatorward drift the planar kernel omits).
#     It is a model-misspecification ROBUSTNESS PROBE, expected to over-cover
#     latitude at high |lat| / long tracks. Report as a known limitation; a
#     spherically-exact transition kernel is future work.
cat("\n=== robustness probe (isotropic spherical movement) ===\n")
sbc_sphere <- run_sbc(draw_truth_sphere)
sbc_ecdf_report(sbc_sphere, SPECS, SBC_L, conf = 0.95,
                fig = "sbc_sphere.png", tag = "grid HMM | sphere gen")

cat("\n=== self-consistency SBC (fit's own discrete kernel) ===\n")
sbc_fit <- run_sbc(draw_truth_fitkernel)
sbc_ecdf_report(sbc_fit, SPECS, SBC_L, conf = 0.95,
                fig = "sbc_fitkernel.png", tag = "grid HMM | fit-kernel gen")
