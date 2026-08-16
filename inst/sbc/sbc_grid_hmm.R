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
# RESOLVED (was: "the engine's spike density is not normalised over the clamped
# light range"). It is normalised -- see `spike_normaliser` in src/rust/src/lib.rs,
# whose form matches `sample_spike()` below exactly, and which was added because
# this script caught the mu-dependence. Generator and fit agree on the emission,
# so a failure here is not attributable to the spike support.
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
AREA      <- TRUE           # cell-area weighting; MUST match TwilightFreeGrid's
                            # default. The fit-kernel generator below reads it, so
                            # generator and fit cannot drift apart again.
DEPLOY    <- c(lon = 150, lat = -50)
T0        <- as.numeric(as.POSIXct("2024-06-15", tz = "UTC"))  # austral winter: lat well identified
NDAYS     <- 40

# DESIGN CAVEAT, and it matters more than it looks. Forty DAILY knots, both
# endpoints anchored, at 50 degrees south in strong austral winter light, is a
# regime where latitude is sharply identified and the movement prior barely moves
# the posterior. A term can be worth nothing here and a great deal on a 12-hourly
# track of several hundred knots with a free retrieval endpoint near an equinox.
# That is not hypothetical: the cell-area term was retired on the strength of a
# null result from this design and later found to be worth 0.4 to 2.2 degrees of
# latitude bias on real tracks. Do NOT generalise a negative from this script
# without re-running it at the knot count, endpoint freedom and season of the
# application you mean to generalise to.
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
# kernel, so generator and fit share a movement model by construction. It reads
# AREA, so it follows the engine's cell-area setting rather than assuming one --
# an earlier version hard-coded "no cos-lat weighting" in a comment and in the
# code, and silently became a generator/fit MISMATCH the day the engine's default
# flipped to TRUE. Self-consistency SBC: a failure here is in the implementation
# (e.g. the un-normalised spike emission), not in the movement model, because
# there is no movement model to be wrong about.
draw_truth_fitkernel <- function() {
  dt_days <- diff(KNOTS) / 86400
  idx <- which.min(gcdist_km(DEPLOY["lon"], DEPLOY["lat"], CELLS_LON, CELLS_LAT))
  lon <- numeric(k_steps); lat <- numeric(k_steps)
  lon[1] <- CELLS_LON[idx]; lat[1] <- CELLS_LAT[idx]
  for (k in 2:k_steps) {
    sig <- DIFFUSION * sqrt(dt_days[k - 1])
    d <- gcdist_km(lon[k - 1], lat[k - 1], CELLS_LON, CELLS_LAT)
    w <- exp(-(d * d) / (2 * sig * sig))               # exactly the engine's kernel
    if (AREA) w <- w * pmax(cos(CELLS_LAT * pi / 180), 1e-6)   # ...including its measure
    idx <- sample.int(length(w), 1, prob = w)
    lon[k] <- CELLS_LON[idx]; lat[k] <- CELLS_LAT[idx]
  }
  # The RANKED truth is jittered within its cell; the CHAIN is not.
  #
  # This generator draws cells, so the true position lands exactly on a cell
  # centre, while fit_and_rank() ranks it against posterior draws jittered by
  # U(+/- CELL/2). Whenever the posterior concentrates on the truth's own cell the
  # rank is then DETERMINISTICALLY 0.5, and the rank ECDF acquires a step at 0.5
  # with a deficit in both tails -- which is exactly the shape this arm has always
  # produced, and it is a property of the harness, not of the engine. It hits
  # latitude and not longitude because at -50 deg a 1 deg longitude cell is 71 km
  # against sigma = 80, so longitude spreads over several cells while latitude
  # (111 km) can sit inside one. The sphere generator never showed it because
  # sphere_step() returns a continuous position.
  #
  # The model has no sub-cell structure, so the correct target is the position
  # uniform within its cell. Jitter only the ranked mid-knot: the chain, the
  # endpoints and the emission stay exactly on cell centres, so the movement model
  # remains bit-exactly the engine's.
  list(lon = lon, lat = lat,
       mid_lat = lat[KMID] + runif(1, -CELL/2, CELL/2),
       mid_lon = lon[KMID] + runif(1, -CELL/2, CELL/2))
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
    step_hours = STEP_H, diffusion = DIFFUSION,
    area_correction = AREA)
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

# Two generators with different roles. NOTE that AREA = TRUE changes what each
# arm means, and the previously recorded outcomes were obtained at AREA = FALSE:
#   * fit_kernel gen draws from the engine's OWN movement model, whatever that
#     currently is. Standard SBC (data from the fitted model); it isolates the
#     IMPLEMENTATION, and it cannot detect a mismatch between the model and the
#     sphere, because the generator inherits the same mismatch. This is precisely
#     why SBC did not catch the missing cell-area factor for as long as it did.
#   * sphere gen draws isotropic Brownian motion on the sphere, i.e. the physical
#     movement. This is the arm that tests CORRESPONDENCE rather than
#     self-consistency, and it is the one the cell-area factor should help: at
#     AREA = FALSE it was a misspecification probe expected to miscalibrate
#     latitude; at AREA = TRUE the fit carries the sphere's measure, so agreement
#     should IMPROVE. If it now calibrates cleanly, the sphere arm becomes the
#     substantive claim and fit_kernel becomes the control -- a stronger result
#     than the old framing, not a repair of it.
#     (The transition still omits the sphere's tan(lat) drift, which the area
#     factor does not supply, so residual latitude structure is expected.)
# Which arms to run: "sphere", "fitkernel", or both. Set SBC_ARMS in the
# environment to re-run one arm without paying for the other.
ARMS <- strsplit(Sys.getenv("SBC_ARMS", "sphere,fitkernel"), ",")[[1]]

if ("sphere" %in% ARMS) {
  cat("\n=== robustness probe (isotropic spherical movement) ===\n")
  sbc_sphere <- run_sbc(draw_truth_sphere)
  sbc_ecdf_report(sbc_sphere, SPECS, SBC_L, conf = 0.95,
                  fig = "sbc_sphere.png", tag = "grid HMM | sphere gen")
}

if ("fitkernel" %in% ARMS) {
  cat("\n=== self-consistency SBC (fit's own discrete kernel) ===\n")
  sbc_fit <- run_sbc(draw_truth_fitkernel)
  sbc_ecdf_report(sbc_fit, SPECS, SBC_L, conf = 0.95,
                  fig = "sbc_fitkernel.png", tag = "grid HMM | fit-kernel gen")
}
