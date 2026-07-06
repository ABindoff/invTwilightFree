# Spherical-metric SBC. Generate a hierarchical panel from the GREAT-CIRCLE
# movement model (latitude-varying Gaussian steps), at high latitude where the
# flat fixed-reference metric is most wrong, and fit both ways via
# TwilightFreeHier. Self-consistency: the spherical fit should calibrate; the
# flat fit should be miscalibrated on wide-latitude tracks (the motivation).
suppressPackageStartupMessages(devtools::load_all(
  "C:/Users/bindoffa/antigravity_projects/invtwilightfree", quiet = TRUE))
source("C:/Users/bindoffa/antigravity_projects/invtwilightfree/inst/sbc/ecdf_bands.R")

DEPLOY_LAT <- -62; DEPLOY_LON <- 150     # high latitude: cos(lat) varies fast
DATE <- "2024-06-15"; STEP_H <- 12; NDAYS <- 12; N_IND <- 3
CAL <- c(558.5, 5.818); LP <- c(0.03, 64, 0.10)
# large movement so tracks span a wide latitude range (where the metric matters):
# sig2 ~ InvGamma(5, beta), beta ~ Gamma(5, 3.1e-5) -> sig2 mean ~40000 km^2 (~200 km/knot)
A_POP <- 5; G0 <- 5; H0 <- 3.1e-5

# knot grid (uniform STEP_H) and matching 10-min obs times
t0 <- as.numeric(as.POSIXct(DATE, tz = "UTC")); hstep <- STEP_H * 3600
K <- 2L * NDAYS + 1L
KNOTS <- t0 + (0:(K - 1)) * hstep
ut <- seq(t0, KNOTS[K], by = 600)
times <- as.POSIXct(ut, origin = "1970-01-01", tz = "UTC")
kidx <- pmin(pmax(findInterval(ut, KNOTS, left.open = TRUE) + 1L, 1L), K)

sample_spike <- function(mu) {
  lam <- LP[1]; maxl <- LP[2]
  m_lo <- 1 - exp(-lam * mu); m_hi <- 0.5 * (1 - exp(-2 * lam * (maxl - mu)))
  if (runif(1) < m_lo / (m_lo + m_hi)) mu - (-log(1 - runif(1) * m_lo) / lam)
  else                                 mu + (-log(1 - runif(1) * (2 * m_hi)) / (2 * lam))
}
# great-circle step: isotropic Gaussian in km on the local tangent plane
gen_track_light <- function(sig2) {
  s <- sqrt(sig2)
  lon <- numeric(K); lat <- numeric(K); lon[1] <- DEPLOY_LON; lat[1] <- DEPLOY_LAT
  for (k in 2:K) {
    dy <- rnorm(1, 0, s); dx <- rnorm(1, 0, s)
    lat[k] <- lat[k-1] + dy / 111.0
    lon[k] <- lon[k-1] + dx / (111.0 * cos(lat[k-1] * pi / 180))
  }
  z <- solar_zenith(ut, lon[kidx], lat[kidx]); mu <- pmin(pmax(CAL[1] - CAL[2] * z, 0), LP[2])
  light <- vapply(seq_along(ut), function(j)
    if (runif(1) < LP[3]) runif(1, 0, LP[2]) else sample_spike(mu[j]), numeric(1))
  list(df = data.frame(Date = times, Light = light), lat_end = lat[K], lat_rng = diff(range(lat)))
}

REPS <- if (exists("REPS_S")) REPS_S else 60L
common <- list(step_hours = STEP_H, calibration = CAL, likelihood_params = LP,
               a_pop = A_POP, hyperprior = c(G0, H0), mesh_pad = 32, coarse_res = 2.5,
               surrogate_diffusion = 250, sweeps = 3000L, burn = 1200L, thin = 4L)
rank_sph <- integer(0); rank_flat <- integer(0); rng_lat <- numeric(0); tstart <- Sys.time()
for (r in seq_len(REPS)) {
  beta <- rgamma(1, G0, H0); sig2 <- 1 / rgamma(N_IND, A_POP, beta)
  gg <- lapply(sig2, gen_track_light)
  rng_lat <- c(rng_lat, mean(sapply(gg, function(g) g$lat_rng)))
  ids <- sprintf("t%d", seq_len(N_IND))
  data <- setNames(lapply(gg, function(g) g$df), ids)
  locs <- data.frame(id = ids, deploy_lon = DEPLOY_LON, deploy_lat = DEPLOY_LAT,
                     retrieve_lon = NA, retrieve_lat = NA)
  fit_s <- do.call(TwilightFreeHier, c(list(data = data, locations = locs, metric = "spherical", seed = r), common))
  fit_f <- do.call(TwilightFreeHier, c(list(data = data, locations = locs, metric = "flat", seed = r), common))
  dt_day <- STEP_H / 24
  for (i in seq_len(N_IND)) {
    true_kmday <- sqrt(sig2[i]) / sqrt(dt_day)
    rank_sph  <- c(rank_sph,  sum(fit_s$draws$sig2_kmday[, i] < true_kmday))
    rank_flat <- c(rank_flat, sum(fit_f$draws$sig2_kmday[, i] < true_kmday))
  }
  if (r %% 5L == 0L) cat(sprintf("  SBC-sph %d/%d (%.0fs) mean lat range=%.1f deg\n",
                                 r, REPS, as.numeric(Sys.time() - tstart, units = "secs"), mean(rng_lat)))
}
L <- length(fit_s$draws$sig2_kmday[, 1])
cat(sprintf("\nSpherical-metric SBC: %d reps, N=%d, L=%d ; mean latitude range %.1f deg\n",
            REPS, N_IND, L, mean(rng_lat)))
sbc_ecdf_report(data.frame(rank = rank_sph), c(rank = "sigma_i (spherical fit)"),
                L, conf = 0.95, fig = "sbc_spherical.png", tag = "spherical fit")
sbc_ecdf_report(data.frame(rank = rank_flat), c(rank = "sigma_i (flat fit)"),
                L, conf = 0.95, fig = "sbc_flat.png", tag = "flat fit")
