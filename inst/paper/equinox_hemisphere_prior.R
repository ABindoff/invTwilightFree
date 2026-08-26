# =============================================================================
# Equinox hemisphere prior demonstration (manuscript Figure: equinox).
#
# Near an equinox and at low latitude, day length is ~12 h in BOTH hemispheres,
# so the light likelihood for latitude is near-symmetric about the equator:
# the reconstruction is genuinely bimodal and can flip hemisphere. This script
# shows that a single additive `hemisphere_prior()` term resolves the ambiguity
# at negligible cost, WITHOUT distorting a case where light already determines
# latitude.
#
# Fully reproducible: fixed seed, data simulated from the package's own solar
# model + normalised spike-and-slab emission.
#
# Run from the package root:  Rscript inst/paper/equinox_hemisphere_prior.R
# =============================================================================

suppressPackageStartupMessages(library(invTwilightFree))
set.seed(20240320)

# ---- ground truth: a near-stationary animal at LOW latitude (+5 deg) within a
#      few days of the exact March equinox (2024-03-20), with FREE endpoints.
#      The twilight-free likelihood normally resolves latitude even near the
#      equinox from the twilight SLOPE; the sign becomes genuinely ambiguous
#      only in this knife-edge regime (low |lat|, tight window about equinox),
#      which is exactly where a hemisphere prior earns its keep. -----------------
TRUE_LAT <- 5.0
TRUE_LON <- 150.0
DEPLOY   <- c(lon = TRUE_LON, lat = TRUE_LAT)
CAL      <- c(64, 64 / 90)          # intercept, slope (light = clamp(a - b*zenith))
LP       <- c(0.5, 64, 0.05)        # lambda, max_light, prob_slab
STEP_H   <- 24
NDAYS    <- 5                       # tight window straddling the equinox
CELL     <- 1.0

t0    <- as.numeric(as.POSIXct("2024-03-17", tz = "UTC"))
times <- seq(t0, t0 + NDAYS * 86400, by = 600)   # 10-min cadence

# normalised spike sampler (matches inst/sbc/sbc_grid_hmm.R)
sample_spike <- function(mu) {
  lam <- LP[1]; maxl <- LP[2]
  m_lo <- 1 - exp(-lam * mu)
  m_hi <- 0.5 * (1 - exp(-2 * lam * (maxl - mu)))
  if (runif(1) < m_lo / (m_lo + m_hi)) mu + log(1 - runif(1) * m_lo) / lam
  else mu - log(1 - runif(1) * (2 * m_hi)) / (2 * lam)
}
z     <- solar_zenith(times, rep(TRUE_LON, length(times)), rep(TRUE_LAT, length(times)))
mu    <- pmin(pmax(CAL[1] - CAL[2] * z, 0), LP[2])
light <- vapply(seq_along(times),
                function(j) if (runif(1) < LP[3]) runif(1, 0, LP[2]) else sample_spike(mu[j]),
                numeric(1))
dt <- as.POSIXct(times, origin = "1970-01-01", tz = "UTC")

# ---- grid spanning BOTH hemispheres so a flip is possible --------------------
GRID <- makeGrid(lon = TRUE_LON + c(-15, 15), lat = c(-30, 30), cell.size = CELL)

fit_grid <- function(terms = list()) {
  TwilightFreeGrid(dt, light, GRID,
                   # endpoints FREE: no positional anchor, so the hemisphere is
                   # determined by light (and any prior) alone
                   step_hours = STEP_H, diffusion = 40,
                   calibration = CAL, likelihood_params = LP,
                   terms = terms)
}

# marginal latitude posterior at the MIDDLE knot (farthest from both free ends,
# where the ambiguity is strongest)
lat_marginal <- function(fit) {
  post <- grid_posterior(fit)
  k <- round(nrow(post$P) / 2)
  ul <- sort(unique(round(post$lat, 6)))
  p  <- sapply(ul, function(l) sum(post$P[k, abs(post$lat - l) < 1e-6]))
  data.frame(lat = ul, p = p / sum(p))
}

post_mode_lat <- function(m) m$lat[which.max(m$p)]
south_mass    <- function(m) sum(m$p[m$lat < 0])
post_sd_lat   <- function(m) { mu <- sum(m$lat * m$p); sqrt(sum((m$lat - mu)^2 * m$p)) }

# ---- (1) light only ---------------------------------------------------------
fit_light <- fit_grid()
m_light   <- lat_marginal(fit_light)

# ---- (2) light + hemisphere prior ("N" through the window) ------------------
north <- hemisphere_prior(function(d) "N", softness = 1e-3)
fit_prior <- fit_grid(terms = list(
  location_term("hemisphere", source = north, rule = identity_rule())))
m_prior   <- lat_marginal(fit_prior)

cat(sprintf("\nTrue latitude: +%.1f deg (northern)\n", TRUE_LAT))
cat(sprintf("Light only     : mode lat = %+6.1f | SD = %5.1f | P(southern) = %.2f\n",
            post_mode_lat(m_light), post_sd_lat(m_light), south_mass(m_light)))
cat(sprintf("Light + prior  : mode lat = %+6.1f | SD = %5.1f | P(southern) = %.2f\n",
            post_mode_lat(m_prior), post_sd_lat(m_prior), south_mass(m_prior)))

# ---- figure -----------------------------------------------------------------
out_png <- if (exists("OUT_PNG")) OUT_PNG else "equinox_hemisphere_prior.png"
grDevices::png(out_png, width = 1500, height = 650, res = 150)
op <- graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(m_light$lat, m_light$p, type = "h", lwd = 3, col = "grey55",
     xlab = "latitude (deg)", ylab = "posterior probability",
     main = "Light only")
graphics::abline(v = TRUE_LAT, col = "firebrick", lwd = 2, lty = 2)
graphics::abline(v = 0, col = "grey70", lty = 3)
plot(m_prior$lat, m_prior$p, type = "h", lwd = 3, col = "steelblue",
     xlab = "latitude (deg)", ylab = "posterior probability",
     main = "Light + hemisphere prior")
graphics::abline(v = TRUE_LAT, col = "firebrick", lwd = 2, lty = 2)
graphics::abline(v = 0, col = "grey70", lty = 3)
graphics::par(op)
grDevices::dev.off()
cat(sprintf("\nWrote %s\n", out_png))
