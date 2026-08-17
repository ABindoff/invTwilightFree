# GATE before touching the engine: does a CENSORED emission remove the latitude tilt?
#
# THE DIAGNOSIS. `spike_normaliser` (lib.rs:65) is the integral of the raw spike over
# [0, max_light], and the density is divided by it. That is a TRUNCATED likelihood:
# it asserts values outside the sensor range never occur. A light sensor does not
# truncate, it SATURATES -- light above the ceiling is RECORDED AS the ceiling. That
# is CENSORING, and it is a different likelihood.
#
#   truncated  p(y|mu) = raw(y|mu) / Z(mu),  Z(mu) = int_0^ml raw   <- mu-DEPENDENT
#   censored   p(y|mu) = f(y|mu) for 0<y<ml, with constant normaliser
#              C = lam_lo*lam_hi/(lam_lo+lam_hi), plus atoms
#              P(Y<=0)  = (lam_hi/(lam_lo+lam_hi)) exp(-lam_lo*mu)
#              P(Y>=ml) = (lam_lo/(lam_lo+lam_hi)) exp(-lam_hi*(ml-mu))
#              total mass 1 for EVERY mu  <- no mu-dependent normaliser
#
# Under truncation, every interior observation contributes -log Z(mu(cell)) -- a
# DATA-INDEPENDENT term that varies with the cell's latitude through the solar
# geometry. That is the measured tilt. Under censoring interior observations carry
# no such term; the boundary atoms do depend on mu, but that is genuine information
# (a saturated reading really does say the light was at least ml).
#
# So the tilt is NOT an unavoidable property of an asymmetric emission. It is the
# signature of modelling a censored sensor as a truncated one.
#
# PREDICTION, WRITTEN BEFORE RUNNING. Per-knot emission tilt at truth on record
# 2021033/90 is -0.116 nats/deg under the current truncated likelihood. Under
# censoring it should fall to well under 0.02 in magnitude. If it does not, this
# diagnosis is wrong and NO ENGINE CHANGE SHOULD BE MADE.
#
# The R twin is validated against the compiled engine first. If the twin does not
# reproduce `eval_logpk_grid` to numerical precision, nothing below is readable.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
suppressMessages(devtools::load_all(".", quiet = TRUE))
BAT <- readRDS("scratch/nes_calibration/light_battery.rds")
b <- Filter(function(z) z$id == "2021033" && z$offset == 90, BAT)[[1]]

ML    <- b$max_light
LAM   <- 1 / (ML * 0.5)
RATIO <- 2
PSLAB <- 0.10
LAM_HI <- LAM * RATIO

expected_light <- function(z) {
  r <- b$response
  pmin(pmax(r[1] + r[2] / (1 + exp((z - r[3]) / r[4])), 0), ML)
}

# ---- R twin of the CURRENT (truncated) emission -----------------------------
logden_trunc <- function(y, mu) {
  raw <- ifelse(y <= mu, LAM * exp(-LAM * (mu - y)), LAM * exp(-LAM_HI * (y - mu)))
  lo <- (1 - exp(-LAM * mu)) / LAM
  hi <- (1 - exp(-LAM_HI * (ML - mu))) / LAM_HI
  Z  <- pmax((lo + hi) * LAM, 1e-12)
  log((1 - PSLAB) * (raw / Z) + PSLAB / ML)
}

# ---- CENSORED emission ------------------------------------------------------
CC <- LAM * LAM_HI / (LAM + LAM_HI)
logden_cens <- function(y, mu) {
  dens <- ifelse(y <= mu, CC * exp(-LAM * (mu - y)), CC * exp(-LAM_HI * (y - mu)))
  atom_lo <- (LAM_HI / (LAM + LAM_HI)) * exp(-LAM * mu)
  atom_hi <- (LAM  / (LAM + LAM_HI)) * exp(-LAM_HI * (ML - mu))
  # saturated / floored readings are atoms; interior readings are densities
  spike <- ifelse(y >= ML - 1e-9, atom_hi,
           ifelse(y <= 1e-9,      atom_lo, dens))
  log((1 - PSLAB) * spike + PSLAB / ML)
}

# ---- validate the twin against the compiled engine --------------------------
# eval_logpk_grid takes lon/lat as CANDIDATE CELLS and returns one summed
# log-likelihood per cell, over all supplied observations.
set.seed(1)
n <- 4000
ii <- sort(sample(seq_along(b$time), n))
tt <- as.numeric(b$time)[ii]; yy <- b$perfect[ii]
cells_lon <- c(220, 225, 230); cells_lat <- c(40, 45, 50)
eng <- eval_logpk_grid(cells_lon, cells_lat, tt, yy, b$response,
                       c(LAM, ML, PSLAB), RATIO)
twin <- vapply(seq_along(cells_lon), function(c) {
  z <- solar_zenith(tt, rep(cells_lon[c], n), rep(cells_lat[c], n))
  sum(logden_trunc(yy, expected_light(z)))
}, 0)
for (c in seq_along(eng))
  cat(sprintf("  cell (%.0f, %.0f): engine %.6f | twin %.6f | diff %.3e\n",
              cells_lon[c], cells_lat[c], eng[c], twin[c], abs(eng[c] - twin[c])))
if (max(abs(eng - twin)) > 1e-6 * max(1, max(abs(eng)))) {
  cat("*** TWIN DOES NOT MATCH THE ENGINE. Nothing below is readable. STOP.\n")
  quit(status = 1)
}
cat("twin validated.\n\n")

# ---- per-knot latitude tilt, both likelihoods --------------------------------
# knots of 12 h, longitude held at truth (a latitude tilt cannot be produced by
# longitude freedom, and longitude is unbiased in the fits)
tn <- as.numeric(b$time)
kb <- seq(min(tn), max(tn), by = 12 * 3600)
dl <- 0.5                                   # +/- half a degree for the gradient
tilt <- function(fn) {
  out <- numeric(0)
  for (k in seq_len(length(kb) - 1)) {
    sel <- tn > kb[k] & tn <= kb[k + 1]
    if (sum(sel) < 6) next
    t_k <- tn[sel]; y_k <- b$perfect[sel]
    lat0 <- approx(tn, b$lat, mean(t_k), rule = 2)$y
    lon0 <- approx(tn, b$lon, mean(t_k), rule = 2)$y
    up <- sum(fn(y_k, expected_light(solar_zenith(t_k, rep(lon0, sum(sel)),
                                                  rep(lat0 + dl, sum(sel))))))
    dn <- sum(fn(y_k, expected_light(solar_zenith(t_k, rep(lon0, sum(sel)),
                                                  rep(lat0 - dl, sum(sel))))))
    out <- c(out, (up - dn) / (2 * dl))
  }
  out
}

tt_trunc <- tilt(logden_trunc)
tt_cens  <- tilt(logden_cens)

cat(sprintf("knots evaluated: %d\n\n", length(tt_trunc)))
cat("=== per-knot latitude tilt at truth (nats/deg) ===\n")
cat(sprintf("  TRUNCATED (current) : mean %+0.4f   median %+0.4f\n",
            mean(tt_trunc), median(tt_trunc)))
cat(sprintf("  CENSORED  (proposed): mean %+0.4f   median %+0.4f\n",
            mean(tt_cens), median(tt_cens)))
cat(sprintf("\n  reduction in |mean tilt|: %.1f%%\n",
            100 * (1 - abs(mean(tt_cens)) / abs(mean(tt_trunc)))))

cat("\n=== how many observations are actually AT the sensor bounds? ===\n")
cat(sprintf("  at max_light: %.1f%%   at zero: %.1f%%   interior: %.1f%%\n",
            100 * mean(b$perfect >= ML - 1e-9), 100 * mean(b$perfect <= 1e-9),
            100 * mean(b$perfect > 1e-9 & b$perfect < ML - 1e-9)))

cat("\nVERDICT\n")
if (abs(mean(tt_cens)) < 0.02) {
  cat("  GATE PASSED: censoring removes the tilt. Implementing it in the engine is\n")
  cat("  justified, and it is the physically correct model for a saturating sensor.\n")
} else {
  cat("  *** GATE FAILED: censoring does NOT remove the tilt. The diagnosis is wrong.\n")
  cat("  DO NOT change the engine on this basis.\n")
}
