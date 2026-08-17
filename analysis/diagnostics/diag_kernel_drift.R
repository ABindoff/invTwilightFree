# Does the movement kernel itself drift equatorward, and how fast?
#
# WHY THIS, AND WHY IT IS CHEAP. Stage 2 said the area correction is behaving
# correctly (removing it moves estimates NORTH, matching the documented 0.204
# nats/step poleward defect) and is NOT the source of the -0.567 deg. What it also
# showed is that the bias GROWS when the movement prior is widened. This measures
# the prior's one-step behaviour directly. It is pure numerics -- no HMM fits, no
# light, seconds to run -- so it costs essentially nothing and is informative
# whichever way stage 2 lands.
#
# THE ANALYTIC PREDICTION. For a point at geodesic angular distance d and bearing
# theta from (phi0, lambda0):
#     sin(phi) = sin(phi0) cos(d) + cos(phi0) sin(d) cos(theta)
# Averaging over uniform bearing kills the second term, leaving
#     E[sin phi] = sin(phi0) * E[cos d]  <  sin(phi0)
# So an ISOTROPIC kernel on a sphere pulls the mean of sin(latitude) toward zero,
# i.e. toward the EQUATOR, before any likelihood is involved. With E[d^2] = 2 sigma^2
# (two dimensions) in units of Earth radii, E[cos d] ~= 1 - sigma^2/R^2 and
#     drift ~= -tan(phi0) * sigma^2 / R^2   radians per step.
#
# This is CORRECT behaviour for Brownian motion on a sphere -- spherical BM relaxes
# to the uniform distribution, whose E[sin phi] is 0 -- but it is a poor property for
# an ANIMAL movement prior, which should not think a seal prefers the equator. And it
# scales as sigma^2, which is why widening the prior amplified the bias.
#
# Uniform bearing is what the CELL AREA FACTOR buys. Dropping the area factor
# over-weights the small high-latitude cells and pushes back poleward, which is why
# turning the correction off partly cancelled the bias instead of fixing anything.
#
# The kernel below reproduces the engine exactly: haversine distance on R = 6371,
# sigma = diffusion * sqrt(dt), var2 = 2 sigma^2, weight exp(-dist^2/var2), the
# 5-sigma truncation and both coarse pre-filters (lib.rs run_grid_hmm).
#
# The fit-level JUSTIFICATION TEST from stage 2 is evaluated first and reported, so
# the conclusion is gated on it even though the numerics run regardless.
suppressMessages(library(data.table))

R_E <- 6371.0
DT  <- 0.5          # days between knots (12 h)

# ---- justification test from stage 2 ---------------------------------------
gate <- "NOT EVALUATED (area_diffusion.csv absent)"
justified <- NA
f <- "scratch/nes_calibration/area_diffusion.csv"
if (file.exists(f)) {
  A <- fread(f, colClasses = list(character = "id"))
  x <- A[area == TRUE]
  w <- dcast(x, id ~ diff, value.var = "bias_mode")
  if (all(c("110", "440") %in% names(w)) && sum(complete.cases(w)) >= 4) {
    w <- w[complete.cases(w)]
    d <- w[["440"]] - w[["110"]]
    justified <- mean(d) < -0.1 && sum(d < 0) >= ceiling(0.8 * length(d))
    gate <- sprintf(
      "D440 - D110 (area ON) = %+0.3f deg, more negative on %d/%d tracks -> %s",
      mean(d), sum(d < 0), length(d),
      if (isTRUE(justified)) "JUSTIFIED" else "NOT JUSTIFIED")
  } else {
    gate <- sprintf("D=440 arm incomplete (%d tracks with both) -- cannot gate",
                    sum(complete.cases(w)))
  }
}
cat("=== stage 2 justification test ===\n  ", gate, "\n\n", sep = "")

# ---- the engine's kernel, reproduced ---------------------------------------
haversine <- function(lat1, lon1, lat2, lon2) {
  p <- pi / 180
  a <- sin((lat2 - lat1) * p / 2)^2 +
       cos(lat1 * p) * cos(lat2 * p) * sin((lon2 - lon1) * p / 2)^2
  2 * R_E * asin(pmin(1, sqrt(a)))
}

# One-step mean latitude displacement induced by the kernel alone.
step_drift <- function(phi0, D, cell = 1, area = TRUE, pad_sigma = 8) {
  sigma <- D * sqrt(DT)
  var2  <- 2 * sigma^2
  thr   <- 5 * sigma                       # engine's truncation

  # grid centred so phi0 is a cell centre; extent generous relative to sigma
  pad_deg <- max(6, pad_sigma * sigma / 111)
  lats <- seq(phi0 - pad_deg, phi0 + pad_deg, by = cell)
  lons <- seq(-pad_deg / cos(phi0 * pi / 180), pad_deg / cos(phi0 * pi / 180), by = cell)
  g    <- expand.grid(lat = lats, lon = lons)

  d <- haversine(phi0, 0, g$lat, g$lon)
  keep <- d <= thr &
          abs(g$lat - phi0) <= thr / 111 &
          abs(g$lon)        <= thr / (111 * pmax(0.1, pmin(cos(phi0 * pi / 180),
                                                           cos(g$lat * pi / 180))))
  if (!any(keep)) return(NA_real_)
  lat <- g$lat[keep]; dd <- d[keep]
  wgt <- exp(-dd^2 / var2) * (if (area) cos(lat * pi / 180) else 1)
  sum(wgt * lat) / sum(wgt) - phi0
}

# analytic: E[sin phi] = sin(phi0) * E[cos d], E[d^2] = 2 sigma^2
analytic <- function(phi0, D) {
  sigma <- D * sqrt(DT)
  p <- pi / 180
  (asin(sin(phi0 * p) * (1 - (sigma / R_E)^2)) / p) - phi0
}

PHI <- c(35, 40, 45, 50, 55)
DD  <- c(55, 110, 220, 440)

cat("=== one-step mean latitude displacement, degrees (area ON) ===\n")
cat(sprintf("%6s %10s %10s %10s %12s\n",
            "lat", "D", "numeric", "analytic", "x480 knots"))
res <- list()
for (p0 in PHI) for (dv in DD) {
  num <- step_drift(p0, dv, area = TRUE)
  ana <- analytic(p0, dv)
  cat(sprintf("%6.0f %10.0f %10.5f %10.5f %12.2f\n", p0, dv, num, ana, num * 480))
  res[[length(res) + 1]] <- data.table(lat = p0, D = dv, numeric = num,
                                       analytic = ana)
}
R1 <- rbindlist(res)
cat(sprintf("\n  numeric vs analytic: max |diff| = %.6f deg, correlation %.5f\n",
            max(abs(R1$numeric - R1$analytic)), cor(R1$numeric, R1$analytic)))
cat(sprintf("  drift scales as sigma^2? ratio D440/D110 = %.2f (sigma^2 ratio = 16)\n",
            mean(R1[D == 440]$numeric / R1[D == 110]$numeric)))

cat("\n=== area factor OFF, same kernel ===\n")
cat(sprintf("%6s %10s %10s %10s\n", "lat", "D", "area ON", "area OFF"))
for (p0 in PHI) for (dv in c(110, 440)) {
  on  <- step_drift(p0, dv, area = TRUE)
  off <- step_drift(p0, dv, area = FALSE)
  cat(sprintf("%6.0f %10.0f %10.5f %10.5f\n", p0, dv, on, off))
}

cat("\n=== grid resolution: is the discrete sum tracking the continuum? ===\n")
cat(sprintf("%6s %10s %10s %10s\n", "lat", "cell 1.0", "cell 0.5", "analytic"))
for (p0 in PHI) {
  cat(sprintf("%6.0f %10.5f %10.5f %10.5f\n", p0,
              step_drift(p0, 110, cell = 1.0), step_drift(p0, 110, cell = 0.5),
              analytic(p0, 110)))
}

cat("\nREADING\n")
cat("  If numeric ~= analytic, the engine is faithfully implementing spherical\n")
cat("  Brownian motion and the equatorward pull is a PROPERTY OF THE PRIOR, not a\n")
cat("  coding error. The fix is then a modelling one -- a prior whose stationary\n")
cat("  distribution is not uniform-on-the-sphere -- not a patch to the kernel.\n")
cat("  If numeric and analytic DIVERGE, the discretisation is at fault and a finer\n")
cat("  grid should close the gap.\n")
if (isTRUE(justified)) {
  cat("\n  Stage 2 JUSTIFIED this test: the fit-level bias grew with prior width.\n")
} else if (isFALSE(justified)) {
  cat("\n  NOTE: stage 2 did NOT justify this test at fit level. The numbers above\n")
  cat("  still describe the kernel, but do not claim they explain the fitted bias.\n")
}
