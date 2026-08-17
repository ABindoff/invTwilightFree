# THREE PRE-REGISTERED TESTS of the emission-tilt mechanism.
#
# The mechanism (from independent code analysis, not from fitting these numbers):
# the per-observation emission carries a latitude tilt pointing toward the
# LONGER-DAY side, from two components acting on the daytime shoulder of the
# response -- the spike normaliser's gradient (Z shrinks as expected light nears
# the max_light clamp) and the arm asymmetry (lam_hi = 2*lam_lo, so a noiseless
# daytime shoulder charges poleward candidates at 2*lambda and equatorward at
# lambda). Measured per-knot tilt at truth: full -0.116 nats/deg, normaliser
# frozen -0.063, symmetric arms -0.071, BOTH removed -0.004.
#
# PREDICTIONS, WRITTEN BEFORE RUNNING. These are the point of the exercise: each
# is a number this mechanism requires and the alternatives do not.
#
#   1. shade_ratio = 1  ->  bias goes MORE NEGATIVE, about -1.2.
#      This is the counterintuitive one. Any "the asymmetry is the culprit, so
#      symmetrise it" story predicts the bias should IMPROVE. It should not:
#      removing the asymmetry widens the posterior, and the normaliser tilt is
#      then amplified by the larger sum-of-covariances.
#      Out-of-sample support already on disk: the 2026-08-06 shade_ratio sweep on
#      REAL data gave ratio 1.00 -> -1.88, 2.00 -> -0.88, 4.00 -> +1.10.
#
#   2. lam_mult = 2  ->  bias about -0.4, i.e. BETTER.
#      Tilt doubles but sigma^2 quarters, so the product falls.
#
#   3. SOUTHERN HEMISPHERE MIRROR  ->  bias flips to about +0.6.
#      The strongest test. Light is regenerated from the same response at
#      (-lat, lon) on the same dates, so the track is a true reflection. The bias
#      is EQUATORWARD, and equatorward in the southern hemisphere means INCREASING
#      latitude, hence a sign flip in (fitted - true). No apparatus artefact --
#      grid indexing, interpolation, binning -- predicts a hemisphere flip.
#
# Primary record 2021033 offset 90 (the one profiled). 2023032 offset 0 repeated
# as a cross-record check. Baseline arm is the reproduction control.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
suppressMessages(devtools::load_all(".", quiet = TRUE))
BAT <- readRDS("scratch/nes_calibration/light_battery.rds")
OUT <- "scratch/nes_calibration/prereg.csv"
STEP_H <- 12; CELL <- 1; PSLAB <- 0.10; DIFF <- 110

pick <- function(id, off) Filter(function(z) z$id == id && z$offset == off, BAT)[[1]]
RECS <- list(pick("2021033", 90), pick("2023032", 0))

fit_arm <- function(b, arm) {
  r <- b$response; ml <- b$max_light
  lat <- b$lat; lon <- b$lon; y <- b$perfect
  lam_mult <- 1; sratio <- 2

  if (arm == "shade1")   sratio  <- 1
  if (arm == "lam2")     lam_mult <- 2
  if (arm == "mirror") {
    lat <- -b$lat
    z <- solar_zenith(as.numeric(b$time), lon, lat)
    y <- pmin(pmax(r[1] + r[2] / (1 + exp((z - r[3]) / r[4])), 0), ml)
  }

  gg <- rast(xmin = floor(min(lon)) - 12, xmax = ceiling(max(lon)) + 12,
             ymin = max(-88, floor(min(lat)) - 12),
             ymax = min(88, ceiling(max(lat)) + 12),
             resolution = CELL, crs = "EPSG:4326")
  values(gg) <- 1
  lam <- lam_mult / (ml * 0.5)
  n <- length(lat)
  invisible(capture.output(
    f <- TwilightFreeGrid(b$time, y, grid = gg,
      start_lon = lon[1], start_lat = lat[1],
      end_lon = lon[n], end_lat = lat[n],
      step_hours = STEP_H, diffusion = DIFF, calibration = r,
      likelihood_params = c(lam, ml, PSLAB), shade_ratio = sratio,
      area_correction = TRUE)))

  gp <- grid_posterior(f); K <- nrow(gp$P)
  mu <- sd <- rep(NA_real_, K)
  for (k in seq_len(K)) {
    w <- gp$P[k, ]; s <- sum(w)
    if (!is.finite(s) || s <= 0) next
    w <- w / s
    mu[k] <- sum(w * gp$lat); sd[k] <- sqrt(sum(w * (gp$lat - mu[k])^2))
  }
  tk <- as.numeric(f$fit$time); tt <- as.numeric(b$time)
  tlat <- approx(tt, lat, tk, rule = 2)$y
  sup  <- approx(tt, as.numeric(b$supported), tk, rule = 2)$y > 0.999
  ok <- sup & is.finite(mu)
  data.table(id = b$id, offset = b$offset, arm = arm, n = sum(ok),
             bias_mean = mean((mu - tlat)[ok]),
             bias_mode = mean((f$fit$lat - tlat)[ok]),
             post_sd = mean(sd[ok], na.rm = TRUE))
}

ARMS <- c("baseline", "shade1", "lam2", "mirror")
done <- if (file.exists(OUT)) {
  d0 <- fread(OUT, colClasses = list(character = "id"))
  paste(d0$id, d0$offset, d0$arm)
} else character(0)

for (b in RECS) for (a in ARMS) {
  key <- paste(b$id, b$offset, a)
  if (key %in% done) next
  r <- try(fit_arm(b, a), silent = TRUE)
  if (inherits(r, "try-error")) {
    message("FAILED ", key, ": ", conditionMessage(attr(r, "condition"))); next
  }
  fwrite(r, OUT, append = file.exists(OUT))
  message(key, " done")
}

D <- fread(OUT, colClasses = list(character = "id"))
cat("\n=== results ===\n")
print(as.data.frame(D[, .(id, offset, arm, n,
                          bias_mean = round(bias_mean, 3),
                          bias_mode = round(bias_mode, 3),
                          post_sd = round(post_sd, 3))]), row.names = FALSE)

cat("\n=== scorecard (primary record 2021033/90, posterior mean) ===\n")
P <- D[id == "2021033" & offset == 90]
g <- function(a) if (nrow(P[arm == a])) P[arm == a]$bias_mean else NA_real_
pred <- c(shade1 = -1.2, lam2 = -0.4, mirror = +0.6)
cat(sprintf("  baseline           : %+0.3f\n", g("baseline")))
for (a in names(pred)) {
  obs <- g(a)
  cat(sprintf("  %-18s : %+0.3f   predicted %+0.2f   %s\n", a, obs, pred[a],
              if (is.na(obs)) "-" else
              if (a == "mirror") { if (obs > 0) "SIGN FLIP as predicted" else "*** NO FLIP" }
              else if (a == "shade1") { if (obs < g("baseline")) "worse, as predicted" else "*** IMPROVED" }
              else { if (obs > g("baseline")) "better, as predicted" else "*** WORSE" }))
}
cat("\n  the shade1 arm is the discriminating one: a symmetrise-the-arms story\n")
cat("  predicts improvement, this mechanism predicts deterioration.\n")
