# IS THE LATITUDE BIAS A GAUGE FREEDOM? Fisher-information eigenanalysis.
#
# THE HYPOTHESIS (not mine -- posed as a reframing of a day of null results). Weeks
# of axis-aligned hyperparameter sweeps returning null, or resolving bias in one
# place and having it reappear in another, is the signature of a SLOPPY MODEL in
# Sethna's sense: the likelihood has a near-flat direction, a combination such as
# (attenuation x threshold x latitude) that the within-track data cannot see. Every
# one-at-a-time sweep is then a coordinate move that the fit re-projects onto the
# invariant manifold, so each returns null individually. The culprit is not a
# parameter, it is a DIRECTION.
#
# ONE REFINEMENT, which this script is built to test rather than assume. A flat
# direction by itself produces VARIANCE, not BIAS -- it is an amplifier, not a force.
# Measured today: on data drawn from the model the latitude bias is -0.51 deg with
# coverage 0.996; on real data it is -1.54 with coverage 0.64. The gauge freedom is a
# property of the likelihood and is present in both. So the picture that fits
# everything is MISSPECIFICATION SUPPLIES THE PUSH AND THE SOFT DIRECTION SUPPLIES
# THE LEVERAGE. That predicts a soft direction with real latitude weight, which is
# exactly what is measured below.
#
# THE MEASUREMENT. For a one-day window at a known position, the emission
# log-likelihood is a function of
#     theta = (latitude, z50, scale, amp)
# where z50 is the assumed solar-elevation threshold, scale the twilight width, and
# amp the sensor attenuation/gain. Compute the observed information -d2 logL/dtheta2
# at the truth, eigendecompose, and read off:
#   * the eigenvalue SPECTRUM -- sloppy models span many orders of magnitude;
#   * the SOFTEST eigenvector -- if it mixes latitude with z50/amp, latitude is not
#     separately identified from the calibration within a single day;
#   * how the softest eigenvector's latitude component behaves ACROSS SEASON.
#
# PRE-REGISTERED PREDICTIONS:
#   1. Eigenvalue spectrum spans >= 3 orders of magnitude (sloppy).
#   2. The softest eigenvector has substantial latitude weight (|cos| > 0.3 with the
#      latitude axis), i.e. latitude is one of the soft coordinates.
#   3. The latitude weight of the soft direction is LARGEST near the equinox, where
#      the day-length-to-latitude map degenerates.
#   4. THE SHARP ONE: the ratio of the latitude component to the z50 component in the
#      soft eigenvector CHANGES SIGN across the year. A single gauge freedom would
#      then explain the observed seasonal sign flip in latitude bias without needing
#      a separate mechanism for each half-year. If it does not flip, the gauge story
#      does not explain the seasonal structure and needs amending.
#
# Uses the R twin of the engine's emission, validated to 1e-10 against
# eval_logpk_grid. Longitude held at truth: it is well identified (bias +0.2 deg) and
# a latitude gauge cannot be produced by longitude freedom.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
suppressMessages(devtools::load_all(".", quiet = TRUE))
BAT <- readRDS("scratch/nes_calibration/light_battery.rds")
PSLAB <- 0.10; RATIO <- 2

# emission log-likelihood for a window, as a function of theta
mk_logL <- function(tt, y, lon0, ML, lam) {
  lam_hi <- lam * RATIO
  function(th) {
    lat <- th[1]; z50 <- th[2]; sc <- th[3]; amp <- th[4] * ML; flo <- th[5] * ML
    z  <- solar_zenith(tt, rep(lon0, length(tt)), rep(lat, length(tt)))
    mu <- pmin(pmax(flo + amp / (1 + exp((z - z50) / sc)), 0), ML)
    raw <- ifelse(y <= mu, lam * exp(-lam * (mu - y)), lam * exp(-lam_hi * (y - mu)))
    lo <- (1 - exp(-lam * mu)) / lam
    hi <- (1 - exp(-lam_hi * (ML - mu))) / lam_hi
    sum(log((1 - PSLAB) * (raw / pmax((lo + hi) * lam, 1e-12)) + PSLAB / ML))
  }
}

hess <- function(f, th, h) {
  p <- length(th); H <- matrix(0, p, p)
  for (i in seq_len(p)) for (j in i:p) {
    ei <- ej <- numeric(p); ei[i] <- h[i]; ej[j] <- h[j]
    H[i, j] <- H[j, i] <-
      (f(th + ei + ej) - f(th + ei - ej) - f(th - ei + ej) + f(th - ei - ej)) /
      (4 * h[i] * h[j])
  }
  H
}

PNAME <- c("lat", "z50", "scale", "amp", "floor")
res <- list()
for (b in BAT) {
  ML <- b$max_light; lam <- 1 / (ML * 0.5); r <- b$response
  tn <- as.numeric(b$time)
  # one-day windows, sampled every 10 days
  starts <- seq(min(tn), max(tn) - 86400, by = 10 * 86400)
  for (s in starts) {
    sel <- tn >= s & tn < s + 86400
    if (sum(sel) < 30) next
    tt <- tn[sel]; y <- b$perfect[sel]
    lat0 <- approx(tn, b$lat, mean(tt), rule = 2)$y
    lon0 <- approx(tn, lon360(b$lon), mean(tt), rule = 2)$y
    if (!approx(tn, as.numeric(b$supported), mean(tt), rule = 2)$y > 0.999) next
    f  <- mk_logL(tt, y, lon0, ML, lam)
    th <- c(lat0, r[3], r[4], r[2] / ML, r[1] / ML)
    h  <- c(0.25, 0.25, 0.25, 0.01, 0.01)
    H  <- -hess(f, th, h)                       # observed information
    ev <- eigen((H + t(H)) / 2, symmetric = TRUE)
    val <- ev$values
    if (min(val) <= 0 || any(!is.finite(val))) next
    v_soft <- ev$vectors[, which.min(val)]
    res[[length(res) + 1]] <- data.table(
      id = b$id, offset = b$offset,
      decl = solar_declination(as.POSIXct(mean(tt), origin = "1970-01-01", tz = "UTC")),
      lat = lat0,
      cond = max(val) / min(val),
      lat_wt = abs(v_soft[1]),
      lat_z50 = v_soft[1] / v_soft[2],
      lat_amp = v_soft[1] / v_soft[4],
      stiff_lat = abs(ev$vectors[1, which.max(val)]))
  }
}
D <- rbindlist(res)
cat(sprintf("%d one-day windows, %d records\n\n", nrow(D), uniqueN(paste(D$id, D$offset))))

cat("=== PREDICTION 1: is the spectrum sloppy? ===\n")
cat(sprintf("  condition number (max/min eigenvalue): median %.3g, range %.3g - %.3g\n",
            median(D$cond), min(D$cond), max(D$cond)))
cat(sprintf("  orders of magnitude spanned: median %.1f\n", median(log10(D$cond))))
cat(sprintf("  => %s\n", if (median(log10(D$cond)) >= 3) "SLOPPY, as predicted" else
                          "NOT sloppy by the >=3 orders criterion"))

cat("\n=== PREDICTION 2: does the soft direction involve latitude? ===\n")
cat(sprintf("  |latitude component| of softest eigenvector: median %.3f\n",
            median(D$lat_wt)))
cat(sprintf("  |latitude component| of STIFFEST eigenvector: median %.3f\n",
            median(D$stiff_lat)))
cat(sprintf("  => %s\n", if (median(D$lat_wt) > 0.3)
      "latitude IS a soft coordinate, as predicted" else
      "latitude is NOT strongly in the soft direction"))

cat("\n=== PREDICTIONS 3 & 4: season ===\n")
D[, sband := cut(decl, c(-24, -15, -7, 7, 15, 24), include.lowest = TRUE,
                 labels = c("NH winter", "-15..-7", "equinox", "+7..+15", "NH summer"))]
print(as.data.frame(D[, .(n = .N,
                          cond = signif(median(cond), 3),
                          lat_wt = round(median(lat_wt), 3),
                          lat_over_z50 = round(median(lat_z50), 3),
                          lat_over_amp = round(median(lat_amp), 3)),
                      by = sband][order(sband)]), row.names = FALSE)

sgn <- D[, .(m = median(lat_z50)), by = sband]
cat(sprintf("\n  PREDICTION 3 (lat weight peaks at equinox): %s\n",
            if (D[sband == "equinox", median(lat_wt)] >=
                max(D[sband %in% c("NH winter", "NH summer"), median(lat_wt), by = sband]$V1))
              "HOLDS" else "FAILS"))
cat(sprintf("  PREDICTION 4 (lat/z50 ratio changes sign across the year): %s\n",
            if (any(sgn$m > 0, na.rm = TRUE) && any(sgn$m < 0, na.rm = TRUE))
              "HOLDS -- one gauge freedom predicts the seasonal sign structure" else
              "FAILS -- the sign structure needs a separate explanation"))
